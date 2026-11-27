"""

Make a colour image cut out from the Rubin DP2 deep_coadds, optionally with radio (e.g. MeerKAT) contours
overlaid, styled like the output of makeContourOverlay_MKAT.py (no coordinate axes; the object name,
redshift and a scale bar are drawn on the image instead).

Based on Mining-with-MeerKAT_Images_CutOut.ipynb. The coadds are found through the RSP ObsCore table, and
the cutouts come from the RSP image cutout (SODA) service. The cutout can straddle patch boundaries, so for
each band the pieces from every overlapping patch in the best tract are stitched onto that tract's pixel
grid. Stitched cutouts are cached (as FITS, one per band) in CACHE_DIR, so re-styling a plot is quick.

Works on the RSP itself, or anywhere else given an RSP access token (https://data.lsst.cloud/auth/tokens)
saved in ~/.rsp-tap.token.

Run with -h for options. Examples:

    python makeDP2CutOut.py 349.3856494 0.1750374 4 "NGC7603 radio galaxy" --redshift 0.0295
    python makeDP2CutOut.py 349.3856494 0.1750374 4 NGC7603 --bands gri --contours meerkat
    python makeDP2CutOut.py 349.3856494 0.1750374 4 NGC7603 --contours myRadioImage.fits -o NGC7603.jpg

"""

import os
import io
import sys
import argparse
import numpy as np
import astropy.table as atpy
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clipped_stats
from astropy.visualization import make_lupton_rgb
from astropy.utils.data import download_file
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path

CACHE_DIR="DP2CutOutCache"
EDR_URL="https://extragalactic.phys.wits.ac.za/mining-with-meerkat/EDR/data/"   # trailing slash is important
BAND_ORDER="ugrizy"       # Increasing wavelength
# Per-band flux multipliers for colour balance (images are in nJy, so redder bands are brighter).
# g, r, z are the values used in Mining-with-MeerKAT_Images_CutOut.ipynb; the others are interpolated.
BAND_SCALES={'u': 2.5, 'g': 1.6, 'r': 0.9, 'i': 0.7, 'z': 0.5, 'y': 0.45}
FIG_SIZE=8.25             # inches (square), as in makeContourOverlay_MKAT.py
DPI=192

#------------------------------------------------------------------------------------------------------------
def getSession(tokenFile = "~/.rsp-tap.token"):
    """Return an authenticated session for the RSP services.

    """

    try:
        # Works when running on the RSP itself
        from lsst.rsp.utils import get_pyvo_auth
        return get_pyvo_auth()
    except ImportError:
        import requests
        with open(os.path.expanduser(tokenFile)) as inFile:
            token=inFile.read().strip()
        session=requests.Session()
        session.headers["Authorization"]="Bearer %s" % (token)
        return session

#------------------------------------------------------------------------------------------------------------
def getTAPService(session):
    """Return the RSP TAP service.

    """

    try:
        from lsst.rsp import get_tap_service
        return get_tap_service("tap")
    except ImportError:
        import pyvo
        return pyvo.dal.TAPService("https://data.lsst.cloud/api/tap", session = session)

#------------------------------------------------------------------------------------------------------------
def findDP2Coadds(RADeg, decDeg, sizeArcmin, bands, session):
    """Return the DP2 deep_coadd ObsCore rows (one tract only) overlapping a sizeArcmin box around
    RADeg, decDeg.

    """

    # The TAP service does not support INTERSECTS, and its CONTAINS misses patches whose s_region crosses
    # RA = 0. So fetch patches whose centres are near the target (a patch is ~11 arcmin across), then
    # keep those whose outline contains any point of a grid (spacing <= 3 arcmin, much smaller than a
    # patch) covering the box. Positions are handled as offsets from the target, so RA wrapping is OK.
    cosDec=np.cos(np.radians(decDeg))
    searchDeg=sizeArcmin/60/np.sqrt(2) + 0.2
    RAHalfWidth=min(searchDeg/max(cosDec, 1e-3), 180)
    RAMin, RAMax=(RADeg - RAHalfWidth) % 360, (RADeg + RAHalfWidth) % 360
    RACut=("s_ra BETWEEN %.6f AND %.6f" if RAMin < RAMax else "(s_ra >= %.6f OR s_ra <= %.6f)") % (RAMin, RAMax)
    bandList=", ".join("'%s'" % (b) for b in bands)
    query="""SELECT DISTINCT lsst_band, lsst_tract, lsst_patch, s_ra, s_dec, s_region, access_url
             FROM ivoa.ObsCore
             WHERE obs_collection = 'LSST.DP2' AND dataproduct_subtype = 'lsst.deep_coadd'
             AND lsst_band IN (%s) AND s_dec BETWEEN %.6f AND %.6f AND %s""" \
             % (bandList, decDeg - searchDeg, decDeg + searchDeg, RACut)
    tab=getTAPService(session).search(query).to_table()
    if len(tab) == 0:
        return tab

    nGrid=int(np.ceil(sizeArcmin/3)) + 1
    offsets=np.linspace(-sizeArcmin/2, sizeArcmin/2, nGrid)/60
    gridPoints=np.array([(dx, dy) for dx in offsets for dy in offsets])   # (dRA*cosDec, dDec) in deg
    keep=[]
    for region in tab['s_region']:
        coords=np.array(region.split()[2:], dtype = float).reshape(-1, 2)   # "POLYGON ICRS ra1 dec1 ..."
        dRA=(coords[:, 0] - RADeg + 180) % 360 - 180
        poly=Path(np.column_stack([dRA*cosDec, coords[:, 1] - decDeg]))
        keep.append(poly.contains_points(gridPoints).any())
    tab=tab[np.array(keep, dtype = bool)]
    if len(tab) == 0:
        return tab

    # Patches from different tracts are on different pixel grids, so use only the tract whose
    # nearest patch centre is closest to the target
    sep=SkyCoord(tab['s_ra'], tab['s_dec'], unit = "deg").separation(SkyCoord(RADeg, decDeg, unit = "deg"))

    return tab[tab['lsst_tract'] == tab['lsst_tract'][np.argmin(sep)]]

#------------------------------------------------------------------------------------------------------------
def fetchCutout(accessURL, RADeg, decDeg, sizeArcmin, session):
    """Fetch a cutout from one coadd via the SODA service. Returns the IMAGE HDU.

    """

    from pyvo.dal.adhoc import DatalinkResults, SodaQuery
    dl=DatalinkResults.from_result_url(accessURL, session = session)
    sq=SodaQuery.from_resource(dl, dl.get_adhocservice_by_id("cutout-sync"), session = session)
    sq.circle=(RADeg*u.deg, decDeg*u.deg, (sizeArcmin/60/2)*u.deg)   # Returns the bounding square

    return fits.open(io.BytesIO(sq.execute_stream().read()))['IMAGE']

#------------------------------------------------------------------------------------------------------------
def stitchCutouts(hdus, RADeg, decDeg, sizeArcmin):
    """Put cutouts from patches of the same tract onto one sizeArcmin x sizeArcmin array.

    Uses the alternate 'A' WCS in each cutout header, which gives the position on the tract pixel grid.
    Returns (data, wcs). Pixels with no coverage are NaN.

    """

    ref=hdus[0].header
    refWCS=WCS(ref)
    pixScaleArcsec=proj_plane_pixel_scales(refWCS)[0]*3600
    n=int(round(sizeArcmin*60/pixScaleArcsec)) + 1
    tractOrigin=lambda h: (int(h['CRVAL1A'] - h['CRPIX1A']) + 1, int(h['CRVAL2A'] - h['CRPIX2A']) + 1)
    rx0, ry0=tractOrigin(ref)
    px, py=refWCS.world_to_pixel_values(RADeg, decDeg)
    ox, oy=rx0 + int(round(float(px))) - n//2, ry0 + int(round(float(py))) - n//2   # Output origin, tract pixels

    out=np.full((n, n), np.nan, dtype = np.float32)
    for h in hdus:
        x0, y0=tractOrigin(h.header)
        ny, nx=h.data.shape
        xa, xb=max(x0, ox), min(x0 + nx, ox + n)
        ya, yb=max(y0, oy), min(y0 + ny, oy + n)
        if xa >= xb or ya >= yb:
            continue
        dst=out[ya - oy:yb - oy, xa - ox:xb - ox]
        src=h.data[ya - y0:yb - y0, xa - x0:xb - x0]
        fill=np.isnan(dst)
        dst[fill]=src[fill]

    outWCS=refWCS.deepcopy()
    outWCS.wcs.crpix=refWCS.wcs.crpix - np.array([ox - rx0, oy - ry0])

    return out, outWCS

#------------------------------------------------------------------------------------------------------------
def getDP2Images(RADeg, decDeg, sizeArcmin, bands, cacheDir = CACHE_DIR, refetch = False):
    """Return a dictionary of stitched DP2 deep_coadd cutouts (nJy) for each band, and their (common) WCS.
    Cutouts are cached as FITS images in cacheDir. Exits if any band has no DP2 coadd at this position.

    """

    os.makedirs(cacheDir, exist_ok = True)
    cachePath=lambda b: os.path.join(cacheDir, "DP2_%.6f_%.6f_%.2f_%s.fits" % (RADeg, decDeg, sizeArcmin, b))

    images={}
    imgWCS=None
    toFetch=[b for b in bands if refetch or not os.path.exists(cachePath(b))]
    for b in bands:
        if b not in toFetch:
            with fits.open(cachePath(b)) as img:
                images[b]=img[0].data
                imgWCS=WCS(img[0].header)

    if len(toFetch) > 0:
        session=getSession()
        coaddTab=findDP2Coadds(RADeg, decDeg, sizeArcmin, toFetch, session)
        missing=[b for b in toFetch if len(coaddTab) == 0 or b not in coaddTab['lsst_band']]
        if len(missing) > 0:
            print("No DP2 deep_coadd at RA, dec = %.5f, %.5f in band(s) %s - no image made."
                  % (RADeg, decDeg, ", ".join(missing)))
            sys.exit(1)
        print("Tract %d, patches %s" % (coaddTab['lsst_tract'][0], sorted(set(int(p) for p in coaddTab['lsst_patch']))))
        for b in toFetch:
            rows=coaddTab[coaddTab['lsst_band'] == b]
            hdus=[fetchCutout(row['access_url'], RADeg, decDeg, sizeArcmin, session) for row in rows]
            images[b], bandWCS=stitchCutouts(hdus, RADeg, decDeg, sizeArcmin)
            print("%s: %d patch(es), %.0f%% of cutout covered" % (b, len(hdus), 100*np.isfinite(images[b]).mean()))
            fits.PrimaryHDU(images[b], header = bandWCS.to_header()).writeto(cachePath(b), overwrite = True)
            imgWCS=bandWCS

    return images, imgWCS

#------------------------------------------------------------------------------------------------------------
def makeRGB(images, bands, stretch, Q):
    """Make an RGB array from 1, 2 or 3 bands (Lupton et al. 2004 asinh stretch). The reddest band is R and
    the bluest is B. With 2 bands, G is their mean; with 1 band the image is greyscale.

    """

    scaled=[np.nan_to_num(images[b])*BAND_SCALES[b] for b in bands]   # bands are in order blue -> red
    if len(scaled) == 1:
        R=G=B=scaled[0]
    elif len(scaled) == 2:
        B, R=scaled
        G=(B + R)/2
    else:
        B, G, R=scaled

    return make_lupton_rgb(R, G, B, stretch = stretch, Q = Q, minimum = 0)

#------------------------------------------------------------------------------------------------------------
def findMeerKATImage(RADeg, decDeg):
    """Return the URL of the Mining with MeerKAT EDR image whose pointing centre is nearest RADeg, decDeg.

    """

    tab=atpy.Table().read(download_file(EDR_URL + "images_EDR.fits", cache = True, show_progress = False))
    sepDeg=SkyCoord(tab['centre_RADeg'], tab['centre_decDeg'], unit = "deg").separation(SkyCoord(RADeg, decDeg, unit = "deg")).deg
    row=tab[sepDeg.argmin()]
    print("MeerKAT image: %s (centre %.2f deg from target)" % (row['object'], sepDeg.min()))
    if sepDeg.min() > 1.0:
        print("WARNING: target is %.2f deg from the nearest MeerKAT EDR pointing centre" % (sepDeg.min()))

    return EDR_URL + row['path']

#------------------------------------------------------------------------------------------------------------
def loadContourImage(path, RADeg, decDeg, sizeArcmin):
    """Load a FITS image (local path, URL, or 'meerkat' for the nearest Mining with MeerKAT EDR image) and
    return a cutout around RADeg, decDeg (with a margin, so that contours reach the plot edges) as
    (data, wcs). URLs are cached by astropy (in ~/.astropy/cache).

    """

    if path.lower() == "meerkat":
        path=findMeerKATImage(RADeg, decDeg)
    if path.startswith("http://") or path.startswith("https://"):
        path=download_file(path, cache = True, show_progress = True)

    with fits.open(path) as img:
        # First HDU with image data; drop degenerate freq/Stokes axes
        hdu=[h for h in img if h.data is not None and h.data.ndim >= 2][0]
        data=hdu.data.squeeze()
        wcs=WCS(hdu.header).celestial
        cutout=Cutout2D(data, SkyCoord(RADeg, decDeg, unit = "deg"), 1.2*sizeArcmin*u.arcmin, wcs = wcs,
                        mode = "partial", fill_value = np.nan)

    return cutout.data, cutout.wcs

#------------------------------------------------------------------------------------------------------------
def addContours(ax, imgWCS, contourData, contourWCS, sigmaLevels, colour, width, smoothArcsec):
    """Draw contours of contourData (levels are multiples of its 3-sigma-clipped RMS) on ax, whose pixel
    grid is described by imgWCS. Contours are drawn at the optical pixel positions of each contour image
    pixel, so the contour image doesn't need regridding.

    """

    _, _, rms=sigma_clipped_stats(contourData, sigma = 3, maxiters = 10)
    print("Contour image 3-sigma-clipped RMS: %.3g (image units, e.g. Jy/beam)" % (rms))
    levels=np.array(sigmaLevels)*rms

    data=np.nan_to_num(contourData)
    if smoothArcsec > 0:
        pixScaleArcsec=proj_plane_pixel_scales(contourWCS)[0]*3600
        data=ndimage.gaussian_filter(data, smoothArcsec/pixScaleArcsec)

    ny, nx=data.shape
    yy, xx=np.mgrid[0:ny, 0:nx]
    world=contourWCS.pixel_to_world(xx, yy)
    X, Y=imgWCS.world_to_pixel(world)
    ax.contour(X, Y, data, levels = levels, colors = colour, linewidths = width)

    return rms

#------------------------------------------------------------------------------------------------------------
def chooseScaleBarArcmin(sizeArcmin):
    """A round scale bar length, roughly a quarter of the image width.

    """

    choices=[0.1, 0.25, 0.5, 1, 2, 5, 10, 15, 20, 30, 60]
    ok=[c for c in choices if c <= sizeArcmin/4]

    return ok[-1] if len(ok) > 0 else choices[0]

#------------------------------------------------------------------------------------------------------------
def addScaleBar(ax, imgWCS, nPix, scaleBarArcmin, colour = "yellow", fontSize = 16):
    """Draw a scale bar in the bottom right corner, labelled with its angular length.

    """

    pixScaleArcmin=proj_plane_pixel_scales(imgWCS)[0]*60
    lenPix=scaleBarArcmin/pixScaleArcmin
    pad=0.05*nPix
    x1=nPix - pad
    x0=x1 - lenPix
    y=pad
    ax.plot([x0, x1], [y, y], '-', color = colour, lw = 2.0, solid_capstyle = 'butt')

    if scaleBarArcmin >= 1:
        label="%g'" % (scaleBarArcmin)
    else:
        label='%g"' % (scaleBarArcmin*60)
    ax.text((x0 + x1)/2, y + 0.012*nPix, label, ha = 'center', va = 'bottom', color = colour, size = fontSize)

#------------------------------------------------------------------------------------------------------------
def makePlot(rgb, imgWCS, name, outFileName, redshift = None, showName = True, showRedshift = True,
             showScaleBar = True, scaleBarArcmin = None, sizeArcmin = None, contourImage = None,
             contourSigmaLevels = (3, 5, 10, 20, 40, 80, 160), contourColour = "cyan", contourWidth = 1.0,
             contourSmoothArcsec = 0.0, RADeg = None, decDeg = None):
    """Plot the RGB image with no coordinate axes, and optional name/redshift labels, scale bar and
    contours, and save to outFileName (format set by the extension, e.g. .png or .jpg).

    """

    nPix=rgb.shape[0]
    fig=plt.figure(figsize = (FIG_SIZE, FIG_SIZE))
    ax=fig.add_axes([0, 0, 1, 1])
    ax.imshow(rgb, origin = "lower", interpolation = "nearest")
    ax.set_autoscale_on(False)   # Keep the optical image extent when adding contours
    ax.set_axis_off()

    if contourImage is not None:
        contourData, contourWCS=loadContourImage(contourImage, RADeg, decDeg, sizeArcmin)
        addContours(ax, imgWCS, contourData, contourWCS, contourSigmaLevels, contourColour, contourWidth,
                    contourSmoothArcsec)

    if showName:
        label=name.replace("_", " ").replace("/", "\n/")   # Split e.g. supercluster names across lines
        fig.text(0.025, 0.975, label, ha = 'left', va = 'top', size = 24, color = 'yellow')
    if showRedshift and redshift is not None:
        fig.text(0.025, 0.03, "z = %.3f" % (redshift) if redshift < 0.1 else "z = %.2f" % (redshift),
                 ha = 'left', va = 'bottom', size = 24, color = 'white')
    if showScaleBar:
        if scaleBarArcmin is None:
            scaleBarArcmin=chooseScaleBarArcmin(sizeArcmin)
        addScaleBar(ax, imgWCS, nPix, scaleBarArcmin)

    fig.savefig(outFileName, dpi = DPI)
    plt.close(fig)
    print("Wrote %s" % (outFileName))

#------------------------------------------------------------------------------------------------------------
def parseBands(bands):
    """Check the bands string (e.g. 'grz') and return it sorted blue -> red.

    """

    bands=bands.strip()
    bad=[b for b in bands if b not in BAND_ORDER]
    if len(bad) > 0 or len(bands) == 0 or len(set(bands)) != len(bands):
        raise argparse.ArgumentTypeError("bands must be 1-3 distinct letters from '%s' (got '%s')" % (BAND_ORDER, bands))
    if len(bands) > 3:
        raise argparse.ArgumentTypeError("at most 3 bands can be used for an RGB image (got '%s')" % (bands))

    return "".join(sorted(bands, key = BAND_ORDER.index))

#------------------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    parser=argparse.ArgumentParser(description = "Make a Rubin DP2 deep_coadd colour cut out image, optionally with contours overlaid.")
    parser.add_argument("RADeg", type = float, help = "RA of the image centre (decimal degrees)")
    parser.add_argument("decDeg", type = float, help = "Dec of the image centre (decimal degrees)")
    parser.add_argument("sizeArcmin", type = float, help = "Width (and height) of the image (arcmin)")
    parser.add_argument("name", help = "Object name, used for the label and the output file name")
    parser.add_argument("-b", "--bands", type = parseBands, default = "grz",
                        help = "Bands to use (1-3 of ugrizy, e.g. grz or gri; reddest = R). Default: grz")
    parser.add_argument("-z", "--redshift", type = float, default = None,
                        help = "Redshift, shown as a label")
    parser.add_argument("-o", "--output", default = None,
                        help = "Output file name (.png or .jpg). Default: from the name, as a .png")
    parser.add_argument("--name-label", action = argparse.BooleanOptionalAction, default = True,
                        help = "Show the object name on the image (default: yes)")
    parser.add_argument("--redshift-label", action = argparse.BooleanOptionalAction, default = True,
                        help = "Show the redshift on the image, if given (default: yes)")
    parser.add_argument("--scale-bar", action = argparse.BooleanOptionalAction, default = True,
                        help = "Show a scale bar (default: yes)")
    parser.add_argument("--scale-bar-arcmin", type = float, default = None,
                        help = "Scale bar length in arcmin (default: a round number ~1/4 of the image width)")
    parser.add_argument("-c", "--contours", default = None,
                        help = "FITS image (path or URL) to overlay as contours, or 'meerkat' to use the "
                               "nearest Mining with MeerKAT EDR image")
    parser.add_argument("--contour-levels", type = float, nargs = "+", default = [3, 5, 10, 20, 40, 80, 160],
                        help = "Contour levels, in units of the contour image's clipped RMS (default: 3 5 10 20 40 80 160)")
    parser.add_argument("--contour-colour", default = "cyan", help = "Contour colour (default: cyan)")
    parser.add_argument("--contour-width", type = float, default = 1.0, help = "Contour line width (default: 1)")
    parser.add_argument("--contour-smooth-arcsec", type = float, default = 0.0,
                        help = "Gaussian smoothing (sigma, arcsec) applied to the contour image (default: none)")
    parser.add_argument("--stretch", type = float, default = 30, help = "make_lupton_rgb stretch (default: 30)")
    parser.add_argument("--Q", type = float, default = 6, help = "make_lupton_rgb Q (default: 6)")
    parser.add_argument("--cache-dir", default = CACHE_DIR, help = "Where DP2 cutouts are cached (default: %s)" % (CACHE_DIR))
    parser.add_argument("--refetch", action = "store_true", help = "Re-download DP2 cutouts even if cached")
    args=parser.parse_args()

    outFileName=args.output
    if outFileName is None:
        outFileName=args.name.replace(" ", "_").replace("/", "and").replace(".", "p") + ".png"
    if os.path.splitext(outFileName)[-1].lower() not in [".png", ".jpg", ".jpeg"]:
        parser.error("output file name must end in .png or .jpg")

    images, imgWCS=getDP2Images(args.RADeg, args.decDeg, args.sizeArcmin, args.bands, cacheDir = args.cache_dir,
                                refetch = args.refetch)
    rgb=makeRGB(images, args.bands, args.stretch, args.Q)
    makePlot(rgb, imgWCS, args.name, outFileName, redshift = args.redshift, showName = args.name_label,
             showRedshift = args.redshift_label, showScaleBar = args.scale_bar, scaleBarArcmin = args.scale_bar_arcmin,
             sizeArcmin = args.sizeArcmin, contourImage = args.contours, contourSigmaLevels = args.contour_levels,
             contourColour = args.contour_colour, contourWidth = args.contour_width,
             contourSmoothArcsec = args.contour_smooth_arcsec, RADeg = args.RADeg, decDeg = args.decDeg)
