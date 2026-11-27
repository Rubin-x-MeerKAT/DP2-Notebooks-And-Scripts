"""

Figuring out the DP2 footprint, for prioritizing MeerKAT observations to process

Uses where DP2 deep_coadd images actually exist, from the RSP ObsCore table (ivoa.ObsCore,
obs_collection = 'LSST.DP2'), instead of the DP2 Visit table (getDP2Footprint_visitTable.py).

Each coadd patch has an outline (s_region). The patch outlines for each band are rasterised onto a HEALPix
map (pixels whose centres fall inside the patch). The footprint is the sky where every band in
requiredBands has a coadd, stored as a MOC.

Some areas with DP2 visits have no DP2 coadds (e.g. the XMM-LSS field around RA 37 deg, Dec -4 deg).
ObsCore has no depth information for coadds (exposure time is empty), so this footprint can't be cut on the
number of visits.

Given a catalog, writes the subset of it that lies within the footprint to InDP2Coadd_<catalog name>.fits
(a FITS table) in the current directory.

Run with -h for options. Examples:

    python getWithinDP2CoaddFootprint.py DR6_cluster-catalog_v1.0.fits
    python getWithinDP2CoaddFootprint.py images_EDR.fits --ra-col centre_RADeg --dec-col centre_decDeg
    python getWithinDP2CoaddFootprint.py MarcelinaPhDClusters.fits --bands gri

"""

import os
import argparse
import numpy as np
import healpy as hp
import astropy.table as atpy
from astropy.coordinates import Angle
from mocpy import MOC

PATCH_CACHE="DP2_coadd_patches.fits"
MOC_CACHE_TEMPLATE="DP2_coadd_footprint_moc_%s.fits"   # %s is the required bands, e.g. "grz"
NSIDE=1024                # HEALPix resolution (~3.4 arcmin pixels; a patch is ~11 arcmin)
ALL_BANDS=["u", "g", "r", "i", "z", "y"]
REQUIRED_BANDS=["g", "r", "z"]

_coverageMaps={}          # Per-band coverage maps, built once per run (keyed by nside)

#------------------------------------------------------------------------------------------------------------
def getTAPService(tokenFile = "~/.rsp-tap.token"):
    """Return the RSP TAP service.

    """

    try:
        # Works when running on the RSP itself
        from lsst.rsp import get_tap_service
        service=get_tap_service("tap")
    except ImportError:
        # Works from anywhere, given an RSP access token (https://data.lsst.cloud/auth/tokens)
        import pyvo
        import requests
        with open(os.path.expanduser(tokenFile)) as inFile:
            token=inFile.read().strip()
        session=requests.Session()
        session.headers["Authorization"]="Bearer %s" % (token)
        service=pyvo.dal.TAPService("https://data.lsst.cloud/api/tap", session = session)

    return service

#------------------------------------------------------------------------------------------------------------
def runQuery(service, query):
    """Run an async TAP query and return the result as an astropy table.

    """

    job=service.submit_job(query)
    job.run()
    job.wait(phases = ['COMPLETED', 'ERROR'])
    job.raise_if_error()
    tab=job.fetch_result().to_table()
    job.delete()

    return tab

#------------------------------------------------------------------------------------------------------------
def fetchDP2CoaddPatches(cacheFileName = PATCH_CACHE, tokenFile = "~/.rsp-tap.token"):
    """Query the DP2 deep_coadd patches from ObsCore (or read a cached copy). Returns a table with one row
    per (tract, patch), holding the patch outline corners (ra, dec; 4 each) and a has_<band> flag per band.

    ObsCore has one row per patch per band (~925,000 rows). To keep memory use down, the patch outlines
    (the same in every band) are fetched once per patch with DISTINCT, and the band list is fetched
    separately without them. The first run takes about two minutes.

    """

    if os.path.exists(cacheFileName):
        return atpy.Table().read(cacheFileName)

    service=getTAPService(tokenFile = tokenFile)
    where="FROM ivoa.ObsCore WHERE obs_collection = 'LSST.DP2' AND dataproduct_subtype = 'lsst.deep_coadd'"
    polyTab=runQuery(service, "SELECT DISTINCT lsst_tract, lsst_patch, s_region %s" % (where))
    bandTab=runQuery(service, "SELECT lsst_band, lsst_tract, lsst_patch %s" % (where))

    # One row per (tract, patch)
    patchKey=lambda t: np.array(t['lsst_tract'], dtype = np.int64)*100000 + np.array(t['lsst_patch'], dtype = np.int64)
    uKey, first=np.unique(patchKey(polyTab), return_index = True)
    patchTab=atpy.Table()
    patchTab['tract']=np.array(polyTab['lsst_tract'][first], dtype = np.int32)
    patchTab['patch']=np.array(polyTab['lsst_patch'][first], dtype = np.int32)
    # s_region is "POLYGON ICRS ra1 dec1 ra2 dec2 ra3 dec3 ra4 dec4"
    corners=np.array([str(s).split()[2:] for s in polyTab['s_region'][first]], dtype = float).reshape(-1, 4, 2)
    patchTab['ra']=corners[..., 0]
    patchTab['dec']=corners[..., 1]
    del polyTab

    rowPatch=np.searchsorted(uKey, patchKey(bandTab))   # Patch index of each (band, patch) row
    bandArr=np.array(bandTab['lsst_band'], dtype = str)
    for b in ALL_BANDS:
        has=np.zeros(len(uKey), dtype = bool)
        has[rowPatch[bandArr == b]]=True
        patchTab['has_%s' % (b)]=has
    del bandTab
    patchTab.write(cacheFileName, overwrite = True)

    return patchTab

#------------------------------------------------------------------------------------------------------------
def getCoverageMaps(nside = NSIDE):
    """Return a dictionary of per-band HEALPix (nested) boolean maps, True where a DP2 deep_coadd exists.
    A patch counts as covered in a band if a coadd exists for it, even if only partly filled with data.

    """

    if nside in _coverageMaps:
        return _coverageMaps[nside]

    patchTab=fetchDP2CoaddPatches()

    # Pixels inside each patch (hp.ang2vec needs 1-D input, so flatten and reshape)
    cornerVecs=hp.ang2vec(np.ravel(patchTab['ra']), np.ravel(patchTab['dec']), lonlat = True).reshape(-1, 4, 3)
    patchPix=[hp.query_polygon(nside, v, nest = True) for v in cornerVecs]
    allPix=np.concatenate(patchPix)
    pixPatch=np.repeat(np.arange(len(patchPix)), [len(p) for p in patchPix])   # Patch index of each allPix entry

    maps={}
    for b in ALL_BANDS:
        m=np.zeros(hp.nside2npix(nside), dtype = bool)
        m[allPix[np.array(patchTab['has_%s' % (b)])[pixPatch]]]=True
        maps[b]=m
    _coverageMaps[nside]=maps

    return maps

#------------------------------------------------------------------------------------------------------------
def parseBands(requiredBands):
    """Return requiredBands as a list of single band names. Accepts e.g. "grz" or ["g", "r", "z"].
    None or empty means coadds in any band.

    """

    if requiredBands is None:
        return []
    bands=list(requiredBands)
    for b in bands:
        if b not in ALL_BANDS:
            raise ValueError("Unknown band '%s'; available bands: %s" % (b, ALL_BANDS))

    return bands

#------------------------------------------------------------------------------------------------------------
def getDP2CoaddMOC(requiredBands = REQUIRED_BANDS, nside = NSIDE, cacheFileName = None):
    """Build (or read a cached) MOC of the sky where every band in requiredBands has a DP2 deep_coadd.
    If requiredBands is None or empty, the MOC covers the sky with a coadd in any band.

    By default the MOC is cached as DP2_coadd_footprint_moc_<bands>.fits (e.g. _grz, or _any), so that
    different band cuts don't overwrite each other. Delete the cache if nside is changed.

    """

    bands=parseBands(requiredBands)
    if cacheFileName is None:
        label="".join(b for b in ALL_BANDS if b in bands) if len(bands) > 0 else "any"
        cacheFileName=MOC_CACHE_TEMPLATE % (label)

    if os.path.exists(cacheFileName):
        return MOC.from_fits(cacheFileName)

    maps=getCoverageMaps(nside = nside)
    if len(bands) > 0:
        mask=np.all([maps[b] for b in bands], axis = 0)
    else:
        mask=np.any(list(maps.values()), axis = 0)

    ipix=np.where(mask)[0].astype(np.uint64)
    order=hp.nside2order(nside)
    moc=MOC.from_healpix_cells(ipix, np.full(len(ipix), order, dtype = np.uint8), max_depth = order)
    moc.save(cacheFileName, format = "fits", overwrite = True)

    return moc

#------------------------------------------------------------------------------------------------------------
def getTabInDP2Coadds(tabFileName, RACol = 'RADeg', decCol = 'decDeg', requiredBands = REQUIRED_BANDS):
    """Given a path to a catalog table, return the subset that fall within the DP2 coadd footprint (where
    every band in requiredBands has a deep_coadd) as an astropy table.

    requiredBands can be a list (e.g. ["g", "r", "z"]) or a string (e.g. "grz"). None or empty means a
    coadd in any band.

    """

    moc=getDP2CoaddMOC(requiredBands = requiredBands)

    tab=atpy.Table().read(tabFileName)
    lon=Angle(tab[RACol], unit = "deg")
    lat=Angle(tab[decCol], unit = "deg")
    hitMask=moc.contains_lonlat(lon, lat)
    DP2Tab=tab[hitMask]

    return DP2Tab

#------------------------------------------------------------------------------------------------------------
def parseBandsArg(bands):
    """Check a --bands string (e.g. 'grz') and return it sorted blue -> red. 'any' means a coadd in any band.

    """

    bands=bands.strip()
    if bands.lower() == "any":
        return "any"
    bad=[b for b in bands if b not in ALL_BANDS]
    if len(bad) > 0 or len(bands) == 0 or len(set(bands)) != len(bands):
        raise argparse.ArgumentTypeError("bands must be distinct letters from '%s', or 'any' (got '%s')"
                                         % ("".join(ALL_BANDS), bands))

    return "".join(sorted(bands, key = ALL_BANDS.index))

#------------------------------------------------------------------------------------------------------------
# Main

if __name__ == "__main__":

    parser=argparse.ArgumentParser(description = "Write the subset of a catalog that lies within the DP2 deep_coadd footprint to InDP2Coadd_<catalog name>.fits.")
    parser.add_argument("catalogFileName", help = "Catalog file (any format astropy can read, e.g. FITS or CSV)")
    parser.add_argument("-b", "--bands", type = parseBandsArg, default = "".join(REQUIRED_BANDS),
                        help = "Bands that must all have coadds (any of ugrizy, e.g. grz or gri), or 'any' for a"
                               " coadd in any band. Default: %s" % ("".join(REQUIRED_BANDS)))
    parser.add_argument("--ra-col", default = "RADeg", help = "RA column name, in decimal degrees (default: RADeg)")
    parser.add_argument("--dec-col", default = "decDeg", help = "Dec column name, in decimal degrees (default: decDeg)")
    args=parser.parse_args()

    if not os.path.exists(args.catalogFileName):
        parser.error("catalog file not found: %s" % (args.catalogFileName))

    requiredBands=None if args.bands == "any" else args.bands
    DP2Tab=getTabInDP2Coadds(args.catalogFileName, RACol = args.ra_col, decCol = args.dec_col,
                             requiredBands = requiredBands)

    outFileName="InDP2Coadd_%s.fits" % (os.path.splitext(os.path.basename(args.catalogFileName))[0])
    DP2Tab.write(outFileName, overwrite = True)
    print("Wrote %s (%d of the input objects are within the DP2 coadd footprint)" % (outFileName, len(DP2Tab)))
