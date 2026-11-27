DP2 Notebooks And Scripts
=========================

This repository contains Jupyter notebooks and scripts that use the 
`Rubin Science Platform (RSP) <https://data.lsst.cloud/>`__ API that combine the 
`Mining with MeerKAT Early Data Release (EDR) <https://extragalactic.phys.wits.ac.za/mining-with-meerkat/EDR/docs/>`_ 
with Rubin Observatory LSST Data Preview 2 (DP2).

Some of the tools in here have more generic uses - e.g., checking which objects
in a catalog are within the DP2 co-add images; making RGB cut-outs around
given object positions.

This repository is under active development. Some of the tools here may be
cleaned up and put into a helper package in the future to streamline the code.

Notebooks
---------

The notebooks are meant to be run on the `RSP <https://data.lsst.cloud/nb>`__.

``Mining-with-MeerKAT_DP2Footprint.ipynb``
    Selects the MeerKAT EDR images and L-band source catalog entries that lie
    inside the DP2 footprint, with the footprint built from the DP2 ``Visit``
    table.

    Outputs: ``DP2_footprint_moc.fits``, ``InDP2_images_EDR.fits``,
    ``InDP2_survey_catalog_L_EDR.fits`` and a sky plot,
    ``DP2_footprint_MeerKAT_EDR.png``.

``Mining-with-MeerKAT_DP2CoaddFootprint.ipynb``
    Does the same selection, but the footprint is where DP2 ``deep_coadd``
    images actually exist. This is more useful for Early DP2, which does not
    contain the single epoch images (which will only be released in Full DP2
    in December 2026).

    Outputs: ``DP2_coadd_footprint_moc.fits``, ``InDP2Coadd_images_EDR.fits``,
    ``InDP2Coadd_survey_catalog_L_EDR.fits`` and a sky plot,
    ``DP2_coadd_footprint_MeerKAT_EDR.png``.

``Mining-with-MeerKAT_Images_CutOut.ipynb``
    Downloads a MeerKAT image from the EDR web server and makes a cutout
    around a single radio source. It then fetches DP2 ``deep_coadd`` cutouts
    of the same region, and creates an overlay plot with radio contours
    on top of the RGB optical image.


Scripts
-------

The scripts can be run locally on your machine, so long as you have set up
an `RSP <https://data.lsst.cloud/>`__ API access token
(this can be done from the Settings menu on the RSP, and should be stored as
a file ``$HOME/.rsp-tap.token``).

For all scripts you can run them with the ``-h`` flag to see a description
of usage and information on arguments and options.

``makeDP2CutOut.py``
    Makes a cutout RGB image around the given coordinates (optionally
    with a contour overlay from another image, e.g., radio contours from MeerKAT).
    Downloaded FITS images are cached, so repeated runs with the same coords
    and clip size are faster, if you want to adjust anything.

    *Example*::

        python makeDP2CutOut.py 326.6608 -3.3523 4 "ACT-CL J2146.6-0321" -z 1.16

    Output is shown in Fig. 1 below.

    .. figure:: figs/ACT-CL_J2146p6-0321.png
       :alt: DP2 co-add cutout around the z = 1.16 galaxy cluster ACT-CL J2146.6-0321 made by the makeDP2CutOut.py script
       :width: 400px
       :align: center

       **Fig. 1:** DP2 co-add cutout (4 arcmin on a side) around the z = 1.16 galaxy
       cluster ACT-CL J2146.6-0321 made by the ``makeDP2CutOut.py`` script.

    *Contour overlay example*::

        python makeDP2CutOut.py 349.3856494 0.1750374 4 "contours example" -c "https://extragalactic.phys.wits.ac.za/mining-with-meerkat/EDR/data/images/pbcorr_trim_img_1672735587_sdp_l0_1024ch_NGC7603.ms_pcalmask-MFS-image.fits"

    Output is shown in Fig. 2 below (the target matches the example in the ``Mining-with-MeerKAT_Images_CutOut.ipynb`` notebook).

    .. figure:: figs/contours_example.png
       :alt: DP2 co-add cutout around some radio sources detected in a MeerKAT image, shown as the cyan contours, made by the makeDP2CutOut.py script
       :width: 400px
       :align: center

       **Fig. 2:** DP2 co-add cutout around some radio sources detected in a MeerKAT
       image, shown as the cyan contours, made by the ``makeDP2CutOut.py`` script.

``getWithinDP2CoaddFootprint.py``
    Selects the objects in a catalog that lie within the DP2 ``deep_coadd``
    footprint (the sky where every required band has a coadd), and writes them
    to a FITS table named ``InDP2Coadd_<catalog name>.fits`` in the current
    directory. The catalog can be any format astropy can read (e.g., FITS or
    CSV). The required bands are set with ``-b`` (default ``grz``, or ``any``
    for a coadd in any band), and the RA and Dec column names (in decimal
    degrees) with ``--ra-col`` and ``--dec-col`` (defaults ``RADeg`` and
    ``decDeg``).

    The DP2 coadd patch list and the footprint MOC are cached as FITS files
    in the current directory. The first run queries the RSP and takes about
    two minutes.

    *Example*::

        python getWithinDP2CoaddFootprint.py images_EDR.fits --ra-col centre_RADeg --dec-col centre_decDeg -b gri

    This makes a catalog that contains the subset of Mining with MeerKAT EDR MeerKAT
    pointings that land within the DP2 co-add regions.

Notes
-----

* The footprint notebooks need `mocpy <https://cds-astro.github.io/mocpy/>`_.
  If it is missing, they install it into your user area. Restart the kernel
  after the install.

* The MeerKAT catalogs are read directly from the EDR web server, and astropy
  caches them after the first download (the source catalog is about 160 MB).
  The DP2 query results are cached as FITS files next to each notebook, so
  later runs are faster. To re-run a query, delete its cache file.

License
-------

This work is licensed under a
`Creative Commons Attribution 4.0 International License <https://creativecommons.org/licenses/by/4.0/>`_
(CC BY 4.0). See the ``LICENSE`` file for the full text.
