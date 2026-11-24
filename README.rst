EDR-DP2-Notebooks
=================

This repository contains Jupyter notebooks that combine the
`Mining with MeerKAT Early Data Release (EDR) <https://extragalactic.phys.wits.ac.za/mining-with-meerkat/EDR/docs/>`_ with the Rubin Observatory LSST Data Preview 2 (DP2).
They are meant to be run on the Rubin Science Platform (RSP).

Notebooks
---------

``Mining-with-MeerKAT_DP2Footprint.ipynb``
    Selects the MeerKAT EDR images and L-band source catalog entries that lie
    inside the DP2 footprint, with the footprint built from the DP2 ``Visit``
    table. Each visit is treated as a 1.75° radius cone (the LSSTCam field of
    view) and added to per-band HEALPix visit-count maps. The footprint is the
    sky where every band in ``REQUIRED_BANDS`` has at least ``MIN_VISITS``
    visits, and it is stored as a MOC.

    Outputs: ``DP2_footprint_moc.fits``, ``InDP2_images_EDR.fits``,
    ``InDP2_survey_catalog_L_EDR.fits`` and a sky plot,
    ``DP2_footprint_MeerKAT_EDR.png``.

``Mining-with-MeerKAT_DP2CoaddFootprint.ipynb``
    Does the same selection, but the footprint is where DP2 ``deep_coadd``
    images actually exist. The coadd patch outlines are taken from the RSP
    ObsCore table and rasterised onto per-band HEALPix coverage maps. The
    footprint is the sky where every band in ``REQUIRED_BANDS`` has a coadd.
    This footprint leaves out areas that have DP2 visits but no coadds (for
    example, the XMM-LSS field). It cannot be cut on depth, because ObsCore
    has no depth information for coadds.

    Outputs: ``DP2_coadd_footprint_moc.fits``, ``InDP2Coadd_images_EDR.fits``,
    ``InDP2Coadd_survey_catalog_L_EDR.fits`` and a sky plot,
    ``DP2_coadd_footprint_MeerKAT_EDR.png``.

``Mining-with-MeerKAT_Images_CutOut.ipynb``
    Downloads a MeerKAT EDR image straight from the EDR web server and makes a
    cutout around a target position (a radio galaxy, by default). It then
    fetches DP2 ``deep_coadd`` cutouts of the same region in the z, r and g
    bands from the RSP image cutout (SODA) service. Cutouts that cross patch
    boundaries are stitched together, and the three bands are combined into
    an RGB colour image to compare with the radio image.

Notes
-----

* The footprint notebooks need `mocpy <https://cds-astro.github.io/mocpy/>`_.
  If it is missing, they install it into your user area. Restart the kernel
  after the install.
* The MeerKAT catalogs are read directly from the EDR web server, and astropy
  caches them after the first download (the source catalog is about 160 MB).
  The DP2 query results are cached as FITS files next to each notebook, so
  later runs are faster. To re-run a query, delete its cache file.
