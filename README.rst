DP2-Notebooks-And-Code
======================

This repository contains Jupyter notebooks and scripts that use the 
`Rubin Science Platform (RSP) <https://data.lsst.cloud/>`_ API that combine the 
`Mining with MeerKAT Early Data Release (EDR) <https://extragalactic.phys.wits.ac.za/mining-with-meerkat/EDR/docs/>`_ 
with Rubin Observatory LSST Data Preview 2 (DP2).

The notebooks are meant to be run on the `Rubin Science Platform (RSP) <https://data.lsst.cloud/>`_.

This repository is under active development. Some of the tools here may be
cleaned up and put into a helper package in the future to streamline the code.

Notebooks
---------

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
