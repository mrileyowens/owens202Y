<div align="justify">

This file details the structure of this repository and the nature of its files.

## Overview

### ...

### Data

The software of this repository uses data from two primary sources: the Magellan Echellette (MagE) spectrograph mounted on ... 

## Repository structure

This repository contains three main folders: `docs`, `figs`, and `results`. The `docs` folder contains the project documentation, such as this file and the project's `readme.md`, the `figs` folder contains the figures used in the article, and the `results` folder contains various measurements and formatted tables used in the article.

Outside these folders are the `.ipynb` Jupyter notebooks used to create the figures, measurements, and tables contained in `figs` and `results`.

## Repository files

### `esc.ipynb`

This notebook is responsible for calculating the LyC escape fractions and tabulating the measurements.

The primary function of the notebook is `measure()`, which is the function that actually measures the LyC escape fractions of the MagE slit apertures.

### `galfit.ipynb`

This notebook creates a multi-panel figure showing the results of using GALFIT to remove a foreground, interloping galaxy covering part of the Sunburst Arc.

Because the LyC escape fraction calculation in `esc.ipynb` requires HST photometry from the F275W and F814W filters for each MagE slit aperture, it is necessary to remove an interloping galaxy partly covering one of the MagE slit apertures. Alex Navarre used GALFIT to remove the interloper galaxy from both filters, creating cutouts of the original images centered on the two largest arcs of the Sunburst Arc. The GALFIT data products also include the model galaxy profile.

The notebook fetches those cutouts and assembles them into a multi-panel figure. To do this, it starts by fetching the FITS files of the cutout images, and defines a new, smaller, square cutout area to later apply to these images using Astropy's `SkyCoord` object. The new cutout centers on the interloper galaxy. The notebook also creates source masks of the two filters from the cutouts of the original, unaltered images to later apply. After creating the matplotlib `Figure` and `Axes` objects for the figure, the notebook explicitly sets the vertical (`hspace`) and horizontal spacing (`wspace`) between the subplots to 0 to make the figure more compact.

The following double loop iterates over each filter (F275W and F814W) and data type (original image, model, or residual) combination, each time creating an Astropy `Cutout2D` object of the smaller, new cutout area with the defined `SkyCoord` object and the WCS of the given image. From the aforementioned source masks, the notebook calculates the pixel value color boundaries of the new cutout image, and then plots the new cutout image on the appropriate panel.

Because of difficulties encountered while trying to assign the multi-subplot figure a WCS, `plot()` opts to instead manually set the coordinate limits, ticks, and tick labels, utilizing that the drizzled pixel size of the images is 0.03″. The final size of each panel is 3″ on a side, with ticks and tick labels (except if the tick labels would appear between two panels) placed at -1″, 0″, and +1″ relative to the center on both coordinate axes. The notebook sets the aspect ratio of each panel to guarantee it appears square, and adds annotations to each panel indicating the image's filter and data type. Finally, the notebook adds coordinate axis labels to the figure and explicitly sets the tick parameters of each panel to ensure there are no tick labels for ticks facing other panels, and that each side of every panel has inward-facing ticks. The `plot()` function then saves the figure as `foreground_galaxy_galfit_modeling_results.pdf` to `figs/`.

### `lya_and_lyc_maps.ipynb`

This notebook plots the Lyα and LyC maps of the two largest arcs of the Sunburst Arc in a multi-panel figure.

In the `plot()` function the notebook fetches the common WCS of the HST images and creates two rotated WCSs, such that either arc appears horizontal. Then, for each map to plot, the function gets the image, reprojects it to the rotated WCSs, and plots the rotated image in the corresponding panel. It continues by setting the coordinate limits of the panels and disabling their ticks and labels. The remaining code adds labels indicating the arcs the panels show and the type of map, as well as compasses, scalebars, the footprints of the MagE slit apertures, and labels of the images of source plane clumps. The notebook saves the final figure as `lya_and_lyc_maps.pdf` in the `figs/` folder.

### `lya_nb_m3.ipynb`

### `lya.ipynb`

### `map.ipynb`

This notebook creates the primary on-sky map figure of the article, showing the Sunburst Arc and its cluster field, the MagE slit apertures, and zoom-in insets on the two largest arcs (those targeted by the MagE observations) with labeled clump numbers.

The notebook's primary function `plot()` starts by instantiating a figure with the common WCS of the HST images and plots the HST/WFC3 F606W image of the Sunburst Arc and cluster field. It then sets the coordinate boundaries and labels of the figure. To create the zoom-in inset axes, the notebook first defines the position of the inset axes in the main axes coordinate system, creates new WCSs by rotating the main image's WCS using `rotate_wcs()` so that the two largest arcs appear either horizontal or vertical, depending on the set of inset axes they appear in, adds the inset axes to the main panel with the rotated WCSs with `add_inset_axes()`, plots the accordingly rotated main image in the inset axes, and trims the inset axes coordinate limits to focus on the arcs using `set_inset_axes_limits()`. The notebook then disables the ticks and labels of the inset axes. The notebook executes `add_box()` several times to add footprints of the displayed area in the inset axes on the main image, as well as highlight the other arc segments of the Sunburst Arc. For each of these boxes in the main image and the corresponding inset axes, the notebook labels them according to the terminology of the lens model of [Sharon et al. (2022)](https://doi.org/10.3847/1538-4357/ac927a) (ApJ, 941, 203). For each of the MagE slit apertures, the notebook then instantiates the aperture as a regions `RectangleSkyRegion()` object, converts it to a pixel-based representation given the WCS of each panel the aperture appears in, and then plots and labels the footprint. For each image of a source plane clump that appears in the inset axes, `plot()` circles (or boxes, for one image) and labels the clump. Following this, `plot()` concludes by annotating the main axes with a compass, adds scalebars to the main and inset axes, and finally saves the figure as `sunburst_arc_and_cluster_field.pdf` to the `figs` folder.

### `masks.ipynb`

This notebook makes masks of the MagE slit apertures and the two largest arcs of the Sunburst Arc in the world coordinate system (WCS) of the latest reduction of the HST images used in the article. These masks are critical for portions of analyses in other notebooks.

When making photometric measurements to compare to the spectroscopic properties of the MagE slit apertures, it is important to measure those photometries from the footprints of the apertures, which requires masks of the apertures. Both `esc.ipynb` and `nb.ipynb` make photometric measurements that utilize the MagE slit aperture footprints and a mask of the Sunburst Arc. It is more efficient to create and save those masks explicitly and have those notebooks fetch them, rather than compute the masks each time either notebook executes. It also makes it easier to share those masks with collaborators if necessary.

The `masks.ipynb` notebook has a sole function, `make()`, which constructs the aforementioned masks. `make()` starts by fetching the WCS of one of the HST images of our collaboration's latest reduction (all the images of the latest reduction share a common WCS, so the choice of image is arbitrary).

Then, for each MagE slit aperture, the notebook creates a regions `RectangleSkyRegion` object of the aperture (since the apertures are rectangular) from an Astropy `SkyCoord` object of the celestical coordinates of the geometric center of the aperture and the aperture's width and rotation on the sky, converting it to a pixel representation as a `RectanglePixelRegion` using the HST image's WCS. From the pixel representation, the notebook creates a mask of the aperture matching the shape of the HST image with the `to_mask()` method of the new `RectanglePixelRegion` object, and then the `to_image()` method of the resulting `RegionMask` object. Using the WCS of the HST image, the notebook creates a new FITS header from the WCS with the `to_header()` method of the `WCS` object, adding several new keywords to the header describing the mask. Finally, the notebook saves the mask of the slit aperture as a `.fits` file to `results/masks/`, named like `{slit_id}_mask_v5.fits`, where the `v5` indicates the version number we assigned to our collaboration's most recent reduction of the HST images.

The notebook then moves on to producing the aforementioned mask of the two largest arcs of the Sunburst Arc. It first retrieves the precursor mask, originally made by our collaboration to apply to observations of the Sunburst Arc with the Very Large Telescope's Multi Unit Spectroscopic Explorer (MUSE). Because MUSE has a much different field of view and pixel size than the HST images (which we want to apply the mask of the arc to), the notebook reprojects the arc mask to the WCS of the HST images with reproject's `reproject_interp()` function. Because (1) the original mask does not have binary values, and (2) the reprojection's interpolation also leads to non-binary values, it is necessary to explicitly convert the reprojected mask to binary, which the notebook performs by setting any positively-valued pixels to 1, and 0 otherwise. As before, the notebook creates a new FITS header from the common WCS of the HST images and inserts several relevant keywords before saving it as `arc_mask_v5.fits` in `results/masks/`.

### `nb.ipynb`

### `seeing.ipynb`

### `source_plane.ipynb()`

The `source_plane.ipynb()` notebook makes a 4-panel figure showing the source plane reconstruction of the Sunburst Arc&mdash;based on the lens model of [Sharon et al. (2022)](https://doi.org/10.3847/1538-4357/ac927a) (ApJ, 941, 203)&mdash;with 1 of the 4 non-LyC-leaking MagE slit apertures overlaid on each panel. 

To interpret the different non-LyC-leaking Lyα profiles of the MagE slit apertures, it is important to understand the geometry of those slit apertures in the source plane, which could reveal significant geometric and morphological information. To do this, Keren Sharon used her lens model to ray trace the non-LyC-leaking MagE slit apertures to the source plane of the Sunburst Arc. The notebook takes the base `.png` images created by Keren, which show each non-LyC-leaking MagE slit aperture overlaid on the labeled source plane reconstruction of the Sunburst Arc, trims them to be square, centered on the relevant portion of the image, plots them onto the corresponding panels, labels each panel by the ID of the MagE slit aperture ray traced onto the source plane reconstruction in the panel, and removes any ticks or tick labels from the panels before saving the figure as `source_plane_nonleaker_apertures.pdf` in `figs/`.

### `stack.ipynb`

</div>