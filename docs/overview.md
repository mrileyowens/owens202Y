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

<!--

This notebook is responsible for calculating the LyC escape fractions and tabulating the measurements.

The primary function of the notebook is `measure()`, which is the function that actually measures the LyC escape fractions of the MagE slit apertures.

-->

### `galfit.ipynb`

This notebook creates a multi-panel figure showing the results of using GALFIT to remove a foreground, interloping galaxy covering part of the Sunburst Arc.

Because the LyC escape fraction calculation in `esc.ipynb` requires HST photometry from the F275W and F814W filters for each MagE slit aperture, it is necessary to remove an interloping galaxy partly covering one of the MagE slit apertures. Alex Navarre used GALFIT to remove the interloper galaxy from both filters, creating cutouts of the original images centered on the two largest arcs of the Sunburst Arc. The GALFIT data products also include the model galaxy profile.

The notebook fetches those cutouts and assembles them into a multi-panel figure. To do this, it starts by fetching the FITS files of the cutout images, and defines a new, smaller, square cutout area to later apply to these images using Astropy's `SkyCoord` object. The new cutout centers on the interloper galaxy. The notebook also creates source masks of the two filters from the cutouts of the original, unaltered images to later apply. After creating the matplotlib `Figure` and `Axes` objects for the figure, the notebook explicitly sets the vertical (`hspace`) and horizontal spacing (`wspace`) between the subplots to 0 to make the figure more compact.

The following double loop iterates over each filter (F275W and F814W) and data type (original image, model, or residual) combination, each time creating an Astropy `Cutout2D` object of the smaller, new cutout area with the defined `SkyCoord` object and the WCS of the given image. From the aforementioned source masks, the notebook calculates the pixel value color boundaries of the new cutout image, and then plots the new cutout image on the appropriate panel.

Because of difficulties encountered while trying to assign the multi-subplot figure a WCS, `plot()` opts to instead manually set the coordinate limits, ticks, and tick labels, utilizing that the drizzled pixel size of the images is 0.03″. The final size of each panel is 3″ on a side, with ticks and tick labels (except if the tick labels would appear between two panels) placed at -1″, 0″, and +1″ relative to the center on both coordinate axes. The notebook sets the aspect ratio of each panel to guarantee it appears square, and adds annotations to each panel indicating the image's filter and data type. Finally, the notebook adds coordinate axis labels to the figure and explicitly sets the tick parameters of each panel to ensure there are no tick labels for ticks facing other panels, and that each side of every panel has inward-facing ticks. The `plot()` function then saves the figure as `foreground_galaxy_galfit_modeling_results.pdf` to `figs/`.

### `lya_and_lyc_maps.ipynb`

<!--

This notebook plots the Lyα and LyC maps of the two largest arcs of the Sunburst Arc in a multi-panel figure.

In the `plot()` function the notebook fetches the common WCS of the HST images and creates two rotated WCSs, such that either arc appears horizontal. Then, for each map to plot, the function gets the image, reprojects it to the rotated WCSs, and plots the rotated image in the corresponding panel. It continues by setting the coordinate limits of the panels and disabling their ticks and labels. The remaining code adds labels indicating the arcs the panels show and the type of map, as well as compasses, scalebars, the footprints of the MagE slit apertures, and labels of the images of source plane clumps. The notebook saves the final figure as `lya_and_lyc_maps.pdf` in the `figs/` folder.

-->

### `lya_nb_m3.ipynb`

This notebook creates a multi-panel figure of the narrowband Lyα maps of the main image in the footprint of the MagE slit aperture M3, dubbed 'Godzilla' in the literature.

Because the two Lyα map schemes are highly discrepant about the Lyα properties of Godzilla, it is important to show the two maps at the location of Godzilla at various statistical significances, which is the purpose of the notebook's sole function `plot()`. The function opens by fetching the `.fits` files of the two Lyα map schemes (with F390W and F555W as the local continuum estimates), and defines the celestial coordinate center of the cutout area to plot in each panel as an Astropy `SkyCoord` object, as well as the width of the square cutout area (slightly wider than the coordinate boundaries set later). After instantiating the matplotlib `Figure` and `Axes` object of the figure, `plot()` explicitly sets the `hspace` and `wspace` of the subplots to 0 to make the panels contiguous.

In the following loop, for the Lyα map of each continuum estimate, the notebook fetches the `HDUList` of the map's `.fits` file and creates a `Cutout2D` object for both the fluxes and uncertainties of the image using the WCS of the map, previous `SkyCoord` object, and the set width of the cutout. A nested loop for a variety of statistical deviations (e.g., +1σ, -2σ, etc.) then plots the cutouts at those statistical deviations in the appropriate panels. In each panel, the nested loop also draws, in pixel coordinates, (1) a manually placed ellipse around Godzilla and the two nearby, smaller clumps, and (2) the footprint of slit M3's aperture, with its pixel coordinates manually calculated. 

Because it is difficult to apply a WCS projection to plots created from `plt.subplots()`, the nested loop also manually calculates the coordinate limits of the panels so that each is 2″ on a side, since each pixel in the drizzled images is 0.03″ across. Following this, `plot()` adds inward-facing ticks to all the plots, and only adds tick labels when the labels would not be between plots, calculating the coordinates of the ticks in similar fashion as the coordinate limits. The nested loop concludes by setting the aspect ratio of each panel to ensure they appear square, and adding labels indicating the statistical deviation of the panel and the continuum estimate method.

`plot()` concludes by adding coordinate labels to the figure and saving the figure to `figs/` as `lya_nb_m3.pdf`.

### `lya.ipynb`

This notebook is the core of the article's results. It measures various properties of the Lyα profiles of the spectra of the MagE slit apertures, calculates statistical correlations between the measurements, and makes associated tables and plots.

Several key functions make up the notebook: `fit()`, `correlate()`, `tabulate()`, and `plot()`, executed in that order. The function `fit()` measures various Lyα profile properties; some by directly fitting a model to the Lyα profile curve, and others with a methodology that is less curve-sensitive. Then, `correlate()` calculates the statistical correlations between every measurement combination (including the LyC escape fractions). `tabulate()` and `plot()` make LaTeX-formatted tables of measurements, and various plots showcasing the data and measurements, respectively.

`fit()` is the most complex function of the notebook. It starts by defining a flat, ΛCDM cosmology with 30% matter energy density and a Hubble constant of 70 km s<sup>-1</sup> Mpc<sup>-1</sup> that compromises between current late- and early-time measurements. This is important because it will be necessary when calculating the Lyα luminosity later in the function. 

Next, `fit()` opens a `for` loop, iterating over each of the MagE slit apertures (including the stacked spectra of the LyC-leaking and non-LyC-leaking apertures). The loop starts by fetching the spectrum's wavelength bins, flux densities, and flux density uncertainties, trimming them to just the relevant portion of the spectrum (roughly 20 Å in the rest frame to either side of the Lyα transition wavelength). This makes the data less bulky and quicker to process. `fit()` then calculates the luminosity distance based on the spectrum's redshift, and creates several masks of the spectrum that will be helpful later. 

Before initiating the Monte Carlo simulation that makes the various Lyα measurements, `fit()`, based on the current spectrum of the first `for` loop, specifies the function to fit to the Lyα profile and its initial parameter estimates and parameter boundaries. For all the spectra, the model function is a combination of Gaussian and skewed Gaussian functions. The notebook then executes a burn-in MC simulation to improve the initial parameter estimates. The burn-in (and the consecutive, main) MC simulation randomly perturbs the rest-frame spectrum by assuming the flux densities and flux density uncertainties of the spectrum correspond to the mean and standard deviations of Gaussian distributions at each wavelength bin, randomly drawing from each distribution to obtain a new spectrum. After completing the burn-in MC simulation, `fit()` saves the randomly sampled spectra and best-fit model parameters of the simulation to `{slit_id}_mc_sim_burn_in_lya_spectra.txt` and `{slit_id}_mc_sim_burn_in_lya_best_fit_model_parameters.txt`, respectively, in the `results/lya_fits/{slit_id}/` folder.

`fit()` then initiates the main MC simulation in a nested `for` loop, calculating the rest-frame Lyα equivalent width, central escape fraction of Lyα (see [Naidu & Matthee et al. (2022)](https://doi.org/10.1093/mnras/stab3601) (MNRAS, 510, 4582)), Lyα luminosity, best-fit parameters of the model function fit to the Lyα profile, as well as the ratio between the 'minimum' flux density between the redshifted and blueshifted Lyα peaks and the local continuum (often called $f_{\rm{min}}/f_{\rm{cont}}$; see the article for more details about how this work treats this calculation). After the main MC simulation concludes, the notebook calculates the FWHM of each Lyα peak and, if applicable, the velocity separation between the redshifted and blueshifted Lyα peaks by applying the best-fit model parameters found in the MC simulation and the `compute_fwhm_and_peak()` function. `fit()` then concludes by saving the randomly sampled spectra of the MC simulation to a file named like `{slit_id}_mc_lya_best_fit_model_parameters.txt` in `results/lya_fits/{slit_id}/` and saving the Lyα measurements (all the measurements that weren't the best-fit model parameters) to a file named like `{slit_id}_mc_sim_lya_measurements.txt` in `results/lya_fits/{slit_id}/`.

An important goal of the article is to understand&mdash;particularly on a subgalactic scale&mdash;how Lyα and LyC properties relate, as well as how different Lyα properties relate to each other. To help judge this, the article presents statistical correlations between all of the Lyα and LyC properties, which is the primary purpose of `correlate()`. `correlate()` opens by fetching the number of iterations of the main MC simulation and the LyC escape fraction measurements made by `esc.ipynb()`. Then, using a series of `for` loops, `correlate()` calculates the Pearson and type 'b' Kendall correlation coefficients for each unique measurement combination and each iteration, saving the results to the file `mc_sim_lya_measurements_statistical_correlations.txt` in `results/lya_fits/`.

After `fit()` and `correlate()` execute, `lya.ipynb` has made all of its measurements. `tabulate` arranges most of the sets of measurements into LaTeX-formatted `.txt` tables that the article's main `.tex` file can automatically incorporate without manually typesetting the table in the `.tex` file. `tabulate()` starts by creating the header of the Lyα measurements table with a manually created string. Then, using a `for` loop iterating over each MagE slit aperture spectrum (including the stacked spectra), `tabulate()` fetches the Lyα measurements of each slit. The notebook deliberately removes the spectral resolutions randomly sampled for the slit from the measurements, since `fit()` just kept the sampled spectral resolutions for reproduceability because they corrected the measured peak FWHMs for the instrumental resolution FWHM. Following this, the loop adds the slit ID to the table row and opens a nested `for` loop which iterates over the type of Lyα measurement. For each type of measurement, the loop calculates the median of the measurement distribution and the differences of the 16th and 84th percentiles of the measurement distribution compared to the median (where the latter two are a proxy for the upper and lower uncertainty bounds). If the measurement is the Lyα luminosity, the notebook rescales the measurement distribution to eliminate a cumbersome scientific notation factor that would otherwise make tabulating the values more difficult. `tabulate()` accounts for this rescaling later to ensure the presented values are accurate. If the median value of the measurement distribution is a NaN, the loop will simply add a new column to the table row and represent the measurement value as a '&mdash;', indicating that the measurement type is not applicable for the spectrum (e.g., the Lyα luminosity of the stacked spectra, etc.). Otherwise, the loop will apply `round_to_uncertainties()` to the quantities calculated from the measurment distribution, rounding the upper and lower uncertainty estimates to 1 significant figure, and rounding the median measurement value to the most significant digit of the two uncertainties. If the measurement type is the central Lyα peak's FWHM, the loop checks if the median measurement is smaller than the instrumental resolution FWHM, in which case it just states the 84th percentile of the measurement distribution as an upper bound. In either case, it then formats and adds the median and uncertainties to a new column in the table row.

After iterating through all the measurement types for a single spectrum, `tabulate()` adds a LaTeX row break character if appropriate and a line break character in the string. Following iterating through all the spectra, `tabulate()` adds a footer to the table and saves the table as `lya_measurements_table.txt` in `results/tables/`.

Next, `tabulate()` creates the table containing the measured statistical correlations and uncertainties between each pair of measurements, including the LyC escape fraction. After creating the manually wrote table header and fetching the correlation measurements, a `for` loop iterates through each measurement pair, adding an appropriate label to the table row and extracting both sets of statistical correlation measurements of the measurement pair. As in the previous table, `tabulate()` calculates the median and estimated uncertainties of the two statistical correlation distributions, feeding the quantities into `round_to_uncertainties()` and formatting the resulting strings into a new column in the table row. The loop closes by adding a LaTeX column break to the table row if appropriate and a line break to the string. The notebook adds a footer to the table and saves it as `lya_measurements_statistical_correlations_table.txt` in `results/tables/`.

The last table `tabulate()` constructs lists the best-fit model parameters of the functions fit to the Lyα profiles. Starting by creating the header of the table, `tabulate()` then opens a `for` loop that iterates through each MagE spectrum (including the stacked spectra), adding the slit ID to the table row and fetching the spectrum's best-fit model parameters from `{slit_id}_mc_sim_lya_best_fit_model_parameters.txt` from `results/lya_fits/{slit_id}/`. The loop checks if the spectrum's Lyα profile is triple-peaked, adding a new column with a blank entry if the profile is not triple-peaked. Otherwise, for each model parameter of the blueshifted Lyα peak, the loop calculates the median and estimates of the upper and lower uncertainties (calculated as stated previously) from the parameter distribution, rounding them according to `round_to_uncertainties()`, and formats and adds them to the table. The loop then adds a new column and repeats the previous process for the remaining model parameters (i.e., the remaining Lyα peaks and the local continuum). If the parameter has units of a flux density, the notebook rescales the parameter by a fixed and consistent factor (except for the non-stacked spectra) in order to prevent a cumbersome scientific notation factor that would otherwise complicate formatting the parameter as a string. The loop then closes and the notebook adds a footer to the table and saves the table as `lya_best_fit_model_parameters_table.txt` in `results/tables/`.

The function `plot()` is the final function executed by `lya.ipynb`, creating plots of the Lyα profiles and the corner plot of the Lyα and LyC escape fraction measurements. `plot()` starts by obtaining the LyC escape fraction measurements made by `esc.ipynb` and creating the matplotlib `Figure` and `Axes` objects of the figure plotting the Lyα profiles of the individual MagE spectra (but not the stacked spectra). Using a `for` loop, `plot()` gets the ...

...

### `map.ipynb`

This notebook creates the primary on-sky map figure of the article, showing the Sunburst Arc and its cluster field, the MagE slit apertures, and zoom-in insets on the two largest arcs (those targeted by the MagE observations) with labeled clump numbers.

The notebook's primary function `plot()` starts by fetching the data and WCS of the HST WFC3/F606W image (also the common WCS of all the images) of the Sunburst Arc and its cluster field. `plot()` then instantiates a figure with the obtained WCS of the HST images and plots the F606W image of the field. It then sets the coordinate boundaries and labels of the figure. Using `rotate_wcs()`, it creates the rotated WCSs of the inset panels by rotating the WCS of the main panel so that the arcs in the inset panels will appear vertical or horizontal. Then, with `add_inset_axes()`, the notebook creates the two inset axes at the given locations and with the rotated WCSs. After reprojecting the main panel's image to the rotated WCSs, `plot()` adds the rotated images to the inset axes, and establishes their coordinate limits with `set_inset_axes_limits()`. 

The notebook then disables the ticks and labels of the inset axes with a double `for` loop, and executes `add_box()` several times to add footprints of the displayed area in the inset axes on the main image, as well as highlight the other arc segments of the Sunburst Arc. For each of these boxes in the main image and the corresponding inset axes, the notebook labels them according to the terminology of the lens model of [Sharon et al. (2022)](https://doi.org/10.3847/1538-4357/ac927a) (ApJ, 941, 203). For each of the MagE slit apertures, `plot()` instantiates the aperture as a regions `RectangleSkyRegion()` object, converts it to a pixel-based representation given the WCS of each panel the aperture appears in, and then plots and labels the footprint with its slit ID. For each image of a source plane clump that appears in the inset axes, `plot()` circles (or boxes, for one image) and labels the clump. Following this, `plot()` concludes by annotating the main axes with a compass, adds scalebars to the main and inset axes, and finally saves the figure as `sunburst_arc_and_cluster_field.pdf` to `figs/`.

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

The `stack.ipynb` notebook creates stacked spectra of the MagE slit apertures targeting the LyC-leaking region and the non-LyC-leaking region. 

</div>