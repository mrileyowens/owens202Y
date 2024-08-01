This file details the structure of this repository and the nature of its files.

## Data

The software of this repository uses data from two primary sources: the Magellan Echellette (MagE) spectrograph mounted on ... 

## Repository structure

This repository contains three main folders: `docs`, `figs`, and `results`. The `docs` folder contains the project documentation, such as this file and the project's `readme.md`, the `figs` folder contains the figures used in the article, and the `results` folder contains various measurements and formatted tables used in the article.

Outside these folders are the `.ipynb` Jupyter notebooks used to create the figures, measurements, and tables contained in `figs` and `results`.

## Repository files

### `esc.ipynb`

This notebook is responsible for calculating the LyC escape fractions and tabulating the measurements.

The primary function of the notebook is `measure()`, which is the function that actually measures the LyC escape fractions of the MagE slit apertures.

### `masks.ipynb`

This notebook makes masks of the MagE slit apertures and the two largest arcs of the Sunburst Arc in the world coordinate system (WCS) of the latest reduction of the HST images used in the article. These masks are critical for portions of analyses in other notebooks.

The notebook starts by fetching the WCS of one of the HST images of the latest reduction (all the images of the latest reduction share a common WCS, so the choice of image is arbitrary). Then, for each MagE slit aperture, the notebook creates a regions `RectangleSkyRegion()` object of the aperture (since the apertures are rectangular) based on its center, width, and rotation on the sky from a dictionary in the notebook, and converts it to a pixel representation from the HST image's WCS. From the pixel representation, the notebook creates a mask of the aperture matching the shape of the HST image, and creates a new FITS header from the WCS. The notebook then adds several relevant keywords to the header, and saves the header and mask to a new `.fits` file, following the name format `{slit_id}_mask_v5.fits`. The notebook saves the new file to `results/masks/`.