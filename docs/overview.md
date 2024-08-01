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