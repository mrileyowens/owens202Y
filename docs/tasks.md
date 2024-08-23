### Open

- [ ] Clean up all files by removing obsolete code and adding more and any missing comments
- [ ] Double check all file headers and code comments for accuracy and consistency
- [ ] Add all Python packages and verisons to the software credit in the article
- [ ] Determine the units of the Lyα maps and verify the units stated in the documentation and code comments 
- [ ] Update the article's methodology with new details about the MC simulation process
- [ ] Add the descriptions of `nb.ipynb` and `seeing.ipynb` once the results of `seeing.ipynb` are satisfactory

### Closed

- [x] Consider whether it is more appropriate to report fluxes (rather than flux densities) in the seeing simulation (closed 22 August 2024)
- [x] Add description of `limit.ipynb` to `overview.md` (closed 22 August 2024 by #19)
- [x] Update flow chart to reflect that the `galfit.ipynb` output also feeds into the final article (closed 22 August 2024 by #18)
- [x] Add `limit.ipynb` to the flow chart (closed 22 August 2024 by #18)
- [x] Add the calculated limiting magnitude to the article (closed 22 August 2024)
- [x] Consider how to determine the number of significant figures to report for the depth of the F275W image (closed 22 August 2024)
- [x] Optimize the convolutions in `esc.ipynb` and `seeing.ipynb` so that they only calculate the effect of the convolution on the pixels summed for photometric measurements (closed 20 August 2024 by #15)
- [x] Create a virtual environment with the latest Python and package versions and execute the project's code with the environment (closed 15 August 2024 by #14)
- [x] Create a `requirements.txt` file listing the package versions the project uses (closed 15 August 2024 by #14)
- [x] Add horizontal lines separating the LyC-leaking and non-LyC-leaking MagE slit apertures in tables (closed 13 August 2024 by #13)
- [x] Update the convolution in `esc.ipynb` to also include the effect of airmass (closed 13 August 2024 by #12)
- [x] Confirm that all uses of the stacked spectra use the correct and current version produced by `stack.ipynb` (closed 12 August 2024 by #11)
- [x] Update all files to use the current naming scheme of the MagE spectra (closed 11 August 2024 by #10)