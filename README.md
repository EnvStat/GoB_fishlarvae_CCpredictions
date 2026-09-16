# GoB_fishlarvae_CCpredictions

This repository includes the necessary files to reproduce the results from the following paper: Pia, I., Veneranta, L., Siiriä, S.-M., and Vanhatalo, J. (2026). Effects of Climate Change on Reproduction of Sea-spawning Coregonids in the Baltic Sea. Marine Ecology Progress Series, in press.

Instructions for use:

1) Run code in: 'Data_preparation.R' 
  To read the data (available as .txt files in the 'Data' folder) and save them as a single R object in the 'Data' folder.

2) Run code in: 'Models_run.R'
  To get posterior samples from both whitefish and vendace integrated SDMs.
  The Stan samples will be automatically saved in the 'Stan_models' folder.
  The Stan code for both integrated SDMs is also stored in the 'Stan_models' folder.

3) Run code in: 'Outputs.R'
  To reproduce all figures present in the main manuscript and in the supplement.
  To reproduce a specific figure, after running lines 1-47 to load posterior samples, functions, and data, navigate the R menu to the desired figure.
  NB: for some figures there is a suggested 'img.ratio' for the RStudio Plots window, which will allow for proper positioning of extra labels and other graphical elements.
