# Visualization of coverage.

  * Use script ```Coverage_plot.V3.20210706.R```. It generates a set of plots based on peaks proided, so later they can be explored, and selected for the ones showing the biggest difference between cell groups.

## Notes

1. We use +/- 2Kb outside peak to show both peak and its surroundings.


2. Need to set ```Idents``` to tell the function what cell idents it should use for vizualization. ```Idents``` is a factor, so it can be used for ordering of rows with coverage tracks per ident.


3. We also use annotation table to add to the file name respective gene symbol.


  * Also contains some useful code for combining tracks.

  * Check also here: ```SC_analysis_scripts/6.Coverage_plots/README.md```.