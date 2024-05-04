# Multi-well-analysis
Scripts for analyzing imaging and bulk population data from multi-well experiments. Julia code works as a standalone, regardless of OS. Python code will run on a Windows OS, as does the executable pipeline.bat. 

## Notes for using the GUI
1. Images directories and bulk data file paths can be specified multiply if the experiment involves more than one plate. 
2. Condition names should look like you would want them to look in a plot. Relatedly, italics can be used by writing the $ symbol around the characters to be specified (as in LaTeX). Similarly, Greek letters can be specified by writing, e.g., "$\\mu$" = $\\mu$ (without quotes). If you want a greek symbol immediately followed by italics, include a space: "$\Delta vpsL$". Subscripts and superscripts need to also be included within $$. For instance, "$dbfR^{\mathrm{D51V}}$". In that example,
   everything within brackets is a superscript, and "\mathrm" makes the text that it brackets non-italicized.
3. The order of selection for plot data types, etc., should be consistent during plot settings specification. For example, if you are using the "two-axis" plot option, you should select the plot data type corresponding to the left y-axis first, the data type corresponding to the right y-axis second, likewise for their respective normalization methods, etc.
4. GUI widgets need to be manually closed before analysis can proceed.
5. B drive should be connected prior to running the GUI.
6. Specify conditions for all wells in the experiment, even if you do not plan on plotting data from certain wells.
7. Y-limits can be manually specified (for now, only for jitter plots) like the following: 0,3 (lower and upper) or 0, (just lower) or ,3 (just upper)
8. For grouped boxplots: ensure that for each condition, characters before the first space (not within $$) are common to a given group (e.g., "WT +01uM": "WT" is identified as a group). The part that refers to the group name should also be written as it would appear in a plot (e.g., "$\Delta vpsL$ +01uM")  
9. There are defaults set for the color, font, and plot size entries. The default plot size is (3,2.5), which would be a compact line plot. If you know your plot is large horizontally, increase the first number accordingly.
10. The color can be the name of a preferred colorscheme, a single hex color, or a list of hex colors. The default for a jitter/box plot is a single blue-ish color and the default for a two-axis plot is a blue-ish and a salmon color,  while for a heatmap, line plot, or grouped jitter plot, the defaults are colorschemes. You can access useable colorschemes here: https://matplotlib.org/stable/users/explain/colors/colormaps.html.

## To-do's
1. Add a separate statistics script
