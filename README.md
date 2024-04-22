# Multi-well-analysis
Scripts for analyzing imaging and bulk population data from multi-well experiments. Julia code works as a standalone, regardless of OS. Python code will run on a Windows OS, as does the executable pipeline.bat. 

## Notes for using the GUI
1. Images directories and bulk data file paths can be specified multiply if the experiment involves more than one plate. 
2. Condition names should look like you would want them to look in a plot. Relatedly, italics can be used by writing the $ symbol around the characters to be specified (as in LaTeX). Similarly, Greek letters can be specified by writing, e.g., "$\\mu$" = $\\mu$ (without quotes).
3. The order of selection for plot data types, etc., should be consistent during plot settings specification. For example, if you are using the "two-axis" plot option, you should select the plot data type corresponding to the left y-axis first, the data type corresponding to the right y-axis second, likewise for their respective normalization methods, etc.
4. GUI widgets need to be manually closed before analysis can proceed.
5. B drive should be connected prior to running the GUI.
6. Specify conditions for all wells in the experiment, even if you do not plan on plotting data from certain wells.
7. Y-limits can be manually specified (for now, only for jitter plots) like the following: 0,3 (lower and upper) or 0, (just lower) or ,3 (just upper)

## To-do's
1. Allow more user specifications for plot aesthetics (colors, fonts, ...)
2. Add a separate statistics script
