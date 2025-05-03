# Multi-well-analysis
Scripts for analyzing brightfield imaging data from multi-well experiments. The executable BiofilmAssay.bat, located in the main directory, only runs on Windows. To run the code on Unix/Linux, you can just cd into src and run julia BiofilmAssay.jl (which is all the .bat file does anyway, after the first run). 

## Installation and basic use
If you want to use this code for the first time on a new computer, first install the Julia programming language from this site: https://julialang.org/install/.

Then, download this repository, preferably via git, although downloading it as a zip is fine. If running code through the batch file, this is all you have to do. The first time you run the batch file, it will automatically install necessary packages etc. which will cause a long lag (only the first time). If you are running the code directly from src, you will have to do the following in command prompt/terminal:

```
juliaup add 1.11.2
cd /path/to/Multi-well-analysis/src
julia +1.11.2 --project=. -e 'using Pkg; Pkg.instantiate()'
```

Thereafter, you can run the code by `cd /path/to/Multi-well-analysis/src`, then run `julia +1.11.2 --project=. BiofilmAssay.jl`. 

## Notes
Guidelines for image analysis parameters (set in the first GUI):
1. Block diameter: A good guideline would be the linear dimension of the largest object to be segmented (in number of pixels). Preset to 101. I have found that non-intuitive segmentations can occur beyond block diameter ~300 pixels, so increase with caution. Of course you can always guess and check.
2. Threshold: I know this is highly microscope dependent (coincidentally it seems like for the two Cytations, the same threshold can be used...don't know why). For the Cytations, a value of 0.04 is the upper bound that I would try and 0.02 is the lower bound I would try. Smaller values are going to be more "liberal" in the assignment of pixels to biofilm.
