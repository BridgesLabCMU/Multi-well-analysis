#!/bin/bash

python3 ./GUI/gui.py && julia ExtractData.jl && julia CytationAnalysis.jl && julia CytationPlotting.jl
