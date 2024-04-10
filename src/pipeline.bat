@echo off

python3 .\GUI\gui.py
if %errorlevel% neq 0 goto end

julia ExtractData.jl
if %errorlevel% neq 0 goto end

julia CytationAnalysis.jl
if %errorlevel% neq 0 goto end

julia CytationPlotting.jl

:end
