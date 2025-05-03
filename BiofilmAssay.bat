@echo off
REM -------------------------------
REM 1) Go to the folder this script lives in
juliaup add 1.11.2
pushd %~dp0

REM 2) Activate the project & install any missing deps
julia +1.11.2 --project=src -e "using Pkg; Pkg.instantiate()"

REM 3) Run your assay
cd src
julia +1.11.2 --project BiofilmAssay.jl

REM 4) Pause so you can read any output
pause
