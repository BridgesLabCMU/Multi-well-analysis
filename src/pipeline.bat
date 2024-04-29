@echo off

set "batch_path=%~dp0"
set "batch_path=%batch_path:~0,-1%"

if not exist "%batch_path%\experiment_config.json" (
    python3 .\GUI\gui.py
    if %errorlevel% neq 0 goto end
)

julia ExtractDataV2.jl
if %errorlevel% neq 0 goto end

for /f "delims=" %%i in ('powershell -Command "Get-Content '%batch_path%\experiment_config.json' | ConvertFrom-Json | Select -ExpandProperty 'image_analysis'"') do set "image_analysis=%%i"
if "%image_analysis%"=="True" (
	echo Performing image analysis 
	julia CytationAnalysisV2.jl
	if %errorlevel% neq 0 goto end
)

wsl julia "/mnt/c/Users/Imaging Controller/Desktop/Multi-well-analysis/src/CytationPlottingV2.jl"
if %errorlevel% neq 0 goto end

setlocal enabledelayedexpansion
set "batch_path=%~dp0"
set "batch_path=%batch_path:~0,-1%"

for /f "delims=" %%i in ('powershell -Command "Get-Content '%batch_path%\experiment_config.json' | ConvertFrom-Json | Select -ExpandProperty experiment_directory"') do set "experiment_directory=%%i"
for /f "delims=" %%d in ('powershell -Command "Get-Date -Format yyyyMMdd"') do set "date_folder=%%d"

set /a count=0
set "new_folder=%experiment_directory%\%date_folder%"
:loop
if exist "!new_folder!\" (
    set /a count+=1
    set "new_folder=%experiment_directory%\%date_folder%_!count!"
    goto loop
)

mkdir "!new_folder!"

for /f "delims=" %%a in ('powershell -Command "Get-Content '%batch_path%\experiment_config.json' | ConvertFrom-Json | Select -ExpandProperty bulk_data | ForEach-Object { $_ -join '`n' }"') do (
    set "file_path=%%~a"
	set "dirPath=!file_path:~0,-4!"
    if exist "!file_path!" (
		set "file_path=!file_path:/=\!"  
        echo Moving file: "!file_path!"
        move "!file_path!" "!new_folder!"
    )
    if exist "!dirPath!" (
		set "dirPath=!dirPath:/=\!"  
		for %%n in ("!dirPath!") do set "dir_name=%%~nxn"
        echo Moving folder: "!dirPath!"
		%systemroot%\System32\robocopy "!dirPath!" "!new_folder!\!dir_name!" /E /MOVE
    )
)

for /f "delims=" %%g in ('powershell -Command "Get-Content '%batch_path%\experiment_config.json' | ConvertFrom-Json | Select -ExpandProperty 'good_data_directory'"') do set "good_data=%%g"
echo good_data: !good_data!
for /f "delims=" %%i in ('powershell -Command "Get-Content '%batch_path%\experiment_config.json' | ConvertFrom-Json | Select -ExpandProperty images_directory | ForEach-Object { $_ -join '`n' }"') do (
    set "dir_path=%%~i"
    :: Get the immediate parent directory path
    for %%p in ("!dir_path!\..") do set "parent_path=%%~fp"
    :: Get the name of the parent directory to copy
    for %%n in ("!parent_path!") do set "parent_name=%%~nxn"
    
    :: Copy experiment_config.json into the parent directory
    if exist "!parent_path!" (
        copy "%batch_path%\experiment_config.json" "!parent_path!"
    )
    
    if exist "!parent_path!" (
        echo Copying parent directory: "!parent_path!"
        %systemroot%\System32\xcopy "!parent_path!" "!new_folder!\!parent_name!" /E /I
        :: Copy to "Good data" directory if not empty
        if not "!good_data!"=="" (
            %systemroot%\System32\xcopy "!parent_path!" "!good_data!\!parent_name!" /E /I
        )
    )
)

move "%batch_path%\experiment_config.json" "!new_folder!"
echo "%batch_path%\experiment_config.json"

echo Operation completed.
pause

endlocal

:end
