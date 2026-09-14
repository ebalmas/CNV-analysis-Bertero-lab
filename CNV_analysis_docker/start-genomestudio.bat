@echo off
REM ============================================================================
REM  GenomeStudio in JupyterLab - standalone launcher (Windows)
REM
REM  Double-click and go. No PowerShell, no ExecutionPolicy, no code signing:
REM  .bat files are not subject to those checks.
REM
REM  The container runs in the FOREGROUND: this window shows the JupyterLab
REM  log and stays open. Press Ctrl-C to shut it down (twice if it does not
REM  stop immediately). Closing the window also stops it.
REM
REM  The folder containing this file is shared with the container as "shared".
REM  Inside GenomeStudio you will find it as Z:\home\jovyan\shared
REM ============================================================================

setlocal enabledelayedexpansion

REM --- configuration ----------------------------------------------------------
set "IMAGE=repbioinfo/genomestudio:1.3"
set "CONTAINER=genomestudio"
set "PORT=8888"
set "MEMORY=12g"
set "URL=http://localhost:%PORT%/lab"

REM %~dp0 ends with a backslash, which Docker will not accept in -v
set "SHARED_DIR=%~dp0"
set "SHARED_DIR=%SHARED_DIR:~0,-1%"

REM --- is Docker available? ---------------------------------------------------
where docker >nul 2>&1
if errorlevel 1 (
    echo !! docker not found. Install Docker Desktop for Windows.
    goto :end
)
docker info >nul 2>&1
if errorlevel 1 (
    echo !! Docker engine is not responding.
    echo    Start Docker Desktop, wait until the icon turns green, then retry.
    goto :end
)

REM --- image present? if not, pull it -----------------------------------------
docker image inspect %IMAGE% >nul 2>&1
if errorlevel 1 (
    echo ==^> image not found locally, downloading
    echo     %IMAGE%
    echo     This is several GB. First run takes a while - do not close.
    echo.
    docker pull %IMAGE%
    if errorlevel 1 (
        echo.
        echo !! Download failed. Three likely causes:
        echo    - no connection, or a proxy/firewall in the way
        echo    - image is not on a public registry: run  docker login  first
        echo    - wrong tag: fix the IMAGE= line at the top of this file
        goto :end
    )
    echo.
    echo ==^> download complete
)

REM --- already running? -------------------------------------------------------
docker ps --format "{{.Names}}" | findstr /x "%CONTAINER%" >nul 2>&1
if not errorlevel 1 (
    echo ==^> already running in another window
    echo     open %URL%
    echo     To stop it: press Ctrl-C in that window, or run
    echo       docker stop %CONTAINER%
    goto :end
)
docker rm -f %CONTAINER% >nul 2>&1

REM --- port check -------------------------------------------------------------
REM  Fixed port on purpose: the URL must always be the same, with nothing to
REM  copy or paste. So if 8888 is taken we stop and say so, rather than
REM  silently moving to another port and changing the address.
netstat -an | findstr /c:":%PORT% " | findstr /i listening >nul 2>&1
if not errorlevel 1 (
    echo !! Port %PORT% is already in use by another program.
    echo    Close whatever is using it, then retry. To find out what it is:
    echo      netstat -ano ^| findstr :%PORT%
    goto :end
)

REM --- is the Wine prefix already baked into the image? -----------------------
REM  If it is (image committed after a run), do NOT mount anything on
REM  /home/jovyan: the container starts instantly. Mounting a host folder
REM  there would hide the committed prefix and force a 3 GB copy.
echo ==^> inspecting image
set "HOME_ARGS="
docker run --rm %IMAGE% test -d /home/jovyan/.wine-genomestudio >nul 2>&1
if errorlevel 1 (
    echo     prefix not in image: using a local home ^(one-time 3 GB copy^)
    echo     First run: the desktop appears right away, GenomeStudio later.
    set "HOME_DIR=%SHARED_DIR%\.genomestudio-home"
    if not exist "!HOME_DIR!" mkdir "!HOME_DIR!"
    set "HOME_ARGS=-v "!HOME_DIR!":/home/jovyan"
) else (
    echo     prefix already in image: direct start, no copy
)

REM --- browser opener ---------------------------------------------------------
REM  Runs alongside the container: polls the port and opens the browser once
REM  JupyterLab answers. It has to be a separate process because the docker
REM  run below occupies this window.
set "OPENER=%TEMP%\gs-open-browser.bat"
>  "%OPENER%" echo @echo off
>> "%OPENER%" echo for /l %%%%i in (1,1,120) do (
>> "%OPENER%" echo   timeout /t 1 /nobreak ^>nul
>> "%OPENER%" echo   netstat -an ^| findstr /c:":%PORT% " ^| findstr /i listening ^>nul 2^>^&1
>> "%OPENER%" echo   if not errorlevel 1 (
>> "%OPENER%" echo     start "" "%URL%"
>> "%OPENER%" echo     exit /b
>> "%OPENER%" echo   )
>> "%OPENER%" echo )
start /b "" cmd /c "%OPENER%"

REM --- start in the foreground ------------------------------------------------
REM  No -d: this window stays attached, shows the JupyterLab log, and Ctrl-C
REM  is forwarded to the container. --rm cleans up on exit.
REM
REM  seccomp=unconfined : without it dbus fails and the XFCE desktop stays grey
REM  shm-size=2g        : Chromium crashes with Docker's default 64 MB
REM
REM  Auth is disabled through a config file rather than a command-line flag.
REM  "--IdentityProvider.token=" with nothing after the equals sign is silently
REM  ignored and Jupyter falls back to a generated token; a quoted empty string
REM  does not survive cmd plus Docker quoting. In a config file the empty
REM  string is written as str(), so no quotes are needed at all.
echo.
echo ==^> starting %IMAGE%
echo     %URL%
echo.
echo     In the Launcher click the Desktop icon: GenomeStudio starts by itself.
echo     Files in this folder appear as Z:\home\jovyan\shared
echo.
echo     Press Ctrl-C here to stop the container.
echo.

docker run --rm -it ^
    --name %CONTAINER% ^
    !HOME_ARGS! ^
    --security-opt seccomp=unconfined ^
    --shm-size=2g ^
    --memory=%MEMORY% ^
    -p 127.0.0.1:%PORT%:8888 ^
    -v "%SHARED_DIR%":/home/jovyan/shared ^
    -e HOME=/home/jovyan ^
    -e GS_AUTOSTART=1 ^
    %IMAGE% ^
    bash -c "printf 'c.IdentityProvider.token = str()\nc.ServerApp.password = str()\n' > /tmp/jupyter_noauth.py; exec /home/python4Jup/bin/jupyter lab --config=/tmp/jupyter_noauth.py --ip=0.0.0.0 --port=8888 --no-browser --allow-root --ServerApp.root_dir=/home/jovyan"

REM --- cleanup ----------------------------------------------------------------
del "%OPENER%" >nul 2>&1
docker rm -f %CONTAINER% >nul 2>&1
echo.
echo ==^> stopped

:end
echo.
pause
