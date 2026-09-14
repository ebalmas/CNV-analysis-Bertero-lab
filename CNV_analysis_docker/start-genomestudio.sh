#!/usr/bin/env bash
# =============================================================================
# GenomeStudio in JupyterLab - standalone launcher (Linux / macOS)
#
# Requires neither JupyterHub nor JupyDo: it runs the container on its own
# and opens JupyterLab in the browser.
#
# The container runs in the FOREGROUND: the terminal shows the JupyterLab log
# and stays open. Press Ctrl-C to shut it down (twice if it does not stop
# immediately).
#
# Usage:
#   ./start-genomestudio.sh              start
#   ./start-genomestudio.sh --reset-home discard the local home and start over
#   ./start-genomestudio.sh --port 9000  use a different port
#   ./start-genomestudio.sh --image X    use a different image tag
#
# The folder you launch this script from is shared with the container as
# "shared": inside GenomeStudio it appears as Z:\home\jovyan\shared
#
# ON THE HOME DIRECTORY: if the image already carries the Wine prefix at
# /home/jovyan/.wine-genomestudio (because it was committed after a run),
# nothing is mounted on /home/jovyan and the container starts instantly.
# Mounting a host folder there would hide the committed prefix and force the
# wrapper to re-copy 3 GB from /opt/wine-gs. A local home is used only when
# genuinely needed: prefix missing from the image, or a uid other than 1000
# on Linux.
# =============================================================================
set -euo pipefail

# --- configuration -----------------------------------------------------------
IMAGE="${GS_IMAGE:-repbioinfo/genomestudio:1.3}"
CONTAINER="${GS_CONTAINER:-genomestudio}"
PORT="${GS_PORT:-8888}"
MEM="${GS_MEM:-12g}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED_DIR="$SCRIPT_DIR"
HOME_DIR="$SCRIPT_DIR/.genomestudio-home"

# --- options -----------------------------------------------------------------
RESET_HOME=0
while [ $# -gt 0 ]; do
    case "$1" in
        --reset-home) RESET_HOME=1 ;;
        --port)       PORT="$2"; shift ;;
        --image)      IMAGE="$2"; shift ;;
        -h|--help)    sed -n '2,28p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

URL="http://localhost:${PORT}/lab"

say() { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
err() { printf '\033[1;31m!!\033[0m %s\n' "$*" >&2; }

port_busy() {
    if command -v nc >/dev/null 2>&1; then
        nc -z 127.0.0.1 "$1" >/dev/null 2>&1
    else
        (exec 3<>"/dev/tcp/127.0.0.1/$1") 2>/dev/null
    fi
}

# --- preflight checks --------------------------------------------------------
if ! command -v docker >/dev/null 2>&1; then
    err "docker not found in PATH. Install Docker Desktop (macOS) or docker.io (Linux)."
    exit 1
fi
if ! docker info >/dev/null 2>&1; then
    err "Docker daemon is not responding."
    echo "   macOS: start Docker Desktop.  Linux: sudo systemctl start docker"
    exit 1
fi

# --- image present? if not, pull it ------------------------------------------
if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
    say "image not found locally, downloading: $IMAGE"
    echo "    This is several GB. First run takes a while."
    if ! docker pull "$IMAGE"; then
        err "download failed. Likely causes:"
        echo "   - no connection, or a proxy in the way"
        echo "   - image not on a public registry: run 'docker login' first"
        echo "   - wrong tag: pass --image your/tag, or edit IMAGE at the top"
        echo "   Local images available:"
        docker images --format '     {{.Repository}}:{{.Tag}}' | grep -i genomestudio || echo "     (none)"
        exit 1
    fi
    say "download complete"
fi

# --- already running? --------------------------------------------------------
if docker ps --format '{{.Names}}' | grep -qx "$CONTAINER"; then
    say "already running in another terminal"
    echo "    open $URL"
    echo "    To stop it: press Ctrl-C there, or run  docker stop $CONTAINER"
    exit 0
fi
docker rm -f "$CONTAINER" >/dev/null 2>&1 || true

# --- port check --------------------------------------------------------------
# Fixed port on purpose: the URL must always be the same, with nothing to copy
# or paste. So if it is taken we stop and say so, rather than silently moving
# to another port and changing the address.
if port_busy "$PORT"; then
    err "port $PORT is already in use by another program."
    echo "   Free it and retry, or pass --port <n>. To see what holds it:"
    echo "     Linux: ss -lptn 'sport = :$PORT'"
    echo "     macOS: lsof -nP -iTCP:$PORT -sTCP:LISTEN"
    exit 1
fi

# --- is the Wine prefix already baked into the image? ------------------------
say "inspecting image"
HAS_PREFIX=no
docker run --rm "$IMAGE" test -d /home/jovyan/.wine-genomestudio >/dev/null 2>&1 && HAS_PREFIX=yes

# On Linux the container runs as jovyan (uid 1000): with a different uid, the
# files GenomeStudio writes into "shared" would end up owned by someone else.
# Not an issue on macOS, where Docker Desktop translates bind-mount ownership.
NEED_USER_MAP=0
if [ "$(uname -s)" = "Linux" ] && [ "$(id -u)" != "1000" ]; then
    NEED_USER_MAP=1
fi

USE_LOCAL_HOME=0
if [ "$HAS_PREFIX" != "yes" ]; then
    USE_LOCAL_HOME=1
    say "prefix not in image: using a local home (one-time 3 GB copy)"
    echo "    First run: the desktop appears right away, GenomeStudio later."
elif [ "$NEED_USER_MAP" = "1" ]; then
    USE_LOCAL_HOME=1
    say "your uid is $(id -u), not 1000: using a local home to keep file ownership correct"
else
    say "prefix already in image: direct start, no copy"
fi

HOME_ARGS=()
if [ "$USE_LOCAL_HOME" = "1" ]; then
    if [ "$RESET_HOME" = "1" ] && [ -d "$HOME_DIR" ]; then
        say "removing local home $HOME_DIR"
        rm -rf "$HOME_DIR"
    fi
    mkdir -p "$HOME_DIR"
    HOME_ARGS+=(-v "$HOME_DIR:/home/jovyan")
    [ "$NEED_USER_MAP" = "1" ] && HOME_ARGS+=(--user "$(id -u):$(id -g)")
elif [ "$RESET_HOME" = "1" ]; then
    err "--reset-home has no effect: the prefix lives in the image, not in a local folder."
    err "To regenerate it, run inside the desktop: genomestudio --reset"
fi

# --- browser opener ----------------------------------------------------------
# Runs alongside the container: polls the port and opens the browser once
# JupyterLab answers. Backgrounded because docker run below holds the terminal.
open_when_ready() {
    for _ in $(seq 1 120); do
        sleep 1
        if port_busy "$PORT"; then
            if command -v open >/dev/null 2>&1; then open "$URL"
            elif command -v xdg-open >/dev/null 2>&1; then xdg-open "$URL" >/dev/null 2>&1 || true
            fi
            return
        fi
    done
}
open_when_ready &
OPENER_PID=$!
trap 'kill "$OPENER_PID" 2>/dev/null || true' EXIT

# --- start in the foreground -------------------------------------------------
# No -d: the terminal stays attached, shows the JupyterLab log, and Ctrl-C is
# forwarded to the container. --rm cleans up on exit.
#
# seccomp=unconfined : without it dbus fails and the XFCE desktop stays grey
# shm-size=2g        : Chromium crashes with Docker's default 64 MB
#
# Auth is disabled through a config file rather than a command-line flag.
# "--IdentityProvider.token=" with nothing after the equals sign is silently
# ignored and Jupyter falls back to a generated token. In a config file the
# empty string is written as str(), so no shell quoting is involved.
TTY_ARGS=()
[ -t 0 ] && TTY_ARGS+=(-it)

echo
say "starting $IMAGE"
echo "    $URL"
echo
echo "    In the Launcher click the Desktop icon: GenomeStudio starts by itself."
echo "    Files in this folder appear as Z:\\home\\jovyan\\shared"
echo
echo "    Press Ctrl-C here to stop the container."
echo

docker run --rm \
    ${TTY_ARGS[@]+"${TTY_ARGS[@]}"} \
    --name "$CONTAINER" \
    ${HOME_ARGS[@]+"${HOME_ARGS[@]}"} \
    --security-opt seccomp=unconfined \
    --shm-size=2g \
    --memory="$MEM" \
    -p "127.0.0.1:${PORT}:8888" \
    -v "$SHARED_DIR:/home/jovyan/shared" \
    -e HOME=/home/jovyan \
    -e GS_AUTOSTART=1 \
    "$IMAGE" \
    bash -c "printf 'c.IdentityProvider.token = str()\nc.ServerApp.password = str()\n' > /tmp/jupyter_noauth.py; exec /home/python4Jup/bin/jupyter lab --config=/tmp/jupyter_noauth.py --ip=0.0.0.0 --port=8888 --no-browser --allow-root --ServerApp.root_dir=/home/jovyan"

docker rm -f "$CONTAINER" >/dev/null 2>&1 || true
echo
say "stopped"
