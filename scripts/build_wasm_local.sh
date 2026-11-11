#!/usr/bin/env bash
set -euo pipefail
variant=${1:-single_thread}
workspace=${2:?Pass the Linux toolchain/build workspace}
jobs=${3:-2}
case "$variant" in all|single_thread|threaded) ;; *) echo 'Invalid build variant' >&2; exit 2 ;; esac
case "$workspace" in /*) ;; *) echo 'Workspace must be an absolute Linux path' >&2; exit 2 ;; esac
case "$workspace" in /mnt/c|/mnt/c/*|/) echo 'Choose a build workspace outside C:' >&2; exit 2 ;; esac
if [[ ! "$jobs" =~ ^[1-9][0-9]*$ ]]; then echo 'Build jobs must be a positive integer' >&2; exit 2; fi
workspace=$(realpath -- "$workspace")
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
mkdir -p "$workspace"/{tmp,cache,downloads,binary-cache,build,dist,logs}
export TMPDIR="$workspace/tmp"
export XDG_CACHE_HOME="$workspace/cache"
export VCPKG_DOWNLOADS="$workspace/downloads"
export VCPKG_DEFAULT_BINARY_CACHE="$workspace/binary-cache"
export VCPKG_BINARY_SOURCES="clear;files,$workspace/binary-cache,readwrite"
export VCPKG_ROOT="$workspace/vcpkg"
export EM_CACHE="$workspace/cache/emscripten"
export WASM_BUILD_ROOT="$workspace/build"
export WASM_DIST_DIR="$workspace/dist"
export BUILD_JOBS="$jobs"
export VCPKG_DISABLE_METRICS=1
export PACKAGE_VARIANT_ONLY=ON
if [[ "$variant" == all ]]; then
    export WASM_DIST_DIR="$workspace/variant-parts"
else
    export WASM_DIST_DIR="$workspace/dist"
fi
export WASM_BUILD_ROOT="$workspace/build"
if [[ ! -f "$workspace/emsdk/emsdk_env.sh" || ! -x "$VCPKG_ROOT/vcpkg" ]]; then
    echo "Toolchains are not installed in $workspace" >&2
    exit 1
fi
# The SDK's activation script also supports shells without nounset.
set +u
source "$workspace/emsdk/emsdk_env.sh" >/dev/null
set -u
# EMSDK activation can unset these during its config sanity check; restore all
# large caches/temp paths afterward so downloads and intermediates stay on E:.
export TMPDIR="$workspace/tmp"
export XDG_CACHE_HOME="$workspace/cache"
export VCPKG_DOWNLOADS="$workspace/downloads"
export VCPKG_DEFAULT_BINARY_CACHE="$workspace/binary-cache"
export VCPKG_BINARY_SOURCES="clear;files,$workspace/binary-cache,readwrite"
export VCPKG_ROOT="$workspace/vcpkg"
export EM_CACHE="$workspace/cache/emscripten"
export PACKAGE_VARIANT_ONLY=ON
export CMAKE_MAKE_PROGRAM=/usr/bin/ninja
export GIT_CONFIG_GLOBAL="$workspace/gitconfig"
git config --file "$GIT_CONFIG_GLOBAL" \
    url."file://$workspace/vcpkg".insteadOf https://github.com/microsoft/vcpkg
if [[ -x /snap/cmake/current/bin/cmake ]]; then
    export PATH="/snap/cmake/current/bin:$PATH"
fi
command -v cmake >/dev/null
command -v ninja >/dev/null
log="$workspace/logs/build-$variant-$(date +%Y%m%d-%H%M%S).log"
printf 'Source: %s\nBuild: %s\nArtifacts: %s\nLog: %s\n' "$root" "$WASM_BUILD_ROOT" "$WASM_DIST_DIR" "$log"
bash "$root/scripts/build_wasm.sh" "$variant" 2>&1 | tee "$log"
if [[ "$variant" == all ]]; then
    stamp=$(date +%Y%m%d-%H%M%S)
    assembled="$workspace/build/assembled-$stamp"
    mkdir -p "$assembled"
    cmake "-DROOT=$root" "-DINPUT_DIR=$WASM_DIST_DIR" \
        "-DOUTPUT_DIR=$assembled" -P "$root/scripts/assemble_wasm.cmake"
    if [[ -e "$workspace/dist" ]]; then
        mv -- "$workspace/dist" "$workspace/build/previous-dist-$stamp"
    fi
    mv -- "$assembled" "$workspace/dist"
    printf 'Combined package: %s\n' "$workspace/dist"
fi
