#!/usr/bin/env sh
set -eu
: "${VCPKG_ROOT:?Set VCPKG_ROOT}"
: "${EMSDK:?Activate the pinned Emscripten SDK first}"
sdk_banner=$(emcc --version | head -n 1)
case " $sdk_banner " in
  *' 6.0.13 '*) ;;
  *) echo 'Expected Emscripten 6.0.13; activate the pinned SDK before building' >&2; exit 1 ;;
esac
export VCPKG_MAX_CONCURRENCY="${BUILD_JOBS:-2}"
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
variant=${1:-all}
case "$variant" in all|threaded|single_thread) ;; *) echo 'Expected all, threaded, or single_thread' >&2; exit 2 ;; esac
build_variant() {
  flavor=$1
  threads=$2
  cmake -S "$root" -B "$root/build_wasm_$flavor" -G Ninja \
    "-DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" \
    "-DVCPKG_CHAINLOAD_TOOLCHAIN_FILE=$root/vcpkg-triplets/emscripten-wrapper.cmake" \
    "-DVCPKG_TARGET_TRIPLET=wasm32-emscripten-$flavor" \
    "-DVCPKG_INSTALLED_DIR=$root/build_wasm_$flavor/vcpkg_installed" \
    -DVCPKG_INSTALL_OPTIONS=--allow-unsupported \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \
    "-DWASM_THREADS=$threads"
  cmake --build "$root/build_wasm_$flavor" --parallel "${BUILD_JOBS:-2}"
  cmake -DROOT="$root" -DVARIANT="$flavor" -P "$root/scripts/package_wasm.cmake"
}
if [ "$variant" = all ] || [ "$variant" = threaded ]; then build_variant threaded ON; fi
if [ "$variant" = all ] || [ "$variant" = single_thread ]; then build_variant single_thread OFF; fi
