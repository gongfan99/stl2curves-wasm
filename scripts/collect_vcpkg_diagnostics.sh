#!/usr/bin/env sh
set -eu

variant=${1:?Expected threaded or single_thread}
case "$variant" in threaded|single_thread) ;; *) echo 'Expected threaded or single_thread' >&2; exit 2 ;; esac
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
output="$root/.cache/vcpkg-diagnostics/$variant"
installed="$root/build_wasm_$variant/vcpkg_installed"
mkdir -p "$output"

# Missing tools and partial installs are expected when collecting a failed build.
version() {
  printf '\n$'
  printf ' %s' "$@"
  printf '\n'
  "$@" 2>&1 || true
}
{
  printf 'variant=%s\n' "$variant"
  printf 'commit=%s\n' "${GITHUB_SHA:-unknown}"
  printf 'run_id=%s\nrun_attempt=%s\n' "${GITHUB_RUN_ID:-unknown}" "${GITHUB_RUN_ATTEMPT:-unknown}"
  printf 'runner_os=%s\nrunner_arch=%s\n' "${RUNNER_OS:-unknown}" "${RUNNER_ARCH:-unknown}"
  printf 'image_os=%s\nimage_version=%s\n' "${ImageOS:-unknown}" "${ImageVersion:-unknown}"
  printf 'emscripten_version=%s\n' "${EM_VERSION:-unknown}"
  printf 'cache_exact_hit=%s\n' "${VCPKG_CACHE_HIT:-unknown}"
  if [ -n "${VCPKG_ROOT:-}" ]; then
    version "$VCPKG_ROOT/vcpkg" version
    version git -C "$VCPKG_ROOT" rev-parse HEAD
  fi
  version cmake --version
  version ninja --version
  version emcc --version
  version em++ --version
  version c++ --version
} > "$output/environment.txt"

copy_abi_files() {
  source_dir=$1
  target_dir=$2
  if [ -d "$source_dir" ]; then
    mkdir -p "$target_dir"
    # Preserve relative paths so host/target packages cannot overwrite each other.
    (cd "$source_dir" && find . -type f -name '*vcpkg_abi_info.txt' \
      -exec cp --parents -t "$target_dir" -- {} +)
  fi
}
copy_abi_files "$installed" "$output/installed"
if [ -f "$installed/vcpkg/status" ]; then
  cp "$installed/vcpkg/status" "$output/installed-status.txt"
fi
if [ -n "${VCPKG_ROOT:-}" ]; then
  # Buildtree ABI files may exist even when a package failed before installation.
  copy_abi_files "$VCPKG_ROOT/buildtrees" "$output/buildtrees"
  compiler_logs="$VCPKG_ROOT/buildtrees/detect_compiler"
  if [ -d "$compiler_logs" ]; then
    mkdir -p "$output/detect_compiler"
    (cd "$compiler_logs" && find . -type f -name '*.log' \
      -exec cp --parents -t "$output/detect_compiler" -- {} +)
  fi
fi
if [ -d "$root/.cache/vcpkg-binaries" ]; then
  find "$root/.cache/vcpkg-binaries" -type f -name '*.zip' -printf '%P\n' \
    | sort > "$output/available-archives.txt"
fi
printf 'Vcpkg diagnostics collected in %s\n' "$output"
