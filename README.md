# stl2curves-wasm

C++17 mesh-to-STEP reconstruction with browser/Node.js WebAssembly bindings.
The conversion reference is [stl2curves](https://github.com/Duvdivan/stl2curves)
revision `9f71aaa16eb689e759435d2e1fbf2ba7468cadfb`.

The source includes STL/3MF reading, mesh preparation, analytic fitting, symmetry
proposals, swept and freeform surface candidates, shared-boundary reconstruction,
sewing, cavities, fusion, STEP serialization, round-trip validation, and faceted
recovery. Advanced recognition is heuristic. **Full upstream conversion parity is
an acceptance target, not an established result.** See [porting notes](docs/PORTING.md)
and the unchanged upstream regression limits in `tests/regression_cases.json`.
Current checks and unresolved acceptance gaps are recorded in
[validation status](docs/VALIDATION.md).

There are no Python runtime dependencies. The Studio UI, renderer, FreeCAD add-on,
Python packaging, and disk stage cache are outside this library.

## Build with GitHub Actions

The workflow is manual only: **Actions → Build WASM (manual) → Run workflow**.
It builds the variants in parallel on separate runners, then merges their
temporary artifacts into a single `stl2curves-wasm` artifact. Temporary artifacts
are deleted after the merge succeeds. Download and extract `stl2curves-wasm.zip`;
it contains the complete deployment layout at its root:

```text
loader.mjs
loader.d.mts
package.json
LICENSE
THIRD_PARTY_NOTICES.md
threaded/
  stl2curves.js
  stl2curves.wasm
  stl2curves.d.ts
single_thread/
  stl2curves.js
  stl2curves.wasm
  stl2curves.d.ts
licenses/
```

Any separate `stl2curves.worker.*` files emitted by the SDK are also packaged.
The workflow builds and packages assets; runtime tests are run separately against
downloaded artifacts. See [validation status](docs/VALIDATION.md) for tested builds
and the JavaScript initialization fix required by the first supplied artifacts.

Dependencies use vcpkg registry baseline `817892b0ecbc6af994dd6df999c73faf795adbf2`, including
OpenCascade **8.0.1**. The headless overlay retains that port's source hash and
patches, selects geometry/STEP/STL toolkits, and keeps C++ exception handling.
Emscripten is pinned to **6.0.13**. The workflow caches EMSDK and dependency
binaries and uses the runner's preinstalled vcpkg.

For a local build, activate that SDK and set `VCPKG_ROOT`:

```sh
sh scripts/build_wasm.sh                # both variants
sh scripts/build_wasm.sh threaded
sh scripts/build_wasm.sh single_thread
```

Builds use `build_wasm_threaded/` and `build_wasm_single_thread/`, with separate
vcpkg installations, and package into `dist/`. Generated outputs are ignored.

## Browser loading

`loader.mjs` fetches JavaScript and WASM, imports JavaScript through a blob URL,
supplies `locateFile` and `mainScriptUrlOrBlob`, and exports the initialized module
as `stl2curves`.
It selects `threaded/` when `SharedArrayBuffer` exists and
`crossOriginIsolated === true`; otherwise it selects `single_thread/`.
Initialization errors are reported instead of silently selecting another build.

Run conversion in an **application Worker**. Conversion is synchronous after
module initialization, including in the single-thread variant. A threaded module
prewarms `navigator.hardwareConcurrency` pthreads and schedules at most eight
fitting tasks concurrently, including the calling thread.

```js
// Application Worker; adjust the loader URL to your deployed asset location.
import { stl2curves } from './assets/loader.mjs';

self.onmessage = ({ data: { name, bytes } }) => {
  const result = stl2curves.convert(
    [{ name, data: new Uint8Array(bytes) }],
    { tolerance: 0.01, accuracy: 0.1, details: true }
  );
  self.postMessage(result, result.files.map(file => file.data.buffer));
};
```

Inputs and outputs are owned byte buffers. No Embind handles need deletion.
`convert` returns `{ files, errors }`, allowing successful objects and individual
failures in one batch. Each output includes diagnostics with solid validity,
volume, face/surface counts, repairs, remaining defects, detected features, and
fallback reasons. Detected features may be rejected during reconstruction;
volume, face counts, and `surfaceCounts` describe the exported STEP after re-import.

Serve threaded applications over HTTPS or localhost with:

```http
Cross-Origin-Opener-Policy: same-origin
Cross-Origin-Embedder-Policy: require-corp
```

Content Security Policy must allow the application's blob module/worker loading
and WebAssembly compilation. Keep the deployment subdirectories intact so any
SDK worker bootstrap assets remain resolvable. Call `stl2curves.dispose()` when
finished; it terminates the pthread pool and releases the loader's blob URLs.
For hard cancellation, terminate the outer application Worker and initialize a
new instance. The conversion time limit is cooperative, not a hard timeout for
an OpenCascade call.

## API and defaults

Both variants export the same factory and API. Types are in
`wasm/stl2curves.d.ts`. `buildInfo()` returns version, OCCT version, threading mode,
and the bounded worker count. Configuration errors throw; conversion failures
are returned in `errors` with input/object names, code, and message.

| Option | Default | Meaning |
| --- | --- | --- |
| `tolerance` | `0.01` | Sewing tolerance in millimeters |
| `accuracy` | `0.1` | Fitting accuracy in millimeters |
| `timeLimit` | `600` | Optional-refinement budget in seconds; `0` unlimited |
| `simplify` | `-1` | Automatic; `0` disabled; positive values in millimeters |
| `fuse` | `true` | Fuse reconstructed bodies and object parts |
| `curves` | `true` | Enable curved-surface recognition |
| `blends` | `true` | Enable optional swept/freeform/blend fitting |
| `repair` | `true` | Repair supported mesh defects |
| `trueSize` | `false` | Apply a common inferred design scale |
| `details` | `false` | Include per-feature detection details |
| `merge` | `""` | Also write a combined STEP using this name |

Automatic simplification uses `0.005 mm` above one million triangles.
Diagnostics `scale` is the inferred design factor; it is applied only when
`trueSize` is enabled. `parallelTasks` counts tasks executed by background workers.
Geometry failures retain diagnostics for work completed before the failure,
including repairs and unresolved defects. Input parsing errors may omit them.
3MF build objects become separate outputs with part labels retained through
repair; overlapping parts are reconstructed separately before fusion. Diagnostics
`meshVolume` is the sum of input part volumes before boolean fusion.

Node.js 24 or later can import a variant factory directly:

```js
import createStl2curves from './dist/single_thread/stl2curves.js';
const module = await createStl2curves();
try {
  const result = module.convert([{ name: 'part.stl', data: inputBytes }]);
} finally {
  module.dispose();
}
```

## Validation

Source checks do not load or execute WASM:

```sh
node scripts/check_sources.mjs
sh -n scripts/build_wasm.sh
# After installing the development dependencies:
npm run check:types
npm run test:post
```

Native tests require CMake, a C++17 compiler, and the manifest dependencies:

```sh
cmake -S . -B build_native -G Ninja \
  -DCMAKE_TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" \
  -DSTL2CURVES_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_native
ctest --test-dir build_native --output-on-failure
build_native/generate_fixtures build_native/fixtures
build_native/test_regression build_native/fixtures
```

The fixture generator also writes each original CAD shape as `.stl.brep` for
native comparisons. `inspect_detection INPUT.stl` reports recognized supports
and checks triangle ownership. `inspect_reference INPUT.stl REFERENCE.stl.brep`
reports face properties and sampled deviations from the original CAD geometry.

`STL2CURVES_SYSTEM_DEPS=ON` allows development checks with system OCCT packages;
it is rejected for WASM. Such checks do not establish compatibility with the
pinned OCCT 8.0.1/Emscripten toolchain.

To test compiled artifacts placed in `dist/`:

```sh
node wasm/test_node.mjs threaded
node wasm/test_node.mjs single_thread
# Include the upstream fixtures in each Node run:
STL2CURVES_FIXTURES=build_native/fixtures node wasm/test_node.mjs threaded
STL2CURVES_FIXTURES=build_native/fixtures node wasm/test_node.mjs single_thread
npm install
npx playwright install chromium
node wasm/test_browser.mjs
# Alternatively, use an installed Google Chrome:
STL2CURVES_BROWSER_CHANNEL=chrome node wasm/test_browser.mjs
```

`STL2CURVES_WASM_MODULE` can select another artifact's JavaScript path. Browser
tests exercise the real blob loader in an application Worker, both with and
without isolation headers. Node tests check ownership, repeated conversions,
batch failures, output validity, variant identity, and the optional regression
fixtures. STEP re-import validation runs inside the C++ conversion core.

License: MIT. See [third-party notices](THIRD_PARTY_NOTICES.md).
