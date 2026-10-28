# stl2curves-wasm

C++17 mesh-to-CAD reconstruction with browser/Node.js WebAssembly bindings.
Version 0.2.0 uses a model-first API: STL/3MF → retained OCCT model → explicit
STEP or GLB export. `convert()` no longer returns STEP files.
The conversion reference is [stl2curves](https://github.com/Duvdivan/stl2curves)
revision `9f71aaa16eb689e759435d2e1fbf2ba7468cadfb`.

The source includes STL/3MF reading, mesh preparation, analytic fitting, symmetry
proposals, swept and freeform surface candidates, shared-boundary reconstruction,
sewing, cavities, fusion, STEP serialization, round-trip validation, and faceted
recovery. Advanced recognition is heuristic. **Full upstream conversion parity is
an acceptance target, not an established result.** See [porting notes](docs/PORTING.md)
and the unchanged upstream regression limits in `tests/regression_cases.json`.

There are no Python runtime dependencies. The Studio UI, renderer, FreeCAD add-on,
Python packaging, and disk stage cache are outside this library.

## Build with GitHub Actions

The workflow is manual only: **Actions → Build WASM (manual) → Run workflow**.
Download and extract `stl2curves-wasm.zip`; it contains the complete deployment
layout at its root:

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
downloaded artifacts using the commands below.

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
import { stl2curves } from "./assets/loader.mjs";

const operations = new Set([
  'convert', 'exportStep', 'exportGltf',
  'getShapeProperties', 'getClosestPoints', 'releaseModel',
]);
self.onmessage = ({ data: { requestId, method, args = [] } }) => {
  try {
    if (!operations.has(method)) throw new Error('Unknown operation');
    const result = stl2curves[method](...args);
    const transfer = result.data instanceof Uint8Array ? [result.data.buffer] : [];
    self.postMessage({ requestId, result }, transfer);
  } catch (error) {
    self.postMessage({ requestId, error: error.message });
  }
};
```

Send a conversion request with `method: 'convert'` and
`args: [[{ name, data: new Uint8Array(bytes) }], options]`. Use the returned
model IDs in later export requests, and send picked selection metadata back to
this same Worker for queries. Model handles remain alive until release or disposal.

`convert` returns `{ models, errors }`; it does not serialize or re-import STEP.
Each model contains `modelId`, `name`, `root` (a selection reference), and
diagnostics for its retained OCCT geometry. Failed objects are reported alongside
successful models. A model can contain one solid or a compound of several solids.

Call `exportStep(modelId)` or `exportGltf(modelId)` independently and in either
order. Each call returns a named, owned byte buffer. STEP export performs bounded
recovery and fresh STEP re-import validation on an isolated copy; its diagnostics
describe the exported STEP. GLB meshing also uses isolated copies. Exports never
change model topology, diagnostics or BRepGraph identifiers.

Call `releaseModel(modelId)` when a model is no longer needed; it returns `false`
for an already released or unknown ID. Exported buffers remain valid after release.
No Embind objects require `.delete()`. `dispose()` releases all models and pthreads.
All model operations are synchronous, reject reentrant use, and fail after disposal.
An export failure leaves its model available for queries and other exports.

Serve threaded applications over HTTPS or localhost with:

```http
Cross-Origin-Opener-Policy: same-origin
Cross-Origin-Embedder-Policy: require-corp
```

Content Security Policy must allow the application's blob module/worker loading
and WebAssembly compilation. Keep the deployment subdirectories intact so any
SDK worker bootstrap assets remain resolvable. Call `stl2curves.dispose()` when
finished; it releases retained models, terminates the pthread pool, and releases blob URLs.
For hard cancellation, terminate the outer application Worker and initialize a
new instance. The conversion time limit is cooperative, not a hard timeout for
an OpenCascade call.

## API and defaults

Both variants export the same factory and API. Types are in
`wasm/stl2curves.d.ts`. `buildInfo()` returns version, OCCT version, threading mode,
and the bounded worker count. Configuration errors throw; conversion failures
are returned in `errors` with input/object names, code, and message.

| Option      | Default | Meaning                                                 |
| ----------- | ------- | ------------------------------------------------------- |
| `tolerance` | `0.01`  | Sewing tolerance in millimeters                         |
| `accuracy`  | `0.1`   | Fitting accuracy in millimeters                         |
| `timeLimit` | `600`   | Optional-refinement budget in seconds; `0` unlimited    |
| `simplify`  | `-1`    | Automatic; `0` disabled; positive values in millimeters |
| `fuse`      | `true`  | Fuse reconstructed bodies and object parts              |
| `curves`    | `true`  | Enable curved-surface recognition                       |
| `blends`    | `true`  | Enable optional swept/freeform/blend fitting            |
| `repair`    | `true`  | Repair supported mesh defects                           |
| `trueSize`  | `false` | Apply a common inferred design scale                    |
| `details`   | `false` | Include per-feature detection details                   |
| `merge`     | `""`    | Also retain a combined compound model with this name              |

Automatic simplification uses `0.005 mm` above one million triangles.
Diagnostics `scale` is the inferred design factor; it is applied only when
`trueSize` is enabled. `parallelTasks` counts tasks executed by background workers.
Geometry failures retain diagnostics for work completed before the failure,
including repairs and unresolved defects. Input parsing errors may omit them.
3MF build objects become separate models with part labels retained through
repair; overlapping parts are reconstructed separately before fusion. Diagnostics
`meshVolume` is the sum of input part volumes before boolean fusion.

Node.js 24 or later can import a variant factory directly:

```js
import createStl2curves from "./dist/single_thread/stl2curves.js";
const module = await createStl2curves();
try {
  const result = module.convert([{ name: "part.stl", data: inputBytes }]);
  for (const model of result.models) {
    const step = module.exportStep(model.modelId);
    const glb = module.exportGltf(model.modelId);
    const properties = module.getShapeProperties(model.root);
    module.releaseModel(model.modelId);
  }
} finally {
  module.dispose();
}
```


## Selectable GLB and exact queries

`exportGltf(modelId, options?)` produces binary glTF (`.glb`) without textures or
compression. The defaults are `linearDeflection: 0.1` mm,
`angularDeflection: 0.5` radians, and `faces`, `edges`, `vertices` all `true`.
Faces are meshes, edges are line segments, and vertices are point objects.
Compound/solid/shell groups are preserved. Degenerate geometry retains selectable
metadata-only nodes and produces warnings; ordinary face meshing failures throw.

Each selectable node's `extras` become Three.js `object.userData`:

```ts
interface GraphId { kind: number; index: number }
interface UsagePathStep { node: GraphId; ref?: GraphId; stepIndex: number }
interface ShapeRef { modelId: number; brepgraph_id?: GraphId; usagePath?: UsagePathStep[] }
```

Use the picked object's `userData` directly with `getShapeProperties(selection)`
or `getClosestPoints(selectionA, selectionB)`, sending it back to the same
application Worker when the WASM module lives there. Both functions resolve the
original OCCT geometry, independently of display tessellation. ID-only selections
work when placement is unambiguous; use the exported usage path for repeated
instances. Omitting `brepgraph_id` selects the whole model.

Graph node kinds include solid `0`, shell `1`, face `2`, wire `3`, edge `4`,
vertex `5`, compound `6`, and compsolid `7`. These are different from OCCT
`TopAbs_ShapeEnum` and `GeomAbs_CurveType`, used in property results.

`getClosestPoints(a, b)` supports selections within one model or across models
and returns `{ p1: [x,y,z], p2: [x,y,z], distance }`. Equal minima may return any
one minimizing pair. `getShapeProperties(ref)` returns
`{ properties: [{ name, value: number[] }] }`: `shape_type` for every shape,
`curve_type` for edges, `start_point`, `end_point`, `length` for bounded lines and
circles, and `center`, `radius`, `angle_degrees` for circles. Unsupported,
unbounded and degenerate edges provide applicable type metadata only.

Query coordinates and dimensions use model millimeters with Z up. The GLB root
maps `(x,y,z)` to `(x,z,-y)/1000` in meters with Y up, and records that convention
in `userData.coordinates`. Transform returned points using this root matrix to
display them. Viewer transforms applied later are not automatically included in
CAD queries.

IDs are scoped to the live model in its WASM instance and are never reused during
that instance's lifetime. Releasing a model invalidates its selections. A saved
GLB cannot recover CAD geometry after a reload; persistent reopening and graph
editing are outside this API. The model API requires OpenCascade 8.0.1 or later.

## Validation

Source checks do not load or execute WASM:

```sh
node scripts/check_sources.mjs
sh -n scripts/build_wasm.sh
# After installing the development dependencies:
npm run check:types
npm run test:post
npm run test:package
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

`STL2CURVES_SYSTEM_DEPS=ON` allows development checks with system OCCT 8.0.1+ packages;
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
fixtures. STEP re-import validation runs inside `exportStep()`. The suites also load GLB
with Three.js, raycast faces, and pass their metadata into exact geometry queries.

To check Three.js interoperability using a native build without loading WASM:

```sh
build_native/test_models build_native/cube.glb
node tests/test_native_gltf.mjs build_native/cube.glb build_native/picked.json
build_native/test_models build_native/cube.glb build_native/picked.json
```

License: MIT. See [third-party notices](THIRD_PARTY_NOTICES.md).
