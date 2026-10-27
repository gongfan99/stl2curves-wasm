# Validation status

On 2026-10-26, user-supplied Emscripten 6.0.13/OpenCascade 8.0.1 artifacts were
tested on Windows with Node.js 24.15.0 and headless browsers. No local WASM
compilation or workflow dispatch was performed during this validation.

The original JavaScript failed initialization in both variants because our
`--post-js` wrapper called `addOnPostRun`, which the generated SDK output did
not define. In this SDK, post-js runs after awaited runtime initialization.
`wasm/post.js` now wraps initialized Embind exports directly and skips wrapping
in pthread bootstrap instances. Tests used local copies of the supplied assets
with only this post-js section replaced; the downloaded originals were preserved.
The source fix still needs to be included in the next workflow build.

| Runtime check | single_thread | threaded |
| --- | --- | --- |
| Node API, ownership, errors, repeated conversions and disposal | pass | pass |
| Node 3MF units/names, malformed archives and mesh defects | pass | pass |
| All seven original regression volume and face limits, one module per variant | pass | pass |
| Chrome blob loader and application Worker | pass | pass |
| Playwright Chromium 141.0.7390.37 blob loader and application Worker | pass | pass |
| Background fitting | no background tasks, as expected | pthread tasks observed |

Browser runs exercised ordinary and cross-origin-isolated pages, automatic
variant selection, owned STEP buffers, STEP re-import validation, 3MF conversion,
repeated conversions, structured failures and twelve disconnected analytic bodies.
Tests also assert that conversion is rejected after disposal in Node. Separate
JavaScript lifecycle tests cover an absent SDK callback helper, error recovery,
busy guards, idempotent disposal and pthread bootstrap without Embind exports.

The unchanged WASM binaries have these SHA-256 digests:

| Variant | SHA-256 |
| --- | --- |
| single_thread | `01c5bc7de1381b443b075fd8ebe2b52ab2d8675734ebd0da894b82f979797670` |
| threaded | `4772e99662b2737aa618920f6f20d8e07fb9fd26bb69f030f6f6175119052b42` |

Current native acceptance checks use OpenCascade 8.0.1, built from the archive
whose SHA-512 matches the pinned vcpkg port. The headless native build disables
X11 and OpenGL. Other development dependencies are GCC 13.3, Eigen 3.4,
nanoflann 1.5.4, libzip 1.7.3, and pugixml 1.14. These native checks alone do not
establish compatibility with Emscripten or the complete pinned vcpkg dependency
set. Earlier development checks also used Ubuntu's OpenCascade 7.6.3 packages.

Both native scheduling configurations pass the core and advanced fitting tests.
Coverage includes in-memory STL, mesh repair, symmetry, nonlinear fitting,
worker exception propagation, arbitrary helical pitch, rotated left-handed threads,
closed noncircular extrusion profiles, swept pipes,
freeform fitting, 3MF references/transforms/units/helper filtering, cavities,
closed single-face spheres, overlapping parts, conversion options, repeated
conversions, partial batch failures, design sizing, multipart sizing and component
positions after STEP re-import, and STEP round trips.

JavaScript/JSON syntax, TypeScript declarations and examples, shell syntax, and
the manual-only workflow structure have been checked without loading WASM.

The seven generated designs retain the original volume-error and maximum-face
limits in `tests/regression_cases.json`. All seven pass with OpenCascade 8.0.1
in both native scheduling configurations, with identical reported results:

| Case | Exported STEP volume error (mm³) | Faces / maximum | Result |
| --- | ---: | ---: | --- |
| fine | -0.000016 | 29 / 29 | pass |
| coarse | -0.000016 | 29 / 29 | pass |
| hard | -0.039908 | 31 / 31 | pass |
| lid | +0.000924 | 45 / 45 | pass |
| sphere | -0.000033 | 34 / 34 | pass |
| corners | -0.048705 | 26 / 30 | pass |
| corners2 | -0.019346 | 26 / 34 | pass |

All seven outputs retain reconstructed surfaces, report no faceted fallback,
and pass fresh STEP re-import validity checks. Blend junctions retain exact
analytic region boundaries; very slender transition facets are reassigned only
when a neighboring surface explains their vertices and preserves orientation.
The limits have not been relaxed. These fixtures do not establish full upstream
parity for every input; broader conversion parity remains a validation target.

Node and browser suites exercise both variants, the actual blob loader,
isolation-based variant selection, owned buffers, repeated conversions, STEP
re-import validation, 3MF units and archive failures, and actual background fitting
on disconnected cylinders.
The optional Node fixture suite applies all seven original regression limits.
Run these against compiled artifacts using the commands in the README.
