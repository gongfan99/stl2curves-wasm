# Validation status

WASM has not been compiled or executed during implementation. The manual build
workflow has not been triggered. Browser loading, Emscripten 6.0.13, and the
WASM build of OpenCascade 8.0.1 still require artifact-based validation.

Current native acceptance checks use OpenCascade 8.0.1, built from the archive
whose SHA-512 matches the pinned vcpkg port. The headless native build disables
X11 and OpenGL. Other development dependencies are GCC 13.3, Eigen 3.4,
nanoflann 1.5.4, libzip 1.7.3, and pugixml 1.14. These checks do not establish
compatibility with the Emscripten build or the complete pinned vcpkg dependency
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
parity for every input; broader conversion parity and WASM runtime acceptance
remain validation targets.

Prepared Node and browser suites exercise both variants, the actual blob loader,
isolation-based variant selection, owned buffers, repeated conversions, STEP
re-import validation, 3MF units and archive failures, and actual background fitting
on disconnected cylinders.
The optional Node fixture suite applies all seven original regression limits.
Run these only after compiled artifacts are supplied and testing is requested,
using the commands in the README.
