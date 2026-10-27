# Porting reference and acceptance

The reference is the sibling `stl2curves` checkout at
`9f71aaa16eb689e759435d2e1fbf2ba7468cadfb`. This repository has no Python runtime
dependency. Python is permitted only as a development oracle when comparing
against that checkout.

| Python subsystem | C++ implementation |
| --- | --- |
| STL reading, welding, repair, simplification, bodies | `cpp/mesh.cpp` |
| 3MF archive and slicer metadata | `cpp/read3mf.cpp` |
| Facets, analytic regions, radii and design size | `cpp/detect.cpp` |
| Pipes, threads, extrusion, freeform and blends | `cpp/advanced.cpp` |
| Shared boundaries, sewing, cavities, fusion, checks | `cpp/build.cpp` |
| Conversion orchestration and STEP | `cpp/converter.cpp`, `cpp/step.cpp` |
| Numeric optimization and worker execution | `cpp/numeric.cpp` |

The C++ algorithms use Eigen dense/sparse systems and robust finite-difference
least squares instead of NumPy/SciPy. Advanced surface candidates are measured
against mesh corners before use. Unrecognized regions retain planar geometry.
Recognition is heuristic: a detected feature is not proof that the final solid
retained it; inspect `diagnostics.surfaceCounts` and `diagnostics.fallbacks`.

Collection types use the current `NCollection` templates rather than relying on
deprecated headers being installed. `cpp/occt_collections.hpp` keeps older native
development builds compatible where handle-managed arrays still have distinct
RTTI classes. WASM dependencies remain pinned to OpenCascade 8.0.1.

3MF vertices and transform translations are converted to millimeters using the
owning model's unit; linear matrix entries remain dimensionless. This follows
the [3MF unit and transformation rules](https://github.com/3MFConsortium/spec_core/blob/master/3MF%20Core%20Specification.md#34-model),
including for component and build-item translations in nonmillimeter models.

Each 3MF part stores only the vertices referenced by its triangles. A shared
measurement mesh uses explicit vertex and triangle offsets for assembly sizing,
then restores the scaled coordinates to each part. Part storage therefore scales
with referenced geometry rather than copying the full assembly into every part.
Measurement-only arrays are released before BRep reconstruction begins.

Native fixtures reproduce the seven upstream CAD designs and their original
acceptance limits in `tests/regression_cases.json`. Do not weaken these limits
to obtain a passing suite. Full parity requires validity, volume, face-count,
and feature-specific regressions to pass, rather than merely producing STEP.

Time limits are cooperative budgets for optional fitting and recovery. They do
not interrupt an OpenCascade call already executing. Python process termination
has no direct equivalent for an individual C++ pthread. Applications requiring
hard cancellation must terminate the outer application Worker and create a new
module instance; pending results from that Worker are then abandoned.

The new implementation does not reproduce every upstream repair and recovery
heuristic verbatim. Self-intersecting/nonmanifold meshes can still fail instead
of yielding a solid. Full conversion parity must remain an explicit acceptance
gate; build success and basic tests alone do not establish it.

WASM compilation is performed through the manual GitHub workflow. Supplied
artifacts have been tested in Node and browsers; see `VALIDATION.md` for scope,
results and the JavaScript initialization correction required by the first build.
