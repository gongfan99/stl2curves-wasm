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

Time limits are cooperative budgets for optional fitting and recovery. The solid
reconstruction budget starts after analysis, so expensive fitting cannot consume
the time available to recover valid curved patches. Recovery preserves shared
boundaries and uses invalid faces and seams to identify which supports to rebuild
as facets, including defective fitted planes next to otherwise valid curves.
The boundary builder follows Python `build.py`'s `_plane_face`, `_rejoin`, and
`_without_collapsed`: probe repairs before modifying shared geometry, join exact
boundary copies before sewing, and remove short edges collapsed by sewing while
preserving genuine closed circles. Recovery follows `convert.py`'s `culprits`
ordering by reverting suspect blends before their analytic neighbors, then
retrying patches that may have been blamed for a neighbor's failure.
These budgets do not interrupt an OpenCascade call already executing. Python process termination
has no direct equivalent for an individual C++ pthread. Applications requiring
hard cancellation must terminate the outer application Worker and create a new
module instance; pending results from that Worker are then abandoned.

Freeform fitting follows `freeform.py`'s weighted corner/interior samples and
accuracy scaling (tuned at 0.02 mm, capped at 0.1 mm). Regions split at 20-degree
creases and, when necessary, by the area-weighted principal variation of their
normals into connected pieces. Precise reconstructed pieces are retained; the
looser accuracy-scaled fit recovers areas those pieces do not explain. Height-field
trim charts use their affine base coordinates, as `Freeform.uvh()` does, rather
than moving each chart vertex independently to a nearest surface point.

As in Python `attempt()` and `sewing_cap()`, fitted-patch sewing adapts to measured
mesh-to-support gaps, capped by 0.2 mm and twice the median mesh edge, with the
configured tolerance as its lower bound. UV trims, repaired 3D outlines, and bounded
kernel splitting are tried before retaining triangulated trims on the same support.
The native kernel-cut fallback is restricted to small supports because an in-process
OpenCascade call cannot use Python's worker timeout.

Final cleanup follows `_tidy` and `_bad_merges`: merge at 0.001 mm and 0.1 degrees,
keep seams of invalid merged faces apart, and check volume integration within 0.1%
before the unchanged corrected-volume and closed-solid checks. Cleanup copies
topology while preserving shared surface identity; copying geometry independently
for every face prevents the kernel from merging B-spline faces on one support.

GLB export retries faces omitted by the default body tessellation with a finer,
area-scaled deflection. This handles tiny valid trims left by reconstruction
without removing faces or changing the retained model's graph identities.

The new implementation does not reproduce every upstream repair and recovery
heuristic verbatim. Self-intersecting/nonmanifold meshes can still fail instead
of yielding a solid. Full conversion parity must remain an explicit acceptance
gate; build success and basic tests alone do not establish it.

WASM compilation is performed through the manual GitHub workflow. Runtime tests
run separately against the compiled artifacts using the commands in the README.

The 0.2.0 public API retains reconstructed OCCT models before serialization.
`convert()` performs mesh preparation, reconstruction, fusion and native validity
checks, then builds a BRepGraph. STEP validation/recovery runs only in
`exportStep()` on a deep copy. GLB meshing and exact selection queries share the
retained model's immutable graph identity; model release is explicit in WASM and
RAII-owned in native C++. New graph functionality requires OCCT 8.0.1.
