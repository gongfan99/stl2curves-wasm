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
The source audit, executable Python-oracle checks, and remaining algorithm gaps
are recorded in [Python/C++ parity](PYTHON_PARITY.md).

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
Between-corner bulge checks also follow `Freeform.signed()`: divide the height
residual by the implicit gradient norm. Their accuracy-scaled allowance stays
independent of stricter corner-fit trials; coupling the two over-refines steep
smooth patches into splines that follow the tessellation.
Several precise height-field pieces can also be replaced by one fit at the
requested accuracy when it uses fewer control points overall. This removes
artificial partition boundaries while retaining precise single patches and
protecting exact or swept child surfaces from replacement. A sampled comparison
against the previous fitted interiors also bounds the change by the same
bulge allowance; passing the mesh-corner check alone does not establish this.

As in Python `attempt()` and `sewing_cap()`, fitted-patch sewing adapts to measured
mesh-to-support gaps, capped by 0.2 mm and twice the median mesh edge, with the
configured tolerance as its lower bound. UV trims, repaired 3D outlines, and bounded
kernel splitting are tried before retaining triangulated trims on the same support.
The native kernel-cut fallback is restricted to small supports because an in-process
OpenCascade call cannot use Python's worker timeout.

Fitted outlines first use shared 3D runs grouped by the neighboring face and split
at 30-degree turns, following `_runs`, `_split_sharp`, and `_even`. The run pool
uses natural primitive lines/arcs where available and interpolates other whole
runs through the original mesh points, densifying long steps. Curved faces are
built before their planar neighbors. Touching outlines are split into unpinched
pieces on the same support. Outline, kernel-split, polygon and wider-split
construction follow Python's order, with the existing limit on expensive cuts.
The result must pass closed-solid and corrected-volume checks; otherwise the
established analytic/projected construction remains available.

In that fallback path, spherical outlines also use shared 3D runs. A single
sphere-side edge may bow onto its support while the neighboring freeform keeps
its mesh chord, as in `_run_edge`; sewing can retain the chord within the measured
gap. Chord preference is restricted to isolated sphere/freeform seams. Applying
it to all freeform boundaries changes unrelated regression geometry. Outline
repair and cleanup retain Python's 1 mm maximum face-tolerance gate alongside
the existing validity and corrected-volume checks.

Shared-corner crowding follows Python's fixed `SHARED_CROWD = 0.02` mm,
independently of adaptive sewing. Using twice the sewing tolerance instead
silently expands the exclusion radius to 0.4 mm on loosely sewn fitted parts.

Final cleanup follows `_tidy` and `_bad_merges`: merge at 0.001 mm and 0.1 degrees,
keep seams of invalid merged faces apart, and check volume integration within 0.1%
before the unchanged corrected-volume and closed-solid checks. Cleanup copies
topology while preserving shared surface identity; copying geometry independently
for every face prevents the kernel from merging B-spline faces on one support.

GLB export retries faces omitted by the default body tessellation with a finer,
area-scaled deflection. This handles tiny valid trims left by reconstruction
without removing faces or changing the retained model's graph identities.

Analytic cone seeds use Python `_cone`/`_line_on_direction` with a fixed,
normal-derived axis. Whole-island proposals retain nonlinear fitting, as in
`_whole_region`/`_tube_fit`. Optimizing tiny local cone seeds can incorrectly
claim spherical bands with only two vertex rings. A subsequent fixed-axis pass
follows `_axis_patches`, using established cylinder/cone axes to recover smaller
tapered regions before freeform fitting. The late sphere-to-cone rewrite has
been removed because it undoes this recognition sequence.

After analytic recognition, small planar facets up to 2% of the largest facet
remain available to smooth fitting, as in `blends.add_blends()`. Applying this
before analytic recognition can lose established cylinders and cones. Height
fields first try the requested fitting tolerance, as in `freeform.fit()`;
stricter fits remain fallbacks. Trying strict fits first creates extra pieces
and artificial boundaries on otherwise smooth regions.

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
