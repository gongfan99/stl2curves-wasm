# Python/C++ parity audit

Full conversion parity is **not established**. This document distinguishes
verified corrections from algorithms that still require a substantive port.
The Python reference is the sibling `stl2curves` checkout at
`9f71aaa16eb689e759435d2e1fbf2ba7468cadfb` (version 0.4.0).
The installed reference supplied for testing uses OCP 8.0.1.1. All its Python
module hashes match the local checkout after normalizing CRLF to LF.
The hashes are preserved in `tests/fixtures/python_mesh_manifest.json`.

The single-thread `test_corners2` facet port and geometry-level oracle checks
are documented in [CORNERS2_ALIGNMENT.md](CORNERS2_ALIGNMENT.md). Its 20 curved
ownership groups, fitted surfaces, classifications and 137 boundary edges now
match the checked Python reference. This does not establish full-library parity.

## Corrections covered by executable checks

| Area | Correction |
| --- | --- |
| Duplicate triangles | Opposite windings cancel in pairs; only one unmatched copy survives, with its original winding and position in the input order. |
| Slivers | Split the triangle across a sliver's opposite directed edge instead of deleting the sliver and opening a T-junction. |
| Cracks | Join open-boundary corners within `max(0.001, size * 1e-5)`, retaining Python's representative-selection order. |
| Flaps and knots | Remove smaller coplanar flaps at overused edges; collapse tiny nonmanifold knots only when topology and normal checks improve and the remaining overused-edge count passes the 25% gate. |
| Orientation | Prefer the component's majority winding; turn an inward closed component outward only if it is not contained by another closed component. Preserve cavity facing. |
| Holes | Match directed-edge traversal order, minimum-area dynamic programming through 60 corners, and best-fit-plane ear clipping for larger loops. Remove the old planar-only and 128-corner restrictions. |
| Bodies | Group overlapping bounding boxes so cavity shells stay with their hosts; pair faces by turns into material at overused edges and color touching bodies apart before sewing. Actually invoke this stage in conversion and preserve the original order of corner samples while compacting each group. |
| Automatic simplification | Choose the million-triangle threshold for the repaired input object, then simplify each 3MF part separately. |
| Simplification algorithm | Replace one-pass short-edge deletion with up to 200 quadric-ranked half-edge collapse passes, independent neighborhoods, a two-corner link condition, 20-degree turn and minimum-shape checks, and reattachment of every previously removed corner within the requested error. Original coordinates remain unchanged. |
| Analysis order | Recognize analytic features, infer the common design size, reanalyze scaled meshes, snap supports, then run optional fitting. Previously optional surfaces voted on size and scaling transformed already fitted surfaces. |
| Size measurements | Use the biggest mesh facets, including unclaimed facets, for parallel-face gaps; measure cavity walls within their reconstruction group. Skip design sizing with curves disabled. |
| Snapping | Use area-weighted groups within 0.2%, fixed-radius position refits, and optional millimeter/inch rounding. Snap near-axis supports at the reference's 0.1-degree threshold. |
| Symmetry | Remove the additional default mirrored recognition pass: the reference has `features.MIRROR = False`. Symmetry discovery remains available. |
| Fusion | Use the reference's fixed `1e-5` boolean fuzz, batch operands, validate material volume, and preserve checked bodies when a union fails. Distinguish within-group overlap, 3MF part overlap, and the volume-preserving final join of distinct STL groups. |
| Volume integration | Use Python's adaptive relative precision `VOLUME_EPS = 1e-5`; check volume estimates at that precision. |

The normal native CTest suite includes 69 committed Python mesh/body/simplification fixtures
and five sizing/snapping fixtures. These need no Python runtime to execute.
Large planar holes are compared by oriented boundaries on each plane: different
SVD sign choices can select different internal diagonals on the same exact
planar polygon. All other mesh fixtures compare the oriented triangles and
coordinates directly. Body tests compare the full triangle partition, not just
the number of groups.
The simplification cases match Python's oriented triangles exactly: planar
cubes with 192 and 588 triangles both reduce to 12, and a smooth graph reduces
from 128 to 34. An independent global distance check also verifies every
original corner against the final simplified surface at the requested 0.01 mm.

Additional ignored fixtures in `.cache/python-parity/mesh_cases.txt` compare
repair, splitting and knot handling on the server's Tesla, joystick, and
`conroller.stl` meshes. All 72 cases in that file pass. In particular, the
controller's repaired mesh has 21,244 triangles in both implementations, from
21,270 input triangles.

The final OCCT 8 builds pass all six CTest suites and all seven CAD regressions
in both threaded and single-thread modes. The original volume and face-count
limits are unchanged. GLB export also passes on those seven fixtures; the
spherical fixture reports the existing degenerate-geometry warnings. Final
build, CTest and CAD logs are in `.cache/python-parity/final-*.log`.

## Full-conversion comparisons

The joystick mesh has SHA-256
`57f013a7d74b8b2e70db21b3b6b94eb64a5c30ef6f10ff66676e89982f0cba1d`.
Both implementations produce one valid solid with 31 faces (eight cylinders
and 23 planes), 16 circular edges, and 56 lines. Native BRep volume is
3100.759374570 mm³; Python BRep volume is 3100.759292894 mm³.

The earlier Tesla comparison has representation differences. The fresh Python result has
247 faces, 17 cylinders, four cones, four spheres, 32 B-spline surfaces, and
96 B-spline edges. The native result still has 184 faces, 13 cylinders, four
cones, four spheres, and 35 B-spline surfaces. See `ALIGNMENT.md` for the
earlier stage-isolation experiments. Matching the repaired mesh does not
establish matching detection, trims, or recovery.
The earlier native BRep has 598 B-spline edges and volume 31475.335874110 mm³,
versus the Python BRep's 96 B-spline edges and 31480.716101567 mm³. The native
STEP volume is 31474.129931718 mm³. Counts alone do not establish a shape defect;
fewer faces can be preferable. The segmented rounded boundaries were a concrete
defect. Their correction and geometric regression checks are described in
[Continuous fitted boundaries](ARC_BOUNDARIES.md).

The Python oracle itself fails to build a closed solid from the server's
`conroller.stl`; its successful mesh-stage comparisons must not be reported
as a successful Python full conversion.

## Remaining algorithm differences

| Reference subsystem | Remaining work |
| --- | --- |
| `repair.py` | Port crossing-triangle detection, noisy sliver handling restricted to crossing zones, and bounded self-crossing-speck excision with patch-shift checks. The current repair helpers do not implement those stages. |
| `simplify.py` | The quadric collapse algorithm and removed-corner error tracking are ported and covered by the oracle cases above. Broader large-mesh coverage remains necessary; floating-point ties in collapse ranking are not guaranteed bit-for-bit across Eigen and NumPy. |
| `features.py`, `fillets.py` | Align the complete facet/noise grouping, recognition/growth order, loose/noisy passes, anchored fillet handling, ownership and acceptance rules. Current geometric heuristics and late support refinements are not a line-for-line port. |
| `threads.py`, `pipes.py`, `extrude.py` | Recognition and surface construction still use native heuristics. Passing the synthetic thread/pipe/extrusion tests does not prove parity with the reference's segmentation and recovery. |
| `freeform.py` | The plane-based weighted P-spline fit is implemented; cylindrical height-field bases and all corresponding parameterization/construction paths are still missing. |
| `blends.py` | Port facet-based candidate selection, whole replaced-piece ownership, small N-sided blends, and restoration of original child features. Height-field splitting alone is not this subsystem. |
| `build.py` | Complete the default shared-run path for every feature class, N-sided `BRepFill_Filling` construction and its mesh/edge-gap checks, and the associated fallback/split paths. The established analytic/intersection alternatives still differ from Python's default. |
| `convert.py` | Recovery still searches subsets of fitted supports rather than reproducing Python's carving, piece restoration, half-rings, blend/pipe fallback and blame-guided second chances. Native in-process OCCT calls also cannot reproduce hard process timeouts. |
| `read3mf.py` | The native reader intentionally retains stricter validation and standards-correct unit scaling of translations. Those are documented API/format differences, not evidence of identical behavior on every malformed or nonmillimeter input. |

The retained-model API, explicit STEP/GLB exports, browser UI and absence of a
Python runtime dependency are intentional library differences. This audit does
not change those interfaces. No WASM artifact has been rebuilt by these native
checks.

## Reproduce

On a development machine with the reference Python package installed:

```sh
python tests/python_mesh_oracle.py tests/fixtures --source ../stl2curves
python tests/python_sizing_oracle.py tests/fixtures/python_sizing_cases.txt --source ../stl2curves
```

For additional real meshes, generate into an ignored directory:

```sh
python tests/python_mesh_oracle.py .cache/python-parity --source ../stl2curves --stl part.stl
python tests/python_conversion_oracle.py .cache/python-reference part.stl --source ../stl2curves
```

Build the native test targets with the dependencies described in the README:

```sh
ctest --test-dir build_native_occt8 --output-on-failure
build_native_occt8/test_regression build_native_occt8/fixtures
build_native_occt8/test_python_mesh_parity .cache/python-parity/python_mesh_cases.txt
build_native_occt8/inspect_conversion part.stl .cache/native-part.step
```

`inspect_conversion` saves both native BRep and STEP and reports each format's
volume separately, plus the retained surface/edge counts. Compare BRep to BRep
or STEP to STEP; the reference itself can change volume during STEP round trips.
