# Tesla conversion alignment: isolated stages

For the later source audit, mesh/sizing oracle fixtures, and remaining porting
work, see [Python/C++ parity](PYTHON_PARITY.md).
The later correction for segmented rounded edges is documented in
[Continuous fitted boundaries](ARC_BOUNDARIES.md); the counts below describe
the earlier reconstruction before that correction.

The reference was regenerated on `root@157.157.221.29:36507` with
stl2curves 0.4.0 and OCP 8.0.1.1. The five Python modules used in conversion
match the local original source at commit
`9f71aaa16eb689e759435d2e1fbf2ba7468cadfb` after normalizing line endings.
The STL SHA-256 is
`001c8843d5f064a324d80597d255b755561d30598f23de319e87f4d3b6462368`.
Native diagnostics use OCCT 8.0.1, independently of the browser WASM build.

## Proven stage isolation

| Experiment | Faces | Lines | Circles | B-spline edges | BRep volume, mm³ |
| --- | ---: | ---: | ---: | ---: | ---: |
| Python final BRep, imported and exported natively | 247 | 446 | 57 | 96 | 31480.717742 |
| Python unsewn faces, joined, sewn and tidied natively | 247 | 446 | 57 | 96 | 31480.717742 |
| Python supports and triangle partition, existing native face builder | 247 | 383 | 41 | 549 | 31470.813664 |
| Python supports and partition, experimental Python-style boundary runs | 247 | 446 | 75 | 94 | 31480.749648 |

The existing native builder also creates one other curve in the third row.
The experimental fourth row is not production parity: its edge types still differ.
It does establish that boundary/trim construction, rather than STEP serialization
or the kernel version alone, accounts for much of the extra spline geometry.

Python's STEP round trip has volume 31479.973399 mm³; compare like formats
when measuring differences. The BRep/STEP volume discrepancy is already present
in the reference, rather than being introduced by the native port.

The final reference has 190 planes, 17 cylinders, four cones, four spheres,
and 32 B-spline surfaces. Its final successful attempt uses **57** features:
17 cylinders, four cones, four spheres, 26 freeforms and six blends. Its
`prepare()` output has only 53 features. Comparing solely against preparation
therefore misses the replacements and restored pieces in `finish()`.

## Confirmed implementation differences

- Python's default boundary graph uses shared 3D runs through original mesh
  points. Exact surface-intersection edges are an optional experiment, disabled
  by default. The native builder previously preferred surface intersections.
- Python's crowded-corner radius is fixed at 0.02 mm. Coupling this to sewing
  makes it grow to 0.4 mm on Tesla. The fixed threshold passes all four native
  suites and seven CAD fixtures, but does not change Tesla's final result alone.
- Python builds curved faces before planar neighbors; a native patch list's
  order does not guarantee that ordering. Shared-corner repair makes order
  observable, though changing order alone did not resolve the full conversion.
- Comparing the same 4652 triangles by coordinates, the production native
  partition assigns 207 triangles from Python freeforms to planes. Python
  keeps facets up to 2% of the largest facet available for smooth fitting.
- Python preserves replaced pieces and can turn failing freeform boundaries
  into N-sided blends during recovery. The native height-field-only recovery
  does not implement the same behavior.

## Reproduction artifacts

Ignored diagnostic files are in `.cache/alignment2/`: `oracle-export.py`,
`oracle.cpp`, `reference-build.inc`, `run-oracle.py`, `oracle.log`,
`reference-build.log`, `overlap.py`, and `overlap.log`. The `oracle/` directory
contains the accepted Python surfaces, mesh partition and unsewn faces.
`candidate*-*` files are experiments and must not be treated as accepted fixes.

Full local conversion alignment remains an acceptance requirement. Matching a
stage with Python's fitted supports does not prove that local fitting and
recovery produce those supports.

## Current source change and remaining blocker

The source now tries shared boundary runs for fitted surfaces, keeps small
planar facets available after analytic recognition, and tries the requested
height-field tolerance before stricter fits. The combined prototype passed all
four native suites and seven CAD fixtures. Its Tesla STEP is valid, with one
solid, 184 faces, 35 B-spline surfaces and 598 B-spline edges (previously 204,
52 and 753). Its STEP volume is 31474.116683 mm³. These are improvements, not
full parity with the reference's 247 faces, 32 B-spline surfaces and 96 B-spline
edges. Fresh source/API verification is recorded separately in
`.cache/alignment2/production-*.log`.

Boundary-only distance samples against the regenerated reference:

| Direction | Previous maximum, mm | Current maximum, mm | Current p95, mm |
| --- | ---: | ---: | ---: |
| Reference to native | 0.081668 | 0.073181 | 0.012642 |
| Native to reference | 0.117190 | 0.116264 | 0.008125 |

These are samples of edges and face interiors, not Hausdorff bounds. The
current candidate has 3199 reference probes and 3726 native probes.

The boundary-run path itself reaches a volume of about 31479.25 mm³, but still
leaves two seams open. Four free edges form two pairs with matching endpoints:
one side is a spline and the other a chord, about 2.8 and 3.3 mm long. Healing
produces a valid **shell with zero solids**, so it must not be accepted as a
closed conversion. The current converter correctly rejects it and uses the
older trim alternatives. N-sided blend recovery using the shared outline is
being investigated for these specific failed seams, following Python's
`blend_fallback` behavior. Accepting the shell or weakening the volume check
would hide the mismatch rather than fix it.
