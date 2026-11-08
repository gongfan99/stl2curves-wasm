# test_corners2: facet, surface and boundary alignment

Checks use the single-thread native OCCT 8 build only. The reference is the
original sibling Python checkout at `9f71aaa16eb689e759435d2e1fbf2ba7468cadfb`
(version 0.4.0). The replacement server's installed module hashes match it.
Both conversions use tolerance 0.01 mm, accuracy 0.1 mm and fusion disabled.
The input SHA-256 is
`5920038bc385cdecd140851a29348ed0aa8c44674dc3c9d8b231ad02d83f4b3d`.

## Cause and corrections

The initial 34-versus-30 face discrepancy came from four unclaimed needles
inside the corner patches. Porting `blends._islands` fixed that count but left
96 triangles assigned to different cylinders/freeform patches.

The former native detector fitted triangles directly, rejected small analytic
supports, rounded radii early and grouped transitions with different rules.
The new facet path follows Python's rolling-ball fillets, small-feature growth,
same-surface merging, pinholes, snapping, fixed-radius refits and straddler
ownership. Transition searches replace whole analytic pieces, including their
needles. It reproduces all 16 initial fillet groups and all 20 final curved
patch ownership groups.

| Native implementation | Original logic |
| --- | --- |
| `feature_facets.hpp`, `feature_order.hpp` | Facet geometry and CPython integer-set neighbor traversal |
| `feature_geometry.hpp` | Models, pair hypotheses, growth, refitting, trimming and evidence |
| `feature_fillets.hpp` | Rolling-ball fillets between flat faces and perpendicular walls |
| `feature_merge.hpp`, `feature_cleanup.hpp` | Merging, pinholes, snapping and straddlers |
| `feature_transitions.hpp` | Whole-piece transition selection |
| `feature_corner.hpp`, `reference_least_squares.hpp` | Corner-cylinder trials and unbounded linear-loss TRF rules |
| `feature_reference.inc` | Integration and natural/generous revolved face construction |

Height-field fitting now uses facet-weighted normals, the projected-point SVD
basis, samples grouped by subdivision count and the reference knot domain.
Small intermediate analytic models remain lazy until their transitions replace
them; some valid revolved profiles cannot be constructed as direct OCCT tori.

Reconstruction corrections preserve the trimmed/natural classification and
angular/profile bounds, split the reference's generous revolved surface,
use rim fitting error for splitter fuzz, and retain Python's run seams.
Splitting a full analytic cylinder had displaced one line by 0.000295 mm;
legacy C1 concatenation also reduced 31 spline edges to 29. Adaptive volume
integration now uses Python's `convert.VOLUME_EPS = 1e-5`, rather than 1e-7.

## Executable evidence

`tests/python_corner_oracle.py` generates the full repaired mesh, ownership,
classification, analytic parameters, spline surface/normal samples, boundary
samples and volume from Python. The saved fixture and JSON manifest record
source hashes and options. Native tests need no Python runtime.

The oracle checks:

- Exact triangle ownership for all 20 curved patches, including natural/trimmed classification.
- Analytic frames and radii within 1e-9.
- Four spline surfaces sampled at 81 points each: measured maximum gap about **8.3e-12 mm**.
- All **137 edges**: **55 lines, 51 circles, 31 splines**. Seventeen samples per reference edge give a maximum gap about **3.4e-12 mm** against native reconstruction.
- One valid solid with **30 faces**: 10 planes, 12 cylinders, four tori and four B-spline surfaces, without faceted recovery.
- Volume at the same adaptive integration precision.

The adaptive volume estimates differ by about 0.00047 mm^3 (relative 3.3e-8),
within the reference's 1e-5 integration precision. Matching integration
settings does not require identical quadrature estimates.

Geometry assertions allow 1e-7 mm, independently of the conversion's fitting
tolerance. Ten island cases and six neighbor-order cases remain covered.
Seven CAD inputs also check Python face/solid/surface counts, design volume
and GLB export. The single-thread Tesla whole-rim checks pass as well.

## Scope and artifacts

This establishes the checked final geometry for `test_corners2`, rather than
only its face count. Full-library literal parity is not established. Other
feature families, true-size conversion, disabled blends and unsupported
transition configurations still use existing recognition/recovery paths.
Eigen and NumPy can choose different tiny intermediate analytic pieces;
on this case they are all replaced inside the same four whole transitions.
Bit-for-bit intermediate arithmetic and STEP serialization are not asserted.

Reviewable files:

- `.cache/corners2-alignment/native-literal-corners2.step`
- `.cache/corners2-alignment/python-literal-corners2.step`
- `tests/fixtures/python_corner_reference.txt`
- `tests/fixtures/python_corner_reference.json`

Reproduce the oracle with the original Python dependencies:

```sh
python tests/python_corner_oracle.py test_corners2.stl tests/fixtures/python_corner_reference.txt --source ../stl2curves
```

Run the native checks:

```sh
ctest --test-dir build_native_occt8_single_thread --output-on-failure
build_native_occt8_single_thread/test_regression build_native_occt8/fixtures
```
