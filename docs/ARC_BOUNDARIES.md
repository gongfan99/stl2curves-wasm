# Continuous fitted boundaries

The Tesla adapter's large rounded front rim is one 78.508 mm B-spline edge in
the Python STEP. The previous native STEP represented it with 40 short spline
edges. The corresponding back rim is one 76.547 mm Python edge and was split
into 39 native edges. These were actual STEP topology divisions, not just
display tessellation. The reference stores these rounded outlines as smooth
B-splines.

Python's `build._runs`, `_split_sharp`, `_run_edge` and `_wire` construct one
edge for a boundary run shared with a neighbor. The native UV-chart fallback
can retain the same smooth geometry in separate mesh-sized edge spans.
`ShapeUpgrade_UnifySameDomain` does not concatenate B-spline spans with its
default settings, so merging the faces alone leaves all those divisions.

`join_boundary_runs` restores the continuous boundary after a reconstruction
has already passed its checks. It copies topology and runs OCCT with edge
unification and C1 spline concatenation enabled, and face merging disabled.
The newer facet-reference path constructs Python's runs directly and skips
this legacy concatenation, preserving the seams checked by the
[corner geometry oracle](CORNERS2_ALIGNMENT.md).
OCCT uses vertex tolerances during spline knot removal. Sewing tolerances can
otherwise allow about 0.014 mm of unintended smoothing on the Tesla rims.
The working copy uses half of the requested tolerance (at most 0.005 mm) for
joining, then restores any larger original vertex tolerances through the
operation's history. The accepted model keeps its original sewing allowances.
It retains the result only when all of the following still hold:

- Valid closed solids, with identical face, solid and surface-type counts.
- The existing volume-preservation and fitted-volume checks pass.
- Maximum shape tolerance stays within the greater of its previous value and
  the configured sewing tolerance.
- Samples on every replaced edge stay within twice the joining tolerance of
  their replacement curve.

This stage runs after reconstruction/recovery in `converter.cpp`. Putting
spline concatenation inside the earlier face-cleanup attempts changed which
attempt passed the volume checks; performing it after selection preserves the
selected fitted surfaces and trims.

## Regression checks

`test_boundaries` forces a quadratic roof through the recovery builder. Before
the correction that six-face roof had 72 edges. After joining its boundary runs
it has the expected 12 edges, with four continuous curved rims. The test also
retains the existing analytic circle, topology, volume and STEP checks.

`tests/fixtures/python_tesla_rims.txt` contains 41 samples on each rounded rim
from the installed Python reference. The source STL SHA-256 is
`001c8843d5f064a324d80597d255b755561d30598f23de319e87f4d3b6462368`.
The optional rim check in `inspect_conversion` requires a single curved edge
to cover every sample within 0.01 mm, in both the retained BRep and a fresh
import of its exported STEP. Projection is restricted to each edge's trim;
sharing an underlying long curve while retaining short edge spans cannot pass.

```sh
build_native_occt8/inspect_conversion part.stl out.step tests/fixtures/python_tesla_rims.txt
```

The corrected Tesla conversion retains 184 faces and one valid solid, without
fallback. Both the retained BRep and exported STEP have one edge per tested
rim. Maximum sampled distances to Python are 0.004235 mm for the front rim and
0.004085 mm for the back rim. The joining pass's sampled displacement across
all replaced edges is 0.004272 mm; the maximum shape tolerance is unchanged.
All six native test suites and seven CAD fixtures with GLB exports pass in
both threaded and single-thread builds.

Face counts need not match Python to establish good reconstruction. Geometry,
validity and continuous boundaries are the relevant checks for this issue.
