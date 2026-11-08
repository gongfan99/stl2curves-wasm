"""Record sizing and snapping results from the original Python implementation."""
import argparse
import pathlib
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("--source", type=pathlib.Path)
    args = parser.parse_args()
    if args.source:
        sys.path.insert(0, str(args.source.resolve()))
    import numpy as np
    from stl2curves import features, sizing

    def cylinder(radius, height, origin, tilt=0):
        n = 24
        angle = 2 * np.pi * np.arange(n) / n
        points = np.vstack([np.c_[radius * np.cos(angle), radius * np.sin(angle), np.full(n, z)]
                            for z in (0, height)])
        points = np.vstack([points, [0, 0, 0], [0, 0, height]])
        rotation = np.array([[np.cos(tilt), 0, np.sin(tilt)], [0, 1, 0],
                             [-np.sin(tilt), 0, np.cos(tilt)]])
        points = points @ rotation.T + origin
        tris = []
        for i in range(n):
            j = (i + 1) % n
            tris += [(i, j, j + n), (i, j + n, i + n), (2 * n, j, i), (2 * n + 1, i + n, j + n)]
        return points, np.array(tris)

    def join(meshes):
        points, tris, offset = [], [], 0
        for p, t in meshes:
            points.append(p)
            tris.append(t + offset)
            offset += len(p)
        return np.vstack(points), np.vstack(tris)

    def record(mesh, feature):
        model = feature.model
        kind = model.kind
        if kind == "sphere":
            origin, axis, radius, secondary = model.c, [0, 0, 1], model.r, 0
        else:
            origin, axis = model.a, model.d
            if model.line:
                radius, secondary = model.line
            else:
                radius, z, secondary = model.circle
                origin = origin + z * model.d
        return (kind, np.concatenate([mesh.ftris[i] for i in feature.facets]).tolist(),
                [*origin, *axis, radius, secondary])

    cases = []
    for label, radii, heights, tilt, unit in (
        ("area_weighted_equal_radii", [2, 2.001], [3, 6], 0, ""),
        ("round_group_mm", [.4996, .5002], [2, 4], 0, "mm"),
        ("round_group_inch", [3.174, 3.176], [3, 4], 0, "inch"),
        ("near_axis", [2], [2], .0008, ""),
        ("distinct_radii", [2, 2.02], [3, 4], 0, ""),
    ):
        p, t = join([cylinder(r, h, [10 * i, 0, 0], tilt) for i, (r, h) in enumerate(zip(radii, heights))])
        mesh, patches, _ = features.analyze(p, t)
        before = [record(mesh, f) for f in patches]
        guess = sizing.guess_size([(mesh, patches)])
        features.snap(mesh, patches, unit or None)
        after = [record(mesh, f) for f in patches]
        assert len(before) == len(after), label
        cases.append((label, unit or "-", mesh.pts, mesh.tris, before, after, guess))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as out:
        out.write(f"{len(cases)}\n")
        for name, unit, p, t, before, after, guess in cases:
            out.write(f"{name} {unit} {len(p)} {len(t)} {len(before)} "
                      f"{guess.factor if guess else 1:.17g} {guess.unit if guess else '-'}\n")
            for array in (p, t):
                for row in array:
                    out.write(" ".join(format(float(x), ".17g") for x in row) + "\n")
            for group in (before, after):
                for kind, triangles, parameters in group:
                    out.write(f"{kind} {len(triangles)} " + " ".join(format(float(x), ".17g") for x in parameters) + "\n")
                    out.write(" ".join(map(str, triangles)) + "\n")
    print(f"Generated {len(cases)} sizing and snapping cases")


if __name__ == "__main__":
    main()
