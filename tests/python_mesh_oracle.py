"""Generate mesh parity cases using the original Python package (development only).

Run on a machine with stl2curves installed, or pass --source /path/to/checkout.
The C++ test_python_mesh_parity executable consumes python_mesh_cases.txt.
"""
import argparse
import hashlib
import json
import pathlib
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("--source", type=pathlib.Path)
    parser.add_argument("--stl", nargs="*", type=pathlib.Path, default=[])
    args = parser.parse_args()
    if args.source:
        sys.path.insert(0, str(args.source.resolve()))
    import numpy as np
    import stl2curves
    from stl2curves import repair as module, bodies, simplify

    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                       [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    cube = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4],
                     [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]])
    cases = []

    def add(name, stage, p, t):
        p, t = np.asarray(p, float), np.asarray(t, np.int64).reshape(-1, 3)
        sizes = []
        if stage == "duplicates":
            out = p, module._duplicates(t)[0]
        elif stage == "orient":
            out = p, module._orient(p, t)[0]
        elif stage == "slivers":
            out = p, module.remove_slivers(p, t)[0]
        elif stage == "flaps":
            out = p, module._folded_flaps(p, t)[0]
        elif stage == "knots":
            out = p, module._untie_knots(p, t)[0]
        elif stage == "simplify":
            out = p, simplify.simplify(p, t, .01)
        elif stage == "cracks":
            out = module._close_cracks(p, t, .001)[:2]
        elif stage == "holes":
            out = module._fill_holes(p, t)[:2]
        elif stage == "repair":
            out = module.repair(p, t)[:2]
        elif stage == "bodies":
            groups = bodies.split_bodies(p, t)
            sizes = [len(g) for g in groups]
            out = p, np.vstack(groups) if groups else np.empty((0, 3), int)
        cases.append((name, stage, p, t, *out, sizes))

    for up in range(5):
        for down in range(5):
            if up + down:
                add(f"copies_{up}_{down}", "duplicates", points,
                    [[0, 2, 1]] * up + [[1, 2, 0]] * down)
    for seed in range(12):
        rng = np.random.default_rng(seed)
        t = cube.copy()
        turn = rng.random(len(t)) < .5
        t[turn] = t[turn, ::-1]
        add(f"cube_orientation_{seed}", "orient", points, t)
    add("inside_out", "repair", points, cube[:, ::-1])
    add("first_face_flipped", "repair", points,
        np.vstack([cube[:1, ::-1], cube[1:]]))
    add("cavity", "orient", np.vstack([3 * points, points + 1]),
        np.vstack([cube, cube[:, ::-1] + 8]))
    add("nested_outward_body", "orient", np.vstack([3 * points, points + 1]),
        np.vstack([cube, cube + 8]))
    add("separate_inward_body", "orient", np.vstack([points, points + 3]),
        np.vstack([cube, cube[:, ::-1] + 8]))
    for dropped in ([2], [2, 3], [2, 3, 6], [2, 3, 6, 7]):
        t = np.delete(cube, dropped, axis=0)
        label = "_".join(map(str, dropped))
        add("hole_" + label, "holes", points, t)
        add("repair_hole_" + label, "repair", points, t)
    # A crack wider than load_stl's weld grid but narrower than CRACK.
    p = np.vstack([points, points[[4, 5, 6, 7]] + [0, 0, .0005]])
    t = cube.copy()
    t[2:4] += 4
    add("cap_crack", "cracks", p, t)
    add("repair_cap_crack", "repair", p, t)
    p = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [1, 1, 0]], float)
    add("split_t_junction", "slivers", p, [[1, 0, 2], [2, 0, 3]])
    add("folded_slivers", "slivers", p, [[1, 0, 2], [2, 0, 1]])
    add("same_direction_flap", "flaps", np.vstack([p, [1, -.1, 0]]),
        [[0, 2, 3], [2, 0, 1], [0, 2, 4]])
    # The old native implementation only filled planar holes <= 128 vertices.
    for n in (64, 160):
        angle = 2 * np.pi * np.arange(n) / n
        p = np.vstack([np.c_[np.cos(angle), np.sin(angle), np.ones(n)], [0, 0, 0]])
        t = np.array([[n, (i + 1) % n, i] for i in range(n)])
        add(f"large_hole_{n}", "holes", p, t)
    add("body_cube", "bodies", points, cube)
    for label, shift in (("apart", [3, 0, 0]), ("face_touch", [1, 0, 0]),
                         ("edge_touch", [1, 1, 0]), ("vertex_touch", [1, 1, 1]),
                         ("overlap", [.5, .5, .5])):
        p = np.vstack([points, points + shift])
        t = np.vstack([cube, cube + 8])
        # The reader welds coincident corners before splitting bodies.
        _, inv = np.unique(p, axis=0, return_inverse=True)
        p = np.unique(p, axis=0)
        t = inv[t]
        add("bodies_" + label, "bodies", p, t)
    add("bodies_cavity", "bodies", np.vstack([3 * points, points + 1]),
        np.vstack([cube, cube[:, ::-1] + 8]))
    for distance in (.0005, .001, .003):
        p = np.vstack([points, points[0] + [distance, 0, 0]])
        t = np.vstack([cube, [[0, 8, 1], [8, 0, 3], [0, 8, 4]]])
        add(f"tiny_nonmanifold_knot_{distance}", "knots", p, t)
    for n in (4, 7):
        p, t, where = [], [], {}

        def vertex(q):
            key = tuple(q)
            if key not in where:
                where[key] = len(p)
                p.append(q)
            return where[key]

        # Tessellate every face of the unit cube, welding the integer grid
        # before scaling so exact planar collapse costs have stable ties.
        for axis in range(3):
            u, v = (axis + 1) % 3, (axis + 2) % 3
            for side in (0, n):
                for i in range(n):
                    for j in range(n):
                        corners = []
                        for x, y in ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)):
                            q = [0, 0, 0]
                            q[axis], q[u], q[v] = side, x, y
                            corners.append(vertex(q))
                        a, b, c, d = corners
                        pair = [[a, b, c], [a, c, d]]
                        t.extend(pair if side else [row[::-1] for row in pair])
        add(f"quadric_cube_{n}", "simplify", np.array(p, float) / n, t)
    n = 9
    p = np.array([[2 * i / (n - 1), j / (n - 1), 0]
                  for i in range(n) for j in range(n)])
    p[:, 2] = .03 * np.sin(1.7 * p[:, 0]) + .02 * np.cos(2.3 * p[:, 1]) + .04 * p[:, 0] * p[:, 1]
    t = []
    for i in range(n - 1):
        for j in range(n - 1):
            a, b, c, d = i * n + j, (i + 1) * n + j, (i + 1) * n + j + 1, i * n + j + 1
            t += [[a, b, c], [a, c, d]]
    add("quadric_smooth_graph", "simplify", p, t)
    from stl2curves.features import load_stl
    for path in args.stl:
        p, t = load_stl(path)
        add(path.stem + "_repair", "repair", p, t)
        add(path.stem + "_bodies", "bodies", p, t)
        add(path.stem + "_knots", "knots", p, t)

    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "python_mesh_cases.txt").open("w") as stream:
        stream.write(f"{len(cases)}\n")
        for name, stage, p, t, q, u, sizes in cases:
            stream.write(f"{name} {stage} {len(p)} {len(t)} {len(q)} {len(u)} {len(sizes)}\n")
            for array in (p, t, q, u):
                for row in array:
                    stream.write(" ".join(format(float(x), ".17g") for x in row) + "\n")
            if sizes:
                stream.write(" ".join(map(str, sizes)) + "\n")
    package = pathlib.Path(stl2curves.__file__).parent
    manifest = {"version": stl2curves.__version__, "cases": len(cases), "modules": {
        p.name: hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        for p in sorted(package.glob("*.py"))}}
    (args.output / "python_mesh_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Generated {len(cases)} reference cases from {package}")


if __name__ == "__main__":
    main()
