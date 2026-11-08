"""Export facet grouping and island decisions from the original Python code."""
import argparse
import pathlib
from types import SimpleNamespace

import numpy as np
from stl2curves import blends, features, freeform


def grid(n=8):
    pts = np.array([(x, y, .005*x*x + .007*y*y)
                    for y in range(n) for x in range(n)], float)
    tris, cells = [], []
    for y in range(n-1):
        for x in range(n-1):
            a, b, c, d = y*n+x, y*n+x+1, (y+1)*n+x+1, (y+1)*n+x
            tris += [(a, b, c), (a, c, d)]
            cells += [(x, y), (x, y)]
    return pts, np.array(tris), cells


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=pathlib.Path)
    args = parser.parse_args()
    cases = []
    for name, hole, two, kind in [
        ('one_facet', {(3, 3)}, False, 'freeform'),
        ('four_facets', {(2, 2), (2, 3), (3, 2), (3, 3)}, False, 'freeform'),
        ('five_facets', {(2, 2), (2, 3), (3, 2), (3, 3), (4, 3)}, False, 'freeform'),
        ('two_owners', {(3, 3)}, True, 'freeform'),
        ('open_boundary', {(0, 0)}, False, 'freeform'),
        ('smooth_blend', {(3, 3)}, False, 'blend'),
    ]:
        pts, tris, cells = grid()
        mesh = features.Mesh(pts, tris)
        L = 7.
        quadratic = np.array([0., 0., L*L/3, L*L])
        C = .005*quadratic[:, None] + .007*quadratic[None, :]
        knots = np.array([0.]*4 + [L]*4)
        model = freeform.Freeform('plane', np.zeros(3), np.array([0., 0., 1.]),
                                  np.array([1., 0., 0.]), 0., knots, knots, C)
        owner = np.array([-1 if c in hole else int(two and c[0] >= 4) for c in cells])
        cases.append((name, mesh, owner, kind, model, L))
    octa = np.array([(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
                     (0, 0, 1), (0, 0, -1)], float)
    tris = np.array([(0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
                     (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5)])
    for name, radius, kind in [('sphere_centroid', 1., 'sphere'),
                               ('thin_blend', .01, 'blend'),
                               ('sharp_blend', 1., 'blend')]:
        mesh = features.Mesh(octa*radius, tris)
        owner = np.zeros(8, int)
        owner[0] = -1
        cases.append((name, mesh, owner, kind, features.Sphere(np.zeros(3), radius), radius))
    pts, tris, _ = grid()
    cases.append(('no_owner', features.Mesh(pts, tris), np.full(len(tris), -1),
                  'freeform', None, 7.))

    out = [str(len(cases))]
    for name, mesh, owners, kind, model, size in cases:
        features._mesh_tol = features.TOL + features.FLOAT_TOL * float(np.abs(mesh.pts).max())
        patches, initial = [], []
        for k in range(int(owners.max())+1):
            fids = np.array([f for f, ts in enumerate(mesh.ftris)
                             if np.all(owners[ts] == k)])
            tids = np.concatenate([mesh.ftris[f] for f in fids])
            T = mesh.pts[mesh.tris[tids]]
            mids = ((T + T[:, [1, 2, 0]])/2).reshape(-1, 3)
            worst = float(np.abs(model.signed(mids)).max())
            patches.append(SimpleNamespace(kind='blend' if kind == 'blend' else 'trimmed',
                                           model=model, facets=fids, worst=worst))
            initial.append(tids.tolist())
        grown = blends._islands(mesh, patches)
        expected = [np.concatenate([mesh.ftris[f] for f in p.facets]).tolist() for p in grown]
        out += [f'{name} {kind} {size:.17g}',
                f'{len(mesh.pts)} {len(mesh.tris)} {len(initial)} {len(mesh.ftris)}']
        out += [' '.join(f'{v:.17g}' for v in p) for p in mesh.pts]
        out += [' '.join(map(str, t)) for t in mesh.tris]
        for f in mesh.ftris:
            out.append(str(len(f))+' '+' '.join(map(str, f)))
        for neighbors in mesh.nbrs:
            out.append(str(len(neighbors))+' '+' '.join(map(str, neighbors)))
        for before, after in zip(initial, expected):
            out.append(str(len(before))+' '+' '.join(map(str, before)))
            out.append(str(len(after))+' '+' '.join(map(str, after)))
        print(name, [len(x) for x in initial], '->', [len(x) for x in expected])
    # Mesh.nbrs uses Python integer sets. Exercise resizing and colliding hashes,
    # which the small geometric fixtures above do not all reach.
    orders = [list(range(80)), list(reversed(range(80))),
              [32*i for i in range(100)],
              [1024*i + (i % 7) for i in range(150)],
              [i % 37 for i in range(200)],
              [2147483647 - 64*i for i in range(100)]]
    out.append(str(len(orders)))
    for values in orders:
        expected = set()
        for value in values:
            expected.add(value)
        out.append(str(len(values))+' '+' '.join(map(str, values)))
        out.append(str(len(expected))+' '+' '.join(map(str, expected)))
    args.output.write_text('\n'.join(out)+'\n')


if __name__ == '__main__':
    main()
