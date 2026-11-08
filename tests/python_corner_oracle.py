"""Export complete corner-patch ownership and surface samples from Python."""
import argparse
import pathlib
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=pathlib.Path)
    parser.add_argument('output', type=pathlib.Path)
    parser.add_argument('--source', type=pathlib.Path)
    parser.add_argument('--step', type=pathlib.Path)
    args = parser.parse_args()
    if args.source:
        sys.path.insert(0, str(args.source))
    import numpy as np
    from stl2curves import accuracy, convert, workers
    accuracy.set_accuracy(.1)
    workers.WORKERS = 1
    parts, info, part = convert.prepare(str(args.input), .01, fuse=False)
    assert len(parts) == 1
    mesh, features, _ = parts[0]
    rows = [f'{len(mesh.pts)} {len(mesh.tris)} {len(features)}']
    rows += [' '.join(format(float(x), '.17g') for x in p) for p in mesh.pts]
    rows += [' '.join(map(str, t)) for t in mesh.tris]
    for feature in features:
        model = feature.model
        kind = model.kind
        radius = secondary = 0.
        origin = getattr(model, 'a', getattr(model, 'c', np.zeros(3))).copy()
        axis = getattr(model, 'd', np.array([0., 0., 1.]))
        if kind in ('cylinder', 'cone'):
            radius, secondary = model.line
        elif kind == 'torus':
            radius, height, secondary = model.circle
            origin += height * axis
        elif kind == 'sphere':
            radius = model.r
        triangles = np.concatenate([mesh.ftris[f] for f in feature.facets])
        rows += [f'{kind} {int(feature.kind == "trimmed")}', ' '.join(format(float(x), '.17g') for x in
                               [radius, secondary, *origin, *axis]),
                 str(len(triangles))+' '+' '.join(map(str, triangles))]
        samples = []
        if kind == 'freeform':
            for u in np.linspace(model.tu[3], model.tu[-4], 9):
                for v in np.linspace(model.tv[3], model.tv[-4], 9):
                    p = model.point([u], [v])[0]
                    n = model.normal([p])[0]
                    samples.append([*p, *n])
        rows.append(str(len(samples)))
        rows += [' '.join(format(float(x), '.17g') for x in p) for p in samples]
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.TopAbs import TopAbs_EDGE
    from OCP.TopExp import TopExp
    from OCP.TopoDS import TopoDS
    from OCP.collections import IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher
    shape, _ = convert.finish(parts, info, part, .01)
    edges = IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher()
    TopExp.MapShapes_s(shape, TopAbs_EDGE, edges)
    rows.append(str(edges.Extent()))
    for i in range(1, edges.Extent()+1):
        curve = BRepAdaptor_Curve(TopoDS.Edge(edges(i)))
        rows.append(f'{int(curve.GetType())} 17')
        for t in np.linspace(curve.FirstParameter(), curve.LastParameter(), 17):
            p = curve.Value(float(t))
            rows.append(' '.join(format(float(x), '.17g') for x in (p.X(), p.Y(), p.Z())))
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props, convert.VOLUME_EPS)
    rows.append(format(props.Mass(), '.17g'))
    rows.append(format(convert.VOLUME_EPS, '.17g'))
    args.output.write_text('\n'.join(rows)+'\n')
    import hashlib
    import importlib.metadata
    import json
    import stl2curves
    package = pathlib.Path(stl2curves.__file__).parent
    modules = ['features', 'fillets', 'blends', 'freeform', 'build', 'convert', 'accuracy']
    manifest = {'input': args.input.name,
                'sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
                'version': importlib.metadata.version('stl2curves'),
                'tolerance': .01, 'accuracy': .1, 'fuse': False, 'workers': 1,
                'points': len(mesh.pts), 'triangles': len(mesh.tris),
                'features': len(features), 'edges': edges.Extent(),
                'volume_epsilon': convert.VOLUME_EPS, 'brep_volume': props.Mass(),
                'module_hashes': {name: hashlib.sha256((package/(name+'.py')).read_bytes()
                                                      .replace(b'\r\n', b'\n')).hexdigest()
                                  for name in modules}}
    args.output.with_suffix('.json').write_text(json.dumps(manifest, indent=2)+'\n')
    if args.step:
        from OCP.BRepTools import BRepTools
        BRepTools.Write_s(shape, str(args.step.with_suffix('.brep')))
        convert.write_step(shape, str(args.step))
    print('reference', len(features), [(f.model.kind, len(f.facets)) for f in features])
    workers.stop()


if __name__ == '__main__':
    main()
