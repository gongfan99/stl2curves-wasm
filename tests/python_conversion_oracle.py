"""Export comparable Python conversion results and geometry for parity checks.

This is a development oracle, never a runtime dependency of the C++ library.
"""
import argparse
import collections
import hashlib
import json
import pathlib
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("inputs", nargs="+", type=pathlib.Path)
    parser.add_argument("--source", type=pathlib.Path)
    parser.add_argument("--accuracy", type=float, default=.1)
    parser.add_argument("--time-limit", type=float, default=600)
    parser.add_argument("--no-fuse", action="store_true")
    args = parser.parse_args()
    if args.source:
        sys.path.insert(0, str(args.source.resolve()))
    import stl2curves
    from stl2curves import accuracy, convert, workers
    from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepTools import BRepTools
    from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_SOLID
    from OCP.TopExp import TopExp
    from OCP.TopoDS import TopoDS

    # No process scheduling differences between successive oracle runs.
    workers.WORKERS = 1
    accuracy.set_accuracy(args.accuracy)
    package = pathlib.Path(stl2curves.__file__).parent
    metadata = {"version": stl2curves.__version__, "accuracy": args.accuracy,
                "fuse": not args.no_fuse, "modules": {
                    p.name: hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
                    for p in sorted(package.glob("*.py"))}}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")

    def unique(shape, kind):
        # OCP 8 uses the current NCollection map; use the helper supplied by
        # stl2curves for compatibility with its supported wheel versions.
        from stl2curves.build import IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher
        index = IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher()
        TopExp.MapShapes_s(shape, kind, index)
        return [index.FindKey(i) for i in range(1, index.Extent() + 1)]

    for path in args.inputs:
        name = path.stem
        start = time.monotonic()
        print("CONVERT", path, flush=True)
        try:
            shape, info = convert.stl_to_solid(str(path), .01, fuse=not args.no_fuse,
                                              time_limit=args.time_limit or None)
            convert.write_step(shape, args.output / (name + ".step"))
            BRepTools.Write_s(shape, str(args.output / (name + ".brep")))
            surfaces = collections.Counter(BRepAdaptor_Surface(TopoDS.Face(f)).GetType().name
                                           for f in unique(shape, TopAbs_FACE))
            curves = collections.Counter(BRepAdaptor_Curve(TopoDS.Edge(e)).GetType().name
                                         for e in unique(shape, TopAbs_EDGE))
            report = {"input": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                      "valid": BRepCheck_Analyzer(shape, True, True).IsValid(),
                      "faces": info["faces"], "solids": convert.count(shape, TopAbs_SOLID),
                      "volume": convert.volume(shape), "surfaces": dict(surfaces),
                      "curves": dict(curves), "triangles": info["triangles"],
                      "repairs": info["repairs"], "mesh_left": info.get("mesh_left", []),
                      "restored": dict(collections.Counter(f.model.kind for f in info["restored"])),
                      "skipped": dict(collections.Counter(f.model.kind for f in info["skipped"]))}
        except Exception as error:
            report = {"input": path.name, "error": str(error)}
        report["seconds"] = time.monotonic() - start
        (args.output / (name + ".json")).write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report), flush=True)
    workers.stop()


if __name__ == "__main__":
    main()
