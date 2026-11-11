"""Check packaging into a separate build/deployment workspace."""
import json
import os
import pathlib
import shutil
import subprocess
import tempfile
import unittest

REPO = pathlib.Path(__file__).resolve().parents[1]
CMAKE = os.environ.get('CMAKE_COMMAND') or shutil.which('cmake')


class PackagePaths(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='wasm paths ')
        self.addCleanup(self.temp.cleanup)
        self.root = pathlib.Path(self.temp.name)
        self.source = self.root/'source checkout'
        (self.source/'wasm').mkdir(parents=True)
        for name in ['loader.mjs', 'loader.d.mts']:
            (self.source/'wasm'/name).write_text(name)
        for name in ['LICENSE', 'THIRD_PARTY_NOTICES.md']:
            (self.source/name).write_text(name)

    def build_tree(self, build):
        build.mkdir(parents=True)
        for name in ['stl2curves.js', 'stl2curves.wasm', 'stl2curves.d.ts']:
            (build/name).write_bytes(name.encode())
        notice = build/'vcpkg_installed/wasm32-emscripten-single-thread/share/opencascade/copyright'
        notice.parent.mkdir(parents=True)
        notice.write_text('OCCT copyright')

    def package(self, **options):
        cmd = [CMAKE, f'-DROOT={self.source}', '-DVARIANT=single_thread']
        cmd += [f'-D{k}={v}' for k, v in options.items()]
        subprocess.run(cmd+['-P', str(REPO/'scripts/package_wasm.cmake')], check=True)

    def test_existing_defaults(self):
        self.build_tree(self.source/'build_wasm_single_thread')
        self.package()
        output = self.source/'dist'
        self.assertEqual((output/'single_thread/stl2curves.wasm').read_bytes(), b'stl2curves.wasm')
        self.assertEqual(json.loads((output/'package.json').read_text()), {'type': 'module'})

    def test_external_build_and_output(self):
        build = self.root/'external builds/single'
        output = self.root/'external output/dist'
        self.build_tree(build)
        self.package(BUILD_DIR=build, OUTPUT_DIR=output)
        self.assertFalse((self.source/'dist').exists())
        self.assertTrue((output/'single_thread/stl2curves.wasm').is_file())
        self.assertEqual((output/'licenses/opencascade/copyright').read_text(), 'OCCT copyright')
        self.assertTrue((output/'loader.mjs').is_file())

    def test_ci_variant_package(self):
        build = self.root/'external builds/single'
        output = self.root/'ci artifact'
        self.build_tree(build)
        self.package(BUILD_DIR=build, OUTPUT_DIR=output, VARIANT_ONLY='ON')
        self.assertTrue((output/'single_thread/licenses/opencascade/copyright').is_file())
        self.assertFalse((output/'loader.mjs').exists())


if __name__ == '__main__':
    if not CMAKE:
        raise SystemExit('Set CMAKE_COMMAND or put cmake on PATH')
    unittest.main()
