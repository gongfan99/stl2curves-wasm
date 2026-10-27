import test from 'node:test';
import assert from 'node:assert/strict';
import { mkdtemp, mkdir, writeFile, readFile, readdir, cp, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { spawnSync } from 'node:child_process';

const repo = fileURLToPath(new URL('../', import.meta.url));
const cmake = process.env.CMAKE_COMMAND || 'cmake';
const ports = ['opencascade', 'eigen3', 'nanoflann', 'libzip', 'pugixml', 'zlib'];
const variants = ['threaded', 'single_thread'];
const types = 'export interface Stl2curvesModule {}\n';

function runScript(name, definitions) {
  const result = spawnSync(cmake, [
    ...Object.entries(definitions).map(([key, value]) => `-D${key}=${value}`),
    '-P', path.join(repo, 'scripts', name),
  ], { encoding: 'utf8' });
  if (result.error) throw result.error;
  return result;
}
function success(result) {
  assert.equal(result.status, 0, result.stderr || result.stdout);
}
async function put(root, name, data) {
  const target = path.join(root, name);
  await mkdir(path.dirname(target), { recursive: true });
  await writeFile(target, data);
}
async function fixture(t) {
  const root = await mkdtemp(path.join(tmpdir(), 'stl2curves-package-'));
  t.after(async () => {
    assert.equal(path.dirname(path.resolve(root)), path.resolve(tmpdir()));
    assert.ok(path.basename(root).startsWith('stl2curves-package-'));
    await rm(root, { recursive: true, force: true });
  });
  await mkdir(path.join(root, 'wasm'));
  for (const name of ['loader.mjs', 'loader.d.mts']) {
    await cp(path.join(repo, 'wasm', name), path.join(root, 'wasm', name));
  }
  for (const name of ['LICENSE', 'THIRD_PARTY_NOTICES.md']) {
    await cp(path.join(repo, name), path.join(root, name));
  }
  for (const variant of variants) {
    const build = `build_wasm_${variant}`;
    await put(root, `${build}/stl2curves.js`, `// ${variant} module\n`);
    await put(root, `${build}/stl2curves.wasm`, Buffer.from([0, 97, 115, 109, 1, 0, 0, 0]));
    await put(root, `${build}/stl2curves.d.ts`, types);
    await put(root, `${build}/stl2curves.worker.mjs`, `// ${variant} worker\n`);
    const triplet = `wasm32-emscripten-${variant.replaceAll('_', '-')}`;
    for (const port of ports) {
      await put(root, `${build}/vcpkg_installed/${triplet}/share/${port}/copyright`, `${port} notice\n`);
    }
    success(runScript('package_wasm.cmake', { ROOT: root, VARIANT: variant, VARIANT_ONLY: 'ON' }));
    await cp(path.join(root, 'dist', variant), path.join(root, 'parts', variant), { recursive: true });
  }
  return { root, input: path.join(root, 'parts'), output: path.join(root, 'final') };
}
function assemble(f) {
  return runScript('assemble_wasm.cmake', { ROOT: f.root, INPUT_DIR: f.input, OUTPUT_DIR: f.output });
}

test('variant artifacts contain no root shared files; assembly adds them once', async t => {
  const f = await fixture(t);
  assert.deepEqual((await readdir(path.join(f.root, 'dist'))).sort(), [...variants].sort());
  for (const variant of variants) {
    assert.deepEqual((await readdir(path.join(f.input, variant))).sort(), [
      'licenses', 'stl2curves.d.ts', 'stl2curves.js', 'stl2curves.wasm', 'stl2curves.worker.mjs',
    ]);
  }
  success(assemble(f));
  assert.deepEqual((await readdir(f.output)).sort(), [
    'LICENSE', 'THIRD_PARTY_NOTICES.md', 'licenses', 'loader.d.mts', 'loader.mjs',
    'package.json', 'single_thread', 'threaded',
  ]);
  assert.deepEqual(JSON.parse(await readFile(path.join(f.output, 'package.json'), 'utf8')), { type: 'module' });
  assert.equal(await readFile(path.join(f.output, 'loader.mjs'), 'utf8'),
    await readFile(path.join(repo, 'wasm/loader.mjs'), 'utf8'));
  assert.deepEqual((await readdir(path.join(f.output, 'licenses'))).sort(), [...ports].sort());
  for (const variant of variants) {
    assert.equal(await readFile(path.join(f.output, variant, 'stl2curves.worker.mjs'), 'utf8'), `// ${variant} worker\n`);
    assert.ok(!(await readdir(path.join(f.output, variant))).includes('licenses'));
  }
});

test('conflicting dependency notices fail assembly', async t => {
  const f = await fixture(t);
  await put(f.input, 'single_thread/licenses/opencascade/copyright', 'different notice\n');
  const result = assemble(f);
  assert.notEqual(result.status, 0);
  assert.match(result.stderr, /Dependency notices differ between variants: opencascade/);
});

test('different type declarations fail assembly', async t => {
  const f = await fixture(t);
  await put(f.input, 'single_thread/stl2curves.d.ts', 'export interface OtherModule {}\n');
  const result = assemble(f);
  assert.notEqual(result.status, 0);
  assert.match(result.stderr, /Type declarations differ between variants/);
});

test('missing or empty variant assets fail assembly', async t => {
  const f = await fixture(t);
  await rm(path.join(f.input, 'single_thread/stl2curves.wasm'));
  let result = assemble(f);
  assert.notEqual(result.status, 0);
  assert.match(result.stderr, /Missing single_thread build output: stl2curves.wasm/);
  f.output = path.join(f.root, 'second-final');
  await put(f.input, 'single_thread/stl2curves.wasm', '');
  result = assemble(f);
  assert.notEqual(result.status, 0);
  assert.match(result.stderr, /Empty single_thread build output: stl2curves.wasm/);
});

test('missing dependency notices fail assembly', async t => {
  const f = await fixture(t);
  await rm(path.join(f.input, 'threaded/licenses/libzip/copyright'));
  const result = assemble(f);
  assert.notEqual(result.status, 0);
  assert.match(result.stderr, /Missing threaded dependency notice: libzip/);
});

test('assembly rejects an output directory with existing files', async t => {
  const f = await fixture(t);
  await put(f.output, 'existing.txt', 'preserve me');
  const result = assemble(f);
  assert.notEqual(result.status, 0);
  assert.match(result.stderr, /Assembly output directory must be empty/);
  assert.equal(await readFile(path.join(f.output, 'existing.txt'), 'utf8'), 'preserve me');
});

test('local packaging still adds shared files and root dependency notices', async t => {
  const f = await fixture(t);
  success(runScript('package_wasm.cmake', { ROOT: f.root, VARIANT: 'single_thread' }));
  assert.equal(await readFile(path.join(f.root, 'dist/loader.mjs'), 'utf8'),
    await readFile(path.join(repo, 'wasm/loader.mjs'), 'utf8'));
  assert.equal(await readFile(path.join(f.root, 'dist/licenses/libzip/copyright'), 'utf8'), 'libzip notice\n');
});
