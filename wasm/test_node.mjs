// Run against built WASM artifacts; optionally include the regression fixtures.
import { readFile } from 'node:fs/promises';
import { pathToFileURL } from 'node:url';
import path from 'node:path';
import { assert, exercise } from './test_helpers.mjs';

const variant = process.argv[2];
if (!['threaded', 'single_thread'].includes(variant)) throw new Error('Specify threaded or single_thread');
const moduleUrl = process.env.STL2CURVES_WASM_MODULE
  ? pathToFileURL(path.resolve(process.env.STL2CURVES_WASM_MODULE))
  : new URL(`../dist/${variant}/stl2curves.js`, import.meta.url);
const { default: factory } = await import(moduleUrl);
const module = await factory({
  wasmBinary: new Uint8Array(await readFile(new URL('stl2curves.wasm', moduleUrl))),
});
try {
  const info = exercise(module, variant === 'threaded');
  console.log(`${variant}: API/ownership/3MF/repair/parallel checks passed (${JSON.stringify(info)})`);
  if (process.env.STL2CURVES_FIXTURES) {
    const cases = JSON.parse(await readFile(new URL('../tests/regression_cases.json', import.meta.url), 'utf8'));
    let observedParallelTasks = 0;
    for (const test of cases) {
      const started = performance.now();
      const data = new Uint8Array(await readFile(path.join(process.env.STL2CURVES_FIXTURES, test.name)));
      const result = module.convert([{ name: test.name, data }], { fuse: false, details: true });
      assert(result.errors.length === 0 && result.files.length === 1, `${test.name}: ${JSON.stringify(result.errors)}`);
      const d = result.files[0].diagnostics;
      assert(d.valid && Math.abs(d.volume - test.volume) <= test.allowedError && d.faces <= test.maxFaces,
        `${test.name}: volume=${d.volume}, faces=${d.faces}, fallback=${d.fallback}`);
      observedParallelTasks += d.parallelTasks;
      console.log(`${variant}: ${test.name} passed: volume=${d.volume}, faces=${d.faces}, fallback=${d.fallback}, parallelTasks=${d.parallelTasks}, seconds=${((performance.now() - started) / 1000).toFixed(2)}`);
    }
    if (variant === 'threaded' && module.buildInfo().workers > 1) assert(observedParallelTasks > 0, 'No parallel fitting tasks observed');
  }
  console.log(`${variant} WASM tests passed`);
} finally { module.dispose(); }
let rejected = false;
try { module.convert([]); } catch { rejected = true; }
assert(rejected, 'Disposed module accepted conversion');
