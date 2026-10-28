import { readFile } from 'node:fs/promises';
import vm from 'node:vm';
import test from 'node:test';
import assert from 'node:assert/strict';

const source = await readFile(new URL('../wasm/post.js', import.meta.url), 'utf8');

test('post-js wraps an initialized module without SDK callback helpers', () => {
  let observed;
  const module = { convert(inputs, options) { observed = { inputs, options }; return 'result'; } };
  vm.runInNewContext(source, { Module: module });
  const inputs = [];
  assert.equal(module.convert(inputs), 'result');
  assert.equal(observed.inputs, inputs);
  assert.equal(Object.keys(observed.options).length, 0);
  const options = { curves: false };
  module.convert(inputs, options);
  assert.equal(observed.options, options);
});

test('conversion errors release the busy guard and disposal terminates the pool once', () => {
  let calls = 0, terminations = 0;
  const module = { convert() {
    if (++calls === 1) throw new Error('conversion failed');
    assert.throws(() => module.convert([]), /already running/);
    assert.throws(() => module.dispose(), /cancel an active model operation/);
    return 'recovered';
  } };
  vm.runInNewContext(source, {
    Module: module, PThread: { terminateAllThreads() { terminations++; } },
  });
  assert.throws(() => module.convert([]), /conversion failed/);
  assert.equal(module.convert([]), 'recovered');
  module.dispose();
  module.dispose();
  assert.equal(terminations, 1);
  assert.throws(() => module.convert([]), /disposed/);
});

test('pthread bootstrap does not require application Embind exports', () => {
  const module = {};
  vm.runInNewContext(source, { Module: module });
  assert.equal(module.convert, undefined);
  assert.equal(typeof module.dispose, 'function');
  module.dispose();
});

test('all model APIs reject use after disposal and disposal releases models first', () => {
  const events = [];
  const module = { _clearModels() { events.push('release'); } };
  const methods = ['convert', 'exportStep', 'exportGltf', 'getClosestPoints', 'getShapeProperties', 'releaseModel'];
  for (const name of methods) module[name] = (...args) => args;
  vm.runInNewContext(source, { Module: module, PThread: { terminateAllThreads() { events.push('threads'); } } });
  assert.equal(Object.keys(module.exportGltf(1)[1]).length, 0);
  module.dispose();
  module.dispose();
  assert.deepEqual(events, ['release', 'threads']);
  for (const name of methods) assert.throws(() => module[name](1), /disposed/);
});

test('native error results become JS exceptions after C++ has returned', () => {
  let calls = 0;
  const module = { exportStep() {
    if (calls++ === 0) return { __stl2curvesError: { name: 'TypeError', message: 'Released model' } };
    if (calls === 2) return { __stl2curvesError: { name: 'Error', message: 'STEP failed' } };
    return 'recovered';
  } };
  vm.runInNewContext(source, { Module: module });
  assert.throws(() => module.exportStep(1), error => error.name === 'TypeError' && error.message === 'Released model');
  assert.throws(() => module.exportStep(1), /STEP failed/);
  assert.equal(module.exportStep(1), 'recovered');
});
