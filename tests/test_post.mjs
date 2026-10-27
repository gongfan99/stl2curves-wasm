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
    assert.throws(() => module.dispose(), /cancel an active conversion/);
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
