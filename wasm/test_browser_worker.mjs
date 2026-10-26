import { exercise } from './test_helpers.mjs';
self.onmessage = async ({ data }) => {
  let module;
  try {
    ({ stl2curves: module } = await import('../dist/loader.mjs'));
    const info = exercise(module, data.expected);
    module.dispose();
    self.postMessage({ ok: true, info, isolated: self.crossOriginIsolated });
  } catch (error) {
    module?.dispose();
    self.postMessage({ ok: false, error: String(error?.stack || error) });
  }
};
