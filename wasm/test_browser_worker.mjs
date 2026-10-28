import { exercise, exercisePicking } from './test_helpers.mjs';
self.onmessage = async ({ data }) => {
  let module;
  try {
    ({ stl2curves: module } = await import('../dist/loader.mjs'));
    const info = exercise(module, data.expected);
    const { GLTFLoader } = await import('../node_modules/three/examples/jsm/loaders/GLTFLoader.js');
    const THREE = await import('../node_modules/three/build/three.module.js');
    await exercisePicking(module, new GLTFLoader(), THREE);
    module.dispose();
    self.postMessage({ ok: true, info, isolated: self.crossOriginIsolated });
  } catch (error) {
    module?.dispose();
    self.postMessage({ ok: false, error: String(error?.stack || error) });
  }
};
