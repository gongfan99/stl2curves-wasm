// Browser loader for blob-based deployment of either WASM variant.
const canUseThreads =
  typeof globalThis.SharedArrayBuffer !== 'undefined' &&
  globalThis.crossOriginIsolated === true;
const assetDirectory = canUseThreads ? './threaded/' : './single_thread/';
const assetBase = new URL(assetDirectory, import.meta.url);
const urls = [];

async function fetchBlob(name) {
  const response = await fetch(new URL(name, assetBase));
  if (!response.ok) throw new Error(`Cannot load ${name}: HTTP ${response.status}`);
  return response.blob();
}

let stl2curves;
try {
  const [jsBlob, wasmBlob] = await Promise.all([
    fetchBlob('stl2curves.js'), fetchBlob('stl2curves.wasm'),
  ]);
  const urlJs = URL.createObjectURL(new Blob([jsBlob], { type: 'text/javascript' }));
  const urlWasm = URL.createObjectURL(wasmBlob);
  urls.push(urlJs, urlWasm);
  const { default: factory } = await import(urlJs);
  stl2curves = await factory({
    locateFile(path) {
      if (path.endsWith('.wasm')) return urlWasm;
      // Some SDK versions emit a separate pthread bootstrap script. Packaging
      // includes it, and its URL is resolved against the real asset directory.
      if (/\.worker\.(m?js)$/.test(path)) return new URL(path, assetBase).href;
      throw new Error(`Unknown stl2curves asset: ${path}`);
    },
    mainScriptUrlOrBlob: urlJs,
  });
  const dispose = stl2curves.dispose.bind(stl2curves);
  stl2curves.dispose = () => {
    dispose();
    for (const url of urls.splice(0)) URL.revokeObjectURL(url);
  };
} catch (error) {
  for (const url of urls) URL.revokeObjectURL(url);
  throw error;
}

export { stl2curves };
export default stl2curves;
