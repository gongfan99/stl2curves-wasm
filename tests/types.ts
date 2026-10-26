import createStl2curves, { type ConversionOptions } from '../wasm/stl2curves.js';

async function example(input: Uint8Array) {
  const module = await createStl2curves({ locateFile: (name, prefix) => prefix + name });
  const options: ConversionOptions = { accuracy: 0.1, tolerance: 0.01, trueSize: false };
  const result = module.convert([{ name: 'part.stl', data: input }], options);
  for (const file of result.files) {
    const transferable: ArrayBuffer = file.data.buffer;
    const valid: boolean = file.diagnostics.valid;
    void transferable;
    void valid;
  }
  // @ts-expect-error Files accept binary buffers, not text.
  module.convert([{ name: 'part.stl', data: 'solid text' }]);
  // @ts-expect-error Native handles are not the byte-based public API.
  module.convert([{ name: 'part.stl', data: input }]).delete();
  module.dispose();
}
void example;
