export function cubeStl(open = false) {
  const points = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]];
  const triangles = [[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]];
  const kept = open ? triangles.filter((_, i) => i !== 2 && i !== 3) : triangles;
  return new TextEncoder().encode('solid cube\n' + kept.map(t =>
    'facet normal 0 0 0\nouter loop\n' + t.map(v => 'vertex ' + points[v].join(' ') + '\n').join('') + 'endloop\nendfacet\n'
  ).join('') + 'endsolid cube\n');
}
export function assert(condition, message) { if (!condition) throw new Error(message); }
function storedZip(entries) {
  const encoder = new TextEncoder(), locals = [], central = [];
  let offset = 0;
  for (const [name, text] of entries) {
    const filename = encoder.encode(name), data = encoder.encode(text);
    let crc = 0xffffffff;
    for (const byte of data) {
      crc ^= byte;
      for (let bit = 0; bit < 8; bit++) crc = (crc >>> 1) ^ (crc & 1 ? 0xedb88320 : 0);
    }
    crc = (crc ^ 0xffffffff) >>> 0;
    const local = new Uint8Array(30 + filename.length + data.length), l = new DataView(local.buffer);
    l.setUint32(0, 0x04034b50, true); l.setUint16(4, 20, true); l.setUint16(12, 33, true);
    l.setUint32(14, crc, true); l.setUint32(18, data.length, true); l.setUint32(22, data.length, true);
    l.setUint16(26, filename.length, true); local.set(filename, 30); local.set(data, 30 + filename.length);
    const directory = new Uint8Array(46 + filename.length), d = new DataView(directory.buffer);
    d.setUint32(0, 0x02014b50, true); d.setUint16(4, 20, true); d.setUint16(6, 20, true);
    d.setUint16(14, 33, true); d.setUint32(16, crc, true);
    d.setUint32(20, data.length, true); d.setUint32(24, data.length, true);
    d.setUint16(28, filename.length, true); d.setUint32(42, offset, true); directory.set(filename, 46);
    locals.push(local); central.push(directory); offset += local.length;
  }
  const directoryLength = central.reduce((sum, data) => sum + data.length, 0);
  const output = new Uint8Array(offset + directoryLength + 22);
  let position = 0;
  for (const data of [...locals, ...central]) { output.set(data, position); position += data.length; }
  const end = new DataView(output.buffer, position);
  end.setUint32(0, 0x06054b50, true); end.setUint16(8, entries.length, true); end.setUint16(10, entries.length, true);
  end.setUint32(12, directoryLength, true); end.setUint32(16, offset, true);
  return output;
}
export function inchCube3mf() {
  const points = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]];
  const triangles = [[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]];
  const model = `<model unit="inch" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02"><resources><object id="1" name="inch cube"><mesh><vertices>${
    points.map(([x,y,z]) => `<vertex x="${x}" y="${y}" z="${z}"/>`).join('')
  }</vertices><triangles>${triangles.map(([a,b,c]) => `<triangle v1="${a}" v2="${b}" v3="${c}"/>`).join('')
  }</triangles></mesh></object></resources><build><item objectid="1" transform="1 0 0 0 1 0 0 0 1 2 3 4"/></build></model>`;
  return storedZip([
    ['[Content_Types].xml', '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>'],
    ['_rels/.rels', '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="model" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel" Target="/3D/3dmodel.model"/></Relationships>'],
    ['3D/3dmodel.model', model],
  ]);
}
function cylindersStl(count = 12) {
  const facets = [];
  const around = 32;
  const facet = (a, b, c) => facets.push('facet normal 0 0 0\nouter loop\n' +
    [a, b, c].map(p => `vertex ${p.join(' ')}\n`).join('') + 'endloop\nendfacet\n');
  for (let body = 0; body < count; body++) {
    const ring = (i, z) => [3 * body + Math.cos(2 * Math.PI * i / around), Math.sin(2 * Math.PI * i / around), z];
    for (let i = 0; i < around; i++) {
      const a = ring(i, 0), b = ring(i + 1, 0), c = ring(i, 1), d = ring(i + 1, 1);
      facet(a, b, d); facet(a, d, c);
      facet([3 * body, 0, 0], b, a); facet([3 * body, 0, 1], c, d);
    }
  }
  return new TextEncoder().encode('solid cylinders\n' + facets.join('') + 'endsolid cylinders\n');
}
export function exercise(module, expectedThreaded) {
  const info = module.buildInfo();
  assert(info.threaded === expectedThreaded, 'Incorrect threading variant');
  assert(/^8\.0\.1/.test(info.occtVersion), 'Expected pinned OpenCascade 8.0.1');
  const bytes = cubeStl();
  const first = module.convert([{ name: 'cube.stl', data: bytes }], { details: true });
  assert(first.errors.length === 0, JSON.stringify(first.errors));
  assert(first.files.length === 1, 'Expected one STEP output');
  const file = first.files[0];
  assert(file.data instanceof Uint8Array, 'Output must be a Uint8Array');
  assert(file.data.byteOffset === 0 && file.data.buffer.byteLength === file.data.byteLength,
    'Output must own an independent, exactly sized buffer');
  assert(file.diagnostics.valid && file.diagnostics.faces === 6, 'Expected valid six-face cube');
  assert(Math.abs(file.diagnostics.volume - 1) < 1e-7, 'Incorrect volume');
  assert(new TextDecoder().decode(file.data.subarray(0, 100)).includes('ISO-10303-21'), 'Invalid STEP header');
  const retained = file.data.slice();
  for (let i = 0; i < 3; i++) {
    const next = module.convert([{ name: 'a.stl', data: bytes }, { name: 'bad.stl', data: new Uint8Array([1,2,3]) }], {
      curves: false, blends: false, repair: false, simplify: 0, merge: 'combined',
    });
    assert(next.files.length === 2 && next.errors.length === 1, 'Partial batch success/merge');
  }
  assert(retained.every((v, i) => v === file.data[i]), 'Output aliases freed WASM memory');
  const projectBytes = inchCube3mf();
  const project = module.convert([{ name: 'assembly.3mf', data: projectBytes },
    { name: 'truncated.3mf', data: projectBytes.subarray(0, 9) }]);
  assert(project.errors.length === 1 && project.files.length === 1, '3MF archive failure/partial success');
  assert(project.files[0].name === 'assembly - inch cube.step', '3MF object name was lost');
  assert(project.files[0].diagnostics.valid && project.files[0].diagnostics.faces === 6 &&
    Math.abs(project.files[0].diagnostics.volume - 25.4 ** 3) < 1e-6, '3MF units or STEP geometry');
  const failed = module.convert([{ name: 'open.stl', data: cubeStl(true) }], { repair: false, curves: false });
  assert(failed.files.length === 0 && failed.errors.length === 1 && failed.errors[0].diagnostics?.openEdges === 4,
    'Failed-object diagnostics lost unresolved mesh defects');
  const parallel = module.convert([{ name: 'cylinders.stl', data: cylindersStl() }], { fuse: false });
  assert(parallel.errors.length === 0 && parallel.files.length === 1, JSON.stringify(parallel.errors));
  const diagnostics = parallel.files[0].diagnostics;
  assert(diagnostics.valid && diagnostics.solids === 12, 'Disconnected analytic bodies were lost');
  if (expectedThreaded && info.workers > 1) assert(diagnostics.parallelTasks > 0, 'No fitting tasks ran on pthreads');
  else assert(diagnostics.parallelTasks === 0, 'Sequential configuration ran background tasks');
  let rejected = false;
  try { module.convert([], { tolerance: 0 }); } catch { rejected = true; }
  assert(rejected, 'Invalid options must fail');
  if (!expectedThreaded) assert(file.diagnostics.parallelTasks === 0, 'Single-thread build ran parallel tasks');
  return info;
}
