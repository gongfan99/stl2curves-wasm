export function cubeStl(open = false, offsetX = 0) {
  const points = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]];
  for (const point of points) point[0] += offsetX;
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
export function parseGlb(bytes) {
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  assert(view.getUint32(0, true) === 0x46546c67 && view.getUint32(4, true) === 2 &&
    view.getUint32(8, true) === bytes.byteLength, 'Invalid GLB header');
  assert(view.getUint32(16, true) === 0x4e4f534a, 'Missing GLB JSON');
  return JSON.parse(new TextDecoder().decode(bytes.subarray(20, 20 + view.getUint32(12, true))));
}
export function exercise(module, expectedThreaded) {
  const info = module.buildInfo();
  assert(info.version === '0.2.0' && info.threaded === expectedThreaded, 'Incorrect API version/variant');
  assert(/^8\.0\.1/.test(info.occtVersion), 'Expected pinned OpenCascade 8.0.1');
  const owned = [];
  const convert = (inputs, options) => {
    const result = module.convert(inputs, options);
    owned.push(...result.models.map(model => model.modelId));
    assert(result.files === undefined, 'convert still returns serialized files');
    return result;
  };
  try {
    const bytes = cubeStl();
    const first = convert([{ name: 'cube.stl', data: bytes }], { details: true });
    assert(first.errors.length === 0 && first.models.length === 1, JSON.stringify(first.errors));
    const model = first.models[0];
    assert(model.name === 'cube' && model.root.modelId === model.modelId, 'Model metadata');
    assert(model.diagnostics.valid && model.diagnostics.faces === 6, 'Expected valid six-face model');
    const snapshot = JSON.stringify(model.diagnostics);
    const glb = module.exportGltf(model.modelId);
    const parsed = parseGlb(glb.data);
    const file = module.exportStep(model.modelId);
    for (const output of [file, glb]) {
      assert(output.data instanceof Uint8Array && output.data.byteOffset === 0 &&
        output.data.buffer.byteLength === output.data.byteLength, 'Output must own its exact-sized buffer');
    }
    assert(file.diagnostics.valid && file.diagnostics.faces === 6 &&
      Math.abs(file.diagnostics.volume - 1) < 1e-7, 'Invalid STEP cube');
    assert(new TextDecoder().decode(file.data.subarray(0, 100)).includes('ISO-10303-21'), 'Invalid STEP header');
    const retainedStep = file.data.slice(), retainedGlb = glb.data.slice();
    const second = parseGlb(module.exportGltf(model.modelId).data);
    assert(JSON.stringify(parsed.nodes) === JSON.stringify(second.nodes), 'Export changed graph identities');
    assert(JSON.stringify(model.diagnostics) === snapshot, 'Export changed model diagnostics');
    const counts = { 2: 0, 4: 0, 5: 0 };
    let selectedEdge;
    for (const node of parsed.nodes) {
      assert(node.extras.modelId === model.modelId && Array.isArray(node.extras.usagePath), 'Missing node selection');
      module.getShapeProperties(node.extras);
      if (node.mesh !== undefined) counts[node.extras.brepgraph_id.kind]++;
      if (node.extras.brepgraph_id.kind === 4) selectedEdge = node.extras;
    }
    assert(counts[2] === 6 && counts[4] === 12 && counts[5] === 8, 'Incorrect cube selection counts');
    assert(module.getShapeProperties(selectedEdge).properties.some(p => p.name === 'length' && Math.abs(p.value[0]-1)<1e-8), 'Exact edge properties');
    const moved = convert([{name:'moved.stl',data:cubeStl(false,5)}]).models[0];
    const distance = module.getClosestPoints(model.root,moved.root);
    assert(Math.abs(distance.distance-4)<1e-7 && distance.p1.length===3 && distance.p2.length===3, 'Cross-model exact distance');
    const invalid = structuredClone(selectedEdge);
    invalid.usagePath.at(-1).node.index += 100;
    let rejected = false;
    try { module.getShapeProperties(invalid); } catch { rejected = true; }
    assert(rejected, 'Invalid usage path accepted');
    for (let i = 0; i < 3; i++) {
      const next = convert([{ name: 'a.stl', data: bytes }, { name: 'bad.stl', data: new Uint8Array([1,2,3]) }], {
        curves: false, blends: false, repair: false, simplify: 0, merge: 'combined',
      });
      assert(next.models.length === 2 && next.errors.length === 1, 'Partial batch success/combined model');
      for (const item of next.models) module.exportStep(item.modelId);
    }
    assert(retainedStep.every((v,i)=>v===file.data[i]) && retainedGlb.every((v,i)=>v===glb.data[i]), 'Output aliases freed WASM memory');
    const projectBytes=inchCube3mf();
    const project=convert([{name:'assembly.3mf',data:projectBytes},{name:'truncated.3mf',data:projectBytes.subarray(0,9)}]);
    assert(project.errors.length===1 && project.models.length===1, '3MF partial success');
    assert(project.models[0].name==='assembly - inch cube', '3MF name lost');
    const projectStep=module.exportStep(project.models[0].modelId);
    assert(projectStep.diagnostics.valid && Math.abs(projectStep.diagnostics.volume-25.4**3)<1e-6, '3MF units/STEP validation');
    const failed=convert([{name:'open.stl',data:cubeStl(true)}],{repair:false,curves:false});
    assert(failed.models.length===0 && failed.errors.length===1 && failed.errors[0].diagnostics?.openEdges===4,'Mesh failure diagnostics');
    const parallel=convert([{name:'cylinders.stl',data:cylindersStl()}],{fuse:false});
    assert(parallel.errors.length===0 && parallel.models.length===1,JSON.stringify(parallel.errors));
    const diagnostics=parallel.models[0].diagnostics;
    assert(diagnostics.valid && diagnostics.solids===12,'Disconnected bodies lost');
    if(expectedThreaded && info.workers>1) assert(diagnostics.parallelTasks>0,'No pthread fitting tasks');
    else assert(diagnostics.parallelTasks===0,'Sequential configuration ran background tasks');
    rejected=false;
    try { module.convert([],{tolerance:0}); } catch { rejected=true; }
    assert(rejected,'Invalid options accepted');
    assert(module.releaseModel(model.modelId) && !module.releaseModel(model.modelId),'Release must be idempotent');
    rejected=false;
    try { module.exportStep(model.modelId); } catch { rejected=true; }
    assert(rejected,'Released model accepted');
    rejected=false;
    try { module.getShapeProperties(selectedEdge); } catch { rejected=true; }
    assert(rejected,'Stale selection accepted');
    assert(retainedStep.every((v,i)=>v===file.data[i]),'Releasing a model invalidated owned bytes');
    return info;
  } finally { for(const id of owned) module.releaseModel(id); }
}
export async function exercisePicking(module, loader, THREE) {
  const result=module.convert([{name:'picking.stl',data:cubeStl()}]);
  assert(result.errors.length===0 && result.models.length===1,'Picking model conversion');
  const model=result.models[0];
  try {
    const glb=module.exportGltf(model.modelId);
    const loaded=await loader.parseAsync(glb.data.buffer,'');
    loaded.scene.updateMatrixWorld(true);
    const meshes=[];
    loaded.scene.traverse(object=>{if(object.isMesh) meshes.push(object);});
    assert(meshes.length===6,'GLTFLoader did not preserve selectable faces');
    const ray=new THREE.Raycaster(new THREE.Vector3(.0005,.0005,.01),new THREE.Vector3(0,0,-1));
    const hits=ray.intersectObjects(meshes,false);
    assert(hits.length>0,'Three.js raycast missed the exported cube');
    const selection=hits[0].object.userData;
    assert(selection.modelId===model.modelId && selection.brepgraph_id.kind===2,'Picking lost node extras');
    const properties=module.getShapeProperties(selection);
    assert(properties.properties.some(p=>p.name==='shape_type'&&p.value[0]===4),'Picked face did not resolve to OCCT face');
    assert(module.getClosestPoints(selection,selection).distance<1e-8,'Picked selection exact query');
    const bounds=new THREE.Box3().setFromObject(loaded.scene);
    assert(Math.abs(bounds.max.x-.001)<1e-8 && Math.abs(bounds.max.y-.001)<1e-8 &&
      Math.abs(bounds.min.z+.001)<1e-8,'GLB meter/Y-up bounds');
  } finally { module.releaseModel(model.modelId); }
}
