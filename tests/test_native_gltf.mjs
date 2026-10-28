// Three.js interoperability check for a native-generated GLB; does not load WASM.
import { readFile, writeFile } from 'node:fs/promises';
import { GLTFLoader } from 'three/examples/jsm/loaders/GLTFLoader.js';
import { Box3, Raycaster, Vector3 } from 'three';
import assert from 'node:assert/strict';

const [input, output] = process.argv.slice(2);
if (!input || !output) throw new Error('Usage: node tests/test_native_gltf.mjs INPUT.glb PICKED.json');
const file = await readFile(input);
const data = Uint8Array.from(file);
const { scene } = await new GLTFLoader().parseAsync(data.buffer, '');
scene.updateMatrixWorld(true);
const faces = [], edges = [], vertices = [];
scene.traverse(object => {
  if (object.isMesh) faces.push(object);
  if (object.isLine) edges.push(object);
  if (object.isPoints) vertices.push(object);
});
assert.equal(faces.length, 6);
assert.equal(edges.length, 12);
assert.equal(vertices.length, 8);
for (const face of faces) {
  const geometry = face.geometry;
  assert.ok(geometry.index.array instanceof Uint32Array, 'Use 32-bit triangle indices');
  const a = new Vector3(), b = new Vector3(), c = new Vector3(), normal = new Vector3();
  geometry.computeBoundingBox();
  const outward = geometry.boundingBox.getCenter(new Vector3()).sub(new Vector3(.5,.5,.5)).normalize();
  normal.fromBufferAttribute(geometry.attributes.normal, 0);
  assert.ok(outward.dot(normal) > .99, 'Cube face normal must point outwards');
  for (let i = 0; i < geometry.index.count; i += 3) {
    const first = geometry.index.getX(i);
    a.fromBufferAttribute(geometry.attributes.position, first);
    b.fromBufferAttribute(geometry.attributes.position, geometry.index.getX(i + 1));
    c.fromBufferAttribute(geometry.attributes.position, geometry.index.getX(i + 2));
    normal.fromBufferAttribute(geometry.attributes.normal, first);
    const winding = b.sub(a).cross(c.sub(a)).normalize();
    assert.ok(winding.dot(normal) > .99, 'Triangle winding and face normal disagree');
  }
}
const bounds = new Box3().setFromObject(scene);
assert.ok(Math.abs(bounds.max.x - .001) < 1e-8 && Math.abs(bounds.max.y - .001) < 1e-8 &&
  Math.abs(bounds.min.z + .001) < 1e-8, 'Meter/Y-up scene bounds');
const ray = new Raycaster(new Vector3(.0005, .0005, .01), new Vector3(0, 0, -1));
const hits = ray.intersectObjects(faces, false);
assert.ok(hits.length > 0, 'Raycast must hit an exported face');
assert.ok(Math.abs(hits[0].distance - .01) < 1e-8, 'Raycast should hit the front surface');
const selection = hits[0].object.userData;
assert.equal(selection.brepgraph_id.kind, 2);
assert.ok(selection.modelId > 0 && selection.usagePath.length > 0);
await writeFile(output, JSON.stringify(selection));
console.log('Native GLB loads in Three.js with face/edge/vertex IDs and valid raycast selection');
