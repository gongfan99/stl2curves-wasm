// Exercise the actual single-page app against the local packaged WASM assets.
import { chromium } from 'playwright';
import { createServer } from '../wasm/browser_server.mjs';
import { cubeStl, inchCube3mf, assert } from '../wasm/test_helpers.mjs';
import { readFile, mkdir } from 'node:fs/promises';

const server = createServer();
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
const browser = await chromium.launch({ headless: true });
await mkdir('.cache/app-tests', { recursive: true });
const bytes = cubeStl();
async function assertPerspective(page, label) {
  const state = await page.evaluate(() => window.__viewerState());
  assert(state.sameCamera && state.sameTarget, `${label}: viewers must share their camera and orbit target`);
  for (const key of ['projection', 'world', 'target', 'renderSize']) {
    assert(JSON.stringify(state.views[0][key]) === JSON.stringify(state.views[1][key]), `${label}: different ${key}`);
  }
  const picking = await page.evaluate(() => window.__inspectionState());
  assert(picking.thresholds[0] === picking.dimension * .006 && picking.thresholds[1] === picking.dimension * .008,
    `${label}: selection sensitivity must follow model size`);
}
async function hoverCandidate(page, kind, exclude) {
  const candidates = await page.evaluate(kind => window.__inspectionState().candidates.filter(c => c.kind === kind), kind);
  for (const candidate of candidates) {
    if (candidate.index === exclude || !candidate.onScreen) continue;
    await page.mouse.move(candidate.x, candidate.y);
    const hover = await page.evaluate(() => window.__inspectionState().hover);
    if (hover?.kind === kind && hover.index === candidate.index) return candidate;
  }
  throw new Error(`No visible ${kind} could be picked`);
}
async function assertNearbyPick(page, candidate) {
  await page.mouse.move(candidate.near.x, candidate.near.y);
  const picked = await page.evaluate(() => window.__inspectionState().hover);
  assert(picked?.kind === candidate.kind && picked.index === candidate.index,
    `Near ${candidate.kind}: world-space pick tolerance was lost to parent scaling`);
  await page.mouse.move(candidate.x, candidate.y);
}
async function exerciseInspection(page, mode, expectedSize = 1) {
  assert(await page.locator('#toggle-measurement, #selection-options').count() === 0, 'Inspection must not require toggles');
  await page.locator('.reset-view').first().click();
  const face1 = await hoverCandidate(page, 'face');
  assert(await page.locator('#shape-tooltip').isHidden(), 'Faces should highlight without a tooltip');
  await page.mouse.click(face1.x, face1.y);
  const face2 = await hoverCandidate(page, 'face', face1.index);
  await page.mouse.click(face2.x, face2.y);
  await page.waitForFunction(() => document.querySelector('#measurement-text').textContent.includes('face to face:'));
  assert(await page.locator('#distance-label').textContent() === '0 mm', 'Adjacent cube faces should have zero minimum distance');
  await page.locator('#clear-measurement').click();
  const edge1 = await hoverCandidate(page, 'edge');
  await assertNearbyPick(page, edge1);
  await page.waitForFunction(expected => !document.querySelector('#shape-tooltip').hidden &&
    document.querySelector('#shape-tooltip').textContent === `Length: ${expected.toFixed(2)} mm`, expectedSize);
  assert((await page.locator('#shape-tooltip').textContent()).includes(`${expectedSize.toFixed(2)} mm`), 'Exact edge length in model millimeters');
  await page.mouse.click(edge1.near.x, edge1.near.y);
  const edge2 = await hoverCandidate(page, 'edge', edge1.index);
  await page.mouse.click(edge2.x, edge2.y);
  await page.waitForFunction(() => document.querySelector('#measurement-text').textContent.includes('edge to edge:'));
  const edgeDistance = Number((await page.locator('#distance-label').textContent()).replace(/[^0-9.]/g, ''));
  assert([0, expectedSize, Math.SQRT2 * expectedSize].some(d => Math.abs(d - edgeDistance) < .0001), 'Unexpected cube edge distance');
  await page.locator('#clear-measurement').click();
  const first = await hoverCandidate(page, 'vertex');
  await assertNearbyPick(page, first);
  await page.mouse.click(first.near.x, first.near.y);
  assert((await page.locator('#measurement-text').textContent()).includes('Select the second'), 'First selection');
  await page.mouse.move(first.x, first.y);
  await page.mouse.down();
  await page.mouse.move(first.x + 30, first.y + 10, { steps: 4 });
  await page.mouse.move(first.x, first.y, { steps: 4 });
  await page.mouse.up();
  assert((await page.evaluate(() => window.__inspectionState().selected)).length === 1, 'Orbit drag must not select geometry');
  const second = await hoverCandidate(page, 'vertex', first.index);
  await page.mouse.click(second.near.x, second.near.y);
  await page.waitForFunction(() => document.querySelector('#measurement-text').textContent.includes('vertex to vertex:'));
  const state = await page.evaluate(() => window.__inspectionState());
  const expected = Math.hypot(...first.point.map((value, i) => (value - second.point[i]) / state.mmScale));
  const measured = Number((await page.locator('#distance-label').textContent()).replace(/[^0-9.]/g, ''));
  assert(Math.abs(measured - expected) < .0001, `Exact distance ${measured}, expected ${expected}`);
  assert(state.line.length === 6 && state.selected.length === 2, 'Missing measurement graphics');
  for (const [i, candidate] of [first, second].entries()) {
    assert(candidate.point.every((v, axis) => Math.abs(v - state.line[i * 3 + axis]) < 1e-5), 'Measurement endpoint transform');
  }
  assert((await page.locator('#measurement-text').textContent()).includes('ΔZ:'), 'Missing axis deltas');
  await assertPerspective(page, 'measurement');
  await page.screenshot({ path: `.cache/app-tests/${mode}-measurement-${expectedSize}.png` });
  await page.mouse.click(first.x, first.y);
  assert((await page.evaluate(() => window.__inspectionState().selected)).length === 1, 'Third click should restart measurement');
  await page.locator('#clear-measurement').click();
  assert((await page.evaluate(() => window.__inspectionState().selected)).length === 0, 'Clear should remove selections');
  await page.keyboard.press('Escape');
  assert(await page.locator('#measurement-readout').isVisible(), 'Measurement must stay active');
}
try {
  for (const mode of ['ordinary', 'isolated']) {
    const page = await browser.newPage({ viewport: { width: 1500, height: 1050 }, acceptDownloads: true });
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    // Inspect the actual render state without adding a debug API to the app.
    await page.route('**/app/index.html', async route => {
      const response = await route.fetch();
      const body = (await response.text()).replace('function render() {', `function render() {
        window.__viewerState = () => ({
          sameCamera: views.every(v => v.camera === views[0].camera),
          sameTarget: views.every(v => v.controls.target === views[0].controls.target),
          views: views.map(v => ({
            projection: v.camera.projectionMatrix.toArray(), world: v.camera.matrixWorld.toArray(),
            target: v.controls.target.toArray(), renderSize: v.renderer.getViewport(new THREE.Vector4()).toArray().slice(2),
          })),
        });
        window.__inspectionState = () => {
          const view = views[1];
          const rect = view.renderer.domElement.getBoundingClientRect(), vp = view.renderer.getViewport(new THREE.Vector4());
          return {
            mmScale: normalization?.scale * .001,
            dimension: inspection.dimension, thresholds: [raycaster.params.Line.threshold, raycaster.params.Points.threshold],
            selected: inspection.selected.map(s => s.shape),
            hover: inspection.hover && {kind: selectionKind(inspection.hover.object), index: inspection.hover.object.userData.brepgraph_id.index},
            line: inspection.distanceGroup ? Array.from(inspection.distanceGroup.children[0].geometry.attributes.position.array) : [],
            candidates: inspection.objects.map(object => {
              const position = object.geometry.attributes.position;
              const p = object.isPoints ? new THREE.Vector3().fromBufferAttribute(position, 0) : object.isLine ? new THREE.Vector3().fromBufferAttribute(position, Math.floor((position.count - 1) / 2)).lerp(new THREE.Vector3().fromBufferAttribute(position, Math.floor((position.count - 1) / 2) + 1), .5) : new THREE.Box3().setFromBufferAttribute(position).getCenter(new THREE.Vector3());
              object.localToWorld(p);
              const projected = p.clone().project(sharedCamera);
              const kind = selectionKind(object);
              let offset = new THREE.Vector3(1,0,0).applyQuaternion(sharedCamera.quaternion);
              if (kind === 'edge') {
                const k = Math.floor((position.count - 1) / 2);
                const a = object.localToWorld(new THREE.Vector3().fromBufferAttribute(position, k));
                const b = object.localToWorld(new THREE.Vector3().fromBufferAttribute(position, k + 1));
                offset = new THREE.Vector3().crossVectors(p.clone().sub(sharedCamera.position), b.sub(a)).normalize();
              }
              const near = p.clone().addScaledVector(offset, inspection.dimension * (kind === 'vertex' ? .008 : .006) * .75).project(sharedCamera);
              return {kind:selectionKind(object), index:object.userData.brepgraph_id.index, point:p.toArray(),
                near: {x:rect.left+(vp.x+(near.x+1)*vp.z/2)*rect.width/view.element.clientWidth,
                  y:rect.bottom-(vp.y+(near.y+1)*vp.w/2)*rect.height/view.element.clientHeight},
                onScreen: Math.abs(projected.x)<.95 && Math.abs(projected.y)<.95,
                x:rect.left+(vp.x+(projected.x+1)*vp.z/2)*rect.width/view.element.clientWidth,
                y:rect.bottom-(vp.y+(projected.y+1)*vp.w/2)*rect.height/view.element.clientHeight};
            }),
          };
        };`);
      await route.fulfill({ response, body });
    });
    await page.goto(`http://127.0.0.1:${server.address().port}/${mode}/app/index.html`);
    await page.waitForFunction(() => document.querySelectorAll('canvas.render').length === 2);
    assert(!await page.locator('#fuse').isChecked(), 'Fuse should default to false');
    await page.waitForFunction(() => document.querySelector('#input-name').textContent === 'sample-cube.stl');
    assert(!await page.locator('#convert').isDisabled(), 'Startup sample should be ready to convert');
    assert(await page.locator('#input-empty').isHidden(), 'Startup sample should be visible');
    await page.screenshot({ path: `.cache/app-tests/${mode}-sample.png` });
    await page.locator('#convert').click();
    await page.waitForFunction(() => !document.querySelector('#download').disabled, undefined, { timeout: 120000 });
    assert((await page.locator('#output-stats').textContent()).includes('6 faces'), 'Startup cube conversion');
    await page.locator('#file').setInputFiles({ name: 'cube.stl', mimeType: 'application/octet-stream', buffer: Buffer.from(bytes) });
    await page.waitForFunction(() => document.querySelector('#status-text').textContent.includes('loaded.'));
    assert(await page.locator('#input-stats').textContent() === '12 triangles', 'Source triangle count');
    await page.locator('#convert').click();
    await page.waitForFunction(() => document.querySelector('#download').disabled === false, undefined, { timeout: 120000 });
    assert((await page.locator('#output-stats').textContent()).includes('6 faces'), 'Converted cube face count');
    await assertPerspective(page, 'conversion');
    const runtime = await page.locator('#runtime').textContent();
    assert(runtime.startsWith(mode === 'isolated' ? 'Threaded' : 'Single-thread'), 'Wrong app WASM variant');
    await exerciseInspection(page, mode);
    const axisBefore = await page.locator('#input-view .axis').evaluate(c => c.toDataURL());
    for (const id of ['input-view', 'output-view']) {
      const box = await page.locator(`#${id} canvas.render`).boundingBox();
      await page.mouse.move(box.x + box.width * .5, box.y + box.height * .5);
      await page.mouse.down();
      await page.mouse.move(box.x + box.width * .7, box.y + box.height * .6, { steps: 10 });
      await page.mouse.up();
      const axes = await page.locator('.axis').evaluateAll(nodes => nodes.map(c => c.toDataURL()));
      assert(axes[0] === axes[1] && axes[0] !== axisBefore, `${id}: cameras did not synchronize`);
      await assertPerspective(page, `${id} orbit`);
      await page.mouse.wheel(0, -60);
      await assertPerspective(page, `${id} zoom`);
      await page.mouse.down({ button: 'right' });
      await page.mouse.move(box.x + box.width * .65, box.y + box.height * .55, { steps: 5 });
      await page.mouse.up({ button: 'right' });
      await assertPerspective(page, `${id} pan`);
    }
    await page.screenshot({ path: `.cache/app-tests/${mode}-converted.png` });
    const downloaded = {};
    assert(await page.locator('[data-format]').count() === 2, 'Only GLB and STEP downloads should be offered');
    for (const format of ['glb', 'step']) {
      await page.locator('#download').click();
      const downloadPromise = page.waitForEvent('download', { timeout: 120000 });
      await page.locator(`[data-format="${format}"]`).click();
      const download = await downloadPromise;
      assert(download.suggestedFilename() === `cube.${format}`, 'Download filename');
      downloaded[format] = await readFile(await download.path());
    }
    const gltf = JSON.parse(downloaded.glb.subarray(20, 20 + downloaded.glb.readUInt32LE(12)).toString('utf8'));
    assert(gltf.nodes.some(n => n.extras?.brepgraph_id?.kind === 2), 'glTF selection metadata missing');
    const loadedFaces = await page.evaluate(async bytes => {
      const { GLTFLoader } = await import('three/addons/loaders/GLTFLoader.js');
      const loaded = await new GLTFLoader().parseAsync(new Uint8Array(bytes).buffer, '');
      let faces = 0;
      loaded.scene.traverse(object => { if (object.isMesh) faces++; });
      return faces;
    }, Array.from(downloaded.glb));
    assert(loadedFaces === 6, 'Downloaded GLB cannot be reopened');
    assert(downloaded.glb.readUInt32LE(0) === 0x46546c67, 'Invalid GLB');
    assert(downloaded.step.toString('utf8').includes('ISO-10303-21'), 'Invalid STEP');
    // Repeat conversion with explicit simplification and fusion options.
    await page.locator('#settings').click();
    assert(await page.locator('#settings-popup').isVisible(), 'Settings popup should open');
    await page.locator('#simplify').selectOption('custom');
    await page.locator('#custom-simplify').fill('0.001');
    await page.locator('#accuracy').fill('0.05');
    await page.locator('#fuse').check();
    await page.locator('#settings-done').click();
    assert(await page.locator('#settings-popup').isHidden(), 'Done should close settings');
    await page.locator('#convert').click();
    await page.waitForFunction(() => !document.querySelector('#download').disabled, undefined, { timeout: 120000 });
    assert((await page.locator('#diagnostics-text').textContent()).includes('Fuse: on'), 'Custom options not applied');
    // Use a real browser drop event to replace the previous retained model.
    await page.evaluate(data => {
      const transfer = new DataTransfer();
      transfer.items.add(new File([new Uint8Array(data)], 'inch-cube.3mf', { type: 'model/3mf' }));
      document.querySelector('#input-view').dispatchEvent(new DragEvent('drop', { bubbles: true, dataTransfer: transfer }));
    }, Array.from(inchCube3mf()));
    await page.waitForFunction(() => document.querySelector('#input-name').textContent === 'inch-cube.3mf');
    await page.locator('#convert').click();
    await page.waitForFunction(() => !document.querySelector('#download').disabled, undefined, { timeout: 120000 });
    assert((await page.locator('#output-stats').textContent()).includes('6 faces'), '3MF conversion');
    await exerciseInspection(page, mode, 25.4);
    await page.screenshot({ path: `.cache/app-tests/${mode}-3mf.png` });
    // Malformed input should preserve the previous valid model/downloads.
    await page.locator('#file').setInputFiles({ name: 'broken.3mf', mimeType: 'model/3mf', buffer: Buffer.from('bad zip') });
    await page.waitForFunction(() => document.querySelector('#status').dataset.state === 'error');
    assert(!await page.locator('#download').isDisabled(), 'Invalid replacement destroyed previous model');
    if (process.env.STL2CURVES_APP_FIXTURE) {
      await page.locator('#file').setInputFiles(process.env.STL2CURVES_APP_FIXTURE);
      await page.waitForFunction(() => document.querySelector('#status-text').textContent.includes('loaded.'));
      await page.locator('#convert').click();
      await page.waitForFunction(() => !document.querySelector('#download').disabled, undefined, { timeout: 120000 });
      assert(await page.locator('#output-badge').textContent() === 'Valid CAD', 'Fixture did not produce valid CAD');
      await assertPerspective(page, 'fixture');
      if (process.env.STL2CURVES_APP_REQUIRE_CIRCLE) {
        const candidates = await page.evaluate(() => window.__inspectionState().candidates.filter(c => c.kind === 'edge' && c.onScreen));
        let found = false;
        for (const candidate of candidates) {
          await page.mouse.move(candidate.x, candidate.y);
          try {
            await page.waitForFunction(() => !document.querySelector('#shape-tooltip').hidden && document.querySelector('#shape-tooltip').textContent.includes('Radius:'), undefined, { timeout: 300 });
            found = true; break;
          } catch { /* Try the next visible edge. */ }
        }
        assert(found, 'No exact circle radius tooltip found on fixture');
        assert(/^Radius: \d+\.\d{2} mm$/.test(await page.locator('#shape-tooltip').textContent()), 'Circle tooltip must contain only its radius');
        await page.screenshot({ path: `.cache/app-tests/${mode}-circle-hover.png` });
      }
      await page.screenshot({ path: `.cache/app-tests/${mode}-fixture.png` });
    }
    for (const width of [1501, 1180, 970, 760, 390]) {
      await page.setViewportSize({ width, height: 844 });
      await page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));
      await assertPerspective(page, `resize ${width}`);
      await page.locator('.reset-view').first().click();
      await assertPerspective(page, `fit ${width}`);
    }
    await page.setViewportSize({ width: 390, height: 844 });
    await page.screenshot({ path: `.cache/app-tests/${mode}-mobile.png`, fullPage: true });
    assert(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), 'Mobile horizontal overflow');
    assert(errors.length === 0, `Browser errors: ${errors.join('; ')}`);
    console.log(`${mode}: file/drop, STL/3MF, synchronized cameras, hover, exact measurements, options, repeat conversion, GLB/STEP downloads and responsive layout passed (${runtime})`);
    await page.close();
  }
} finally { await browser.close(); await new Promise(resolve => server.close(resolve)); }
