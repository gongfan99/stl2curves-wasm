import { chromium } from 'playwright';
import { createServer } from '../wasm/browser_server.mjs';
import { cubeStl, assert } from '../wasm/test_helpers.mjs';
import { mkdir } from 'node:fs/promises';
import path from 'node:path';

// The default exercises the rebuilt package end to end. For UI development
// before rebuilding, this explicit mode substitutes only STEP import with a
// cube conversion; GLB export, picking metadata, and measurements still use WASM.
const importDouble = process.env.STEP_VIEW_IMPORT_DOUBLE === '1';
const server = createServer();
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
const browser = await chromium.launch({ headless: true });
await mkdir('.cache/view-tests', { recursive: true });
const cube = Array.from(cubeStl());

async function hover(page, kind, exclude) {
  for (const candidate of (await page.evaluate(() => window.__stepState())).candidates) {
    if (candidate.kind !== kind || candidate.index === exclude || !candidate.onScreen) continue;
    await page.mouse.move(candidate.x, candidate.y);
    const current = (await page.evaluate(() => window.__stepState())).hover;
    if (current?.kind === kind && current.index === candidate.index) return candidate;
  }
  throw new Error(`No visible ${kind} could be picked`);
}
async function opened(page, name) {
  await page.waitForFunction(name => !document.querySelector('#open-file').disabled &&
    (document.querySelector('#output-name').textContent === name || document.querySelector('#status').dataset.state === 'error'), name, { timeout: 120000 });
  assert(await page.locator('#output-name').textContent() === name, await page.locator('#status-text').textContent());
}
try {
  for (const mode of ['ordinary', 'isolated']) {
    const page = await browser.newPage({ viewport: { width: 1400, height: 1000 }, serviceWorkers: 'block' });
    await page.addInitScript(() => { window.coi = { shouldRegister: () => false }; });
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    await page.route('https://cdn.jsdelivr.net/npm/three@0.186.1/**', async route => {
      const asset = new URL(route.request().url()).pathname.split('/three@0.186.1/')[1];
      await route.fulfill({ path: path.resolve('node_modules/three', asset), contentType: 'text/javascript' });
    });
    if (importDouble) await page.route('**/loader.mjs', async route => {
      const url = route.request().url();
      await route.fulfill({ contentType: 'text/javascript', body: `
        import m from ${JSON.stringify(url + '?actual=1')};
        m.importStep = input => {
          if (!new TextDecoder().decode(input.data).startsWith('ISO-10303-21')) throw new Error('Malformed STEP test input');
          return m.convert([{name:'cube.stl',data:new Uint8Array(${JSON.stringify(cube)})}]).models[0];
        };
        export default m;
      ` });
    });
    await page.route('**/app/view.html', async route => {
      const response = await route.fetch();
      const body = (await response.text()).replace('function render() {', `function render() {
        window.__stepState = () => {
          const view = views[0], rect = view.renderer.domElement.getBoundingClientRect();
          return { modelId, selected: inspection.selected.map(item => item.kind),
            hover: inspection.hover && { kind: selectionKind(inspection.hover.object), index: inspection.hover.object.userData.brepgraph_id.index },
            candidates: inspection.objects.map(object => {
              const positions = object.geometry.attributes.position, kind = selectionKind(object);
              let p;
              if (kind === 'edge') {
                const middle = Math.floor((positions.count - 1) / 2);
                p = new THREE.Vector3().fromBufferAttribute(positions, middle).lerp(new THREE.Vector3().fromBufferAttribute(positions, middle + 1), .5);
              } else p = new THREE.Box3().setFromBufferAttribute(positions).getCenter(new THREE.Vector3());
              object.localToWorld(p); p.project(sharedCamera);
              return { kind, index: object.userData.brepgraph_id.index, onScreen: Math.abs(p.x)<.95 && Math.abs(p.y)<.95,
                x: rect.left+(p.x+1)*rect.width/2, y:rect.top+(1-p.y)*rect.height/2 };
            }) };
        };`);
      await route.fulfill({ response, body });
    });
    const base = `http://127.0.0.1:${server.address().port}/${mode}`;
    await page.goto(`${base}/app/view.html`);
    await page.waitForFunction(() => document.querySelectorAll('canvas.render').length === 1);
    assert(await page.evaluate(() => getComputedStyle(document.body).margin === '0px'), 'Shared workspace stylesheet must load');
    assert(await page.locator('#output-empty').isVisible(), 'Empty drop area');
    assert(await page.locator('input[type=file]').count() === 1, 'One upload control');
    await page.screenshot({ path: `.cache/view-tests/${mode}-empty.png` });
    const step = await page.evaluate(async ({ cube, base }) => {
      // Fixture generation uses the single-thread factory even on the isolated
      // page; threaded model operations must run in the application Worker.
      const url = `${base}/app/lib/stl2curves-wasm/v0.1.0/single_thread/stl2curves.js`;
      const factory = (await import(url)).default;
      const m = await factory({ locateFile: name => new URL(name, url).href });
      try {
        const model = m.convert([{ name: 'cube.stl', data: new Uint8Array(cube) }]).models[0];
        return Array.from(m.exportStep(model.modelId).data);
      } finally { m.dispose(); }
    }, { cube, base });
    await page.locator('#file').setInputFiles({ name: 'cube.STP', mimeType: 'application/step', buffer: Buffer.from(step) });
    await opened(page, 'cube.STP');
    assert((await page.locator('#output-stats').textContent()).includes('6 faces'), 'STEP face count');
    assert((await page.locator('#runtime').textContent()).startsWith(mode === 'isolated' ? 'Threaded' : 'Single-thread'), 'Same WASM variant selection as index');
    const face = await hover(page, 'face');
    await page.mouse.click(face.x, face.y);
    const otherFace = await hover(page, 'face', face.index);
    await page.mouse.click(otherFace.x, otherFace.y);
    await page.waitForFunction(() => document.querySelector('#measurement-text').textContent.includes('face to face:'));
    await page.locator('#clear-measurement').click();
    const edge = await hover(page, 'edge');
    await page.waitForFunction(() => document.querySelector('#shape-tooltip').textContent === 'Length: 1.00 mm' && !document.querySelector('#shape-tooltip').hidden);
    await page.mouse.click(edge.x, edge.y);
    const vertex = await hover(page, 'vertex');
    await page.mouse.click(vertex.x, vertex.y);
    await page.waitForFunction(() => document.querySelector('#measurement-text').textContent.includes('edge to vertex:'));
    assert(await page.locator('#distance-label').isVisible(), 'Measurement annotation');
    await page.screenshot({ path: `.cache/view-tests/${mode}-measurement.png` });
    await page.keyboard.press('Escape');
    const canvas = await page.locator('canvas.render').boundingBox();
    await page.mouse.move(canvas.x+canvas.width/2, canvas.y+canvas.height/2);
    await page.mouse.down(); await page.mouse.move(canvas.x+canvas.width*.7, canvas.y+canvas.height*.6, {steps: 8}); await page.mouse.up();
    assert((await page.evaluate(() => window.__stepState())).selected.length === 0, 'Orbit drag must not select');
    await page.locator('.reset-view').click();
    const firstId = (await page.evaluate(() => window.__stepState())).modelId;
    await page.locator('#output-view').evaluate((element, bytes) => {
      const dataTransfer = new DataTransfer(); dataTransfer.items.add(new File([new Uint8Array(bytes)], 'dropped.step'));
      element.dispatchEvent(new DragEvent('dragenter', { bubbles: true, dataTransfer }));
      element.dispatchEvent(new DragEvent('drop', { bubbles: true, dataTransfer }));
    }, step);
    await opened(page, 'dropped.step');
    assert((await page.evaluate(() => window.__stepState())).modelId !== firstId, 'Replacing file creates a new live model');
    await page.locator('#file').setInputFiles({ name: 'invalid.step', mimeType: 'application/step', buffer: Buffer.from('bad STEP') });
    await page.waitForFunction(() => document.querySelector('#status').dataset.state === 'error' && !document.querySelector('#open-file').disabled);
    assert(await page.locator('#output-name').textContent() === 'dropped.step', 'Failed import preserves displayed model');
    await page.locator('#file').setInputFiles({ name: 'wrong.stl', mimeType: 'application/octet-stream', buffer: Buffer.from('bad') });
    await page.waitForFunction(() => document.querySelector('#status-text').textContent.includes('Choose one .step'));
    await page.setViewportSize({ width: 390, height: 844 });
    await page.locator('.reset-view').click();
    assert(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), 'Mobile layout must not overflow');
    await page.screenshot({ path: `.cache/view-tests/${mode}-mobile.png` });
    assert(errors.length === 0, errors.join('\n'));
    await page.close();
    console.log(`${mode}: STEP viewer passed${importDouble ? ' (STEP import test double; native import tested separately)' : ''}`);
  }
} finally {
  await browser.close();
  await new Promise(resolve => server.close(resolve));
}
