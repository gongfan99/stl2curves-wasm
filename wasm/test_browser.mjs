// Do not run until WASM artifacts have been built and testing is requested.
import { chromium } from 'playwright';
import { createServer } from './browser_server.mjs';
import { assert } from './test_helpers.mjs';

const server = createServer();
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
let browser;
try {
  browser = await chromium.launch({ headless: true });
  for (const mode of ['ordinary', 'isolated']) {
    const page = await browser.newPage();
    await page.goto(`http://127.0.0.1:${server.address().port}/${mode}/wasm/test_browser.html`);
    await page.waitForFunction(() => window.__testResult !== undefined, undefined, { timeout: 120000 });
    const result = await page.evaluate(() => window.__testResult);
    assert(result.ok, `${mode}: ${result.error}`);
    assert(result.info.threaded === (mode === 'isolated'), 'Loader selected incorrect variant');
    await page.close();
  }
  console.log('Browser blob-loader tests passed for both variants');
} finally {
  await browser?.close();
  await new Promise(resolve => server.close(resolve));
}
