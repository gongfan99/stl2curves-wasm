import { readFile, readdir } from 'node:fs/promises';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const root = new URL('../', import.meta.url);
for (const name of ['package.json', 'vcpkg.json', 'vcpkg-configuration.json', 'tests/regression_cases.json', 'vcpkg-overlays/opencascade/vcpkg.json']) {
  JSON.parse(await readFile(new URL(name, root), 'utf8'));
}
for (const directory of ['wasm', 'scripts']) {
  for (const name of await readdir(new URL(`${directory}/`, root))) {
    if (!/\.(mjs|js)$/.test(name)) continue;
    const file = new URL(`${directory}/${name}`, root);
    const result = spawnSync(process.execPath, ['--check', fileURLToPath(file)], { encoding: 'utf8' });
    if (result.status !== 0) throw new Error(result.stderr);
  }
}
const workflow = await readFile(new URL('.github/workflows/build_wasm.yml', root), 'utf8');
if (!/workflow_dispatch:/.test(workflow) || /^\s+(push|pull_request):/m.test(workflow)) {
  throw new Error('WASM compilation must remain manually triggered');
}
console.log('JSON and JavaScript syntax checks passed; no WASM was compiled or executed.');
