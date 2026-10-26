import http from 'node:http';
import { readFile } from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const root = path.resolve(fileURLToPath(new URL('../', import.meta.url)));
export function createServer() {
  return http.createServer(async (request, response) => {
    try {
      const url = new URL(request.url, 'http://localhost');
      const isolated = url.pathname.startsWith('/isolated/');
      const name = decodeURIComponent(url.pathname.replace(/^\/(isolated|ordinary)\//, '/'));
      const target = path.resolve(root, '.' + name);
      if (!target.startsWith(root + path.sep) && target !== root) { response.writeHead(403).end(); return; }
      const types = { '.html': 'text/html', '.mjs': 'text/javascript', '.js': 'text/javascript', '.wasm': 'application/wasm' };
      response.setHeader('Content-Type', types[path.extname(target)] || 'application/octet-stream');
      if (isolated) {
        response.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
        response.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
      }
      response.end(await readFile(target));
    } catch { response.writeHead(404).end(); }
  });
}
