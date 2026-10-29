# python3 src/build.py → standalone.html: index.html with src/core.js and HiGHS (JS + WASM) embedded, for offline use (needs `npm install`)
import base64, pathlib
root = pathlib.Path(__file__).resolve().parent.parent
highs = root / 'node_modules' / 'highs' / 'build'
page = (root / 'index.html').read_text(encoding='utf-8')
hj, core = (highs / 'highs.js').read_text(encoding='utf-8'), (root / 'src' / 'core.js').read_text(encoding='utf-8')
wasm = base64.b64encode((highs / 'highs.wasm').read_bytes()).decode()
for s in (hj, core): assert '</script' not in s.lower() and '<!--' not in s
web = ['importScripts(self.HIGHS_CDN + "highs.js");', 'importScripts(self.GSTO_CORE_URL);', '<script type="text/plain" id="wasmB64"></script>']
assert all(page.count(w) == 1 for w in web), 'index.html markers not found'
out = page.replace(web[0], hj, 1).replace(web[1], core, 1).replace(web[2], f'<script type="text/plain" id="wasmB64">{wasm}</script>', 1)
(root / 'standalone.html').write_text(out, encoding='utf-8')
print(f'standalone.html {len(out)/1e6:.2f} MB')
