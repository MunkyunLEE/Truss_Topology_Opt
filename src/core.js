// GSTO core — 3D ground structure truss layout optimization (pin-jointed, one-cell connectivity, no X-crossing)
const GSTO = (() => {
  const OFFS = [                                   // 13 half-neighbour offsets (one-cell connectivity)
    [1,0,0],[0,1,0],[0,0,1],                       // edges
    [1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1], // face diagonals
    [1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1]];          // body diagonals
  const OFF_KIND = [0,0,0,1,1,1,1,1,1,2,2,2,2];   // 0 edge, 1 face diag, 2 body diag

  const sub = (a,b) => [a[0]-b[0],a[1]-b[1],a[2]-b[2]];
  const dot = (a,b) => a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  const amaxOf = arr => { let m = -Infinity; for (const v of arr) if (v > m) m = v; return m; };   // max without argument spreading (large arrays)
  const cross = (a,b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
  const norm = a => Math.hypot(a[0],a[1],a[2]);

  function inBoxOpen(p, bx, tol) {                 // strictly interior of AABB
    for (let a = 0; a < 3; a++) if (!(p[a] > bx.min[a]+tol && p[a] < bx.max[a]-tol)) return false;
    return true;
  }
  function inBoxClosed(p, bx, tol) {               // inside or on boundary of AABB
    for (let a = 0; a < 3; a++) if (p[a] < bx.min[a]-tol || p[a] > bx.max[a]+tol) return false;
    return true;
  }
  function segHitsBoxInterior(p, q, bx, tol) {     // slab test against box shrunk by tol
    let t0 = 0, t1 = 1;
    for (let a = 0; a < 3; a++) {
      const lo = bx.min[a]+tol, hi = bx.max[a]-tol, d = q[a]-p[a];
      if (Math.abs(d) < 1e-14) { if (p[a] <= lo || p[a] >= hi) return false; continue; }
      let ta = (lo-p[a])/d, tb = (hi-p[a])/d;
      if (ta > tb) [ta, tb] = [tb, ta];
      t0 = Math.max(t0, ta); t1 = Math.min(t1, tb);
      if (t0 >= t1) return false;
    }
    return t1 - t0 > 1e-9;
  }
  function pointInSegInterior(m, a, b, tol) {      // node m strictly inside segment ab
    const d = sub(b,a), L = norm(d), u = [d[0]/L,d[1]/L,d[2]/L], w = sub(m,a), t = dot(w,u);
    return t > tol && t < L-tol && norm(cross(w,u)) < tol;
  }
  function segsCross(a, b, c, d, tol) {            // proper interior crossing or collinear overlap
    const d1 = sub(b,a), d2 = sub(d,c), n = cross(d1,d2), nn = norm(n), L1 = norm(d1), L2 = norm(d2);
    if (nn < 1e-12*L1*L2) {                        // parallel
      if (norm(cross(sub(c,a),d1))/L1 > tol) return false;
      const u = [d1[0]/L1,d1[1]/L1,d1[2]/L1], s0 = dot(sub(c,a),u), s1 = dot(sub(d,a),u);
      return Math.min(Math.max(s0,s1), L1) - Math.max(Math.min(s0,s1), 0) > tol;
    }
    if (Math.abs(dot(sub(c,a), n))/nn > tol) return false;  // not coplanar
    const w = sub(c,a), s = dot(cross(w,d2), n)/(nn*nn), t = dot(cross(w,d1), n)/(nn*nn);
    const e1 = tol/L1, e2 = tol/L2;
    return s > e1 && s < 1-e1 && t > e2 && t < 1-e2;
  }

  // ---------- face interpolation nodes: grid lines piercing region box faces ----------
  function faceNodes(P) {
    const { origin, h, nx, ny, nz } = P.grid, N = [nx, ny, nz], tol = 1e-6*h, out = [];
    const regs = [...(P.supports||[]), ...(P.loads||[]), ...(P.obstacles||[])].filter(r => r.faceNodes);
    const lo = origin, hi = origin.map((o, t) => o + N[t]*h);
    const vals = (t, a, b) => {                    // grid coordinates within [a,b] on axis t, plus the limits themselves
      const v = new Set([a, b]);
      for (let k = Math.max(0, Math.ceil((a-origin[t])/h - 1e-9)); k <= Math.min(N[t], Math.floor((b-origin[t])/h + 1e-9)); k++) v.add(origin[t] + k*h);
      return [...v].filter(x => x >= lo[t]-tol && x <= hi[t]+tol);
    };
    regs.forEach(r => {
      const mn = r.min.map((v, t) => Math.max(Math.min(v, r.max[t]), lo[t])), mx = r.max.map((v, t) => Math.min(Math.max(v, r.min[t]), hi[t]));
      for (let a = 0; a < 3; a++) {
        const [b, c] = [0,1,2].filter(t => t !== a), vb = vals(b, mn[b], mx[b]), vc = vals(c, mn[c], mx[c]);
        for (const fa of (mx[a]-mn[a] > tol ? [mn[a], mx[a]] : [mn[a]])) vb.forEach(xb => vc.forEach(xc => { const p = [0,0,0]; p[a] = fa; p[b] = xb; p[c] = xc; out.push(p); }));
      }
    });
    const seen = new Set();                        // dedupe and drop points already on the grid
    return out.filter(p => {
      const key = p.map(v => Math.round(v/tol)).join(',');
      if (seen.has(key)) return false; seen.add(key);
      return !p.every((v, t) => Math.abs((v-origin[t])/h - Math.round((v-origin[t])/h)) < 1e-9);
    });
  }

  // ---------- ground structure ----------
  function buildGroundStructure(P) {
    const { origin, h, nx, ny, nz } = P.grid;      // origin [mm], spacing h [mm], cell counts
    const obs = P.obstacles || [];                 // [{min:[..],max:[..]}]
    const kinds = P.options?.kinds || [1,1,1];    // allowed member kinds: [edge h, face diag √2h, body diag √3h]
    const tol = 1e-6*h;
    const X = [], ijk = [], gid = new Int32Array((nx+1)*(ny+1)*(nz+1)).fill(-1);
    const gidx = (i,j,k) => (i<0||j<0||k<0||i>nx||j>ny||k>nz) ? -1 : gid[(k*(ny+1)+j)*(nx+1)+i];
    for (let k = 0; k <= nz; k++) for (let j = 0; j <= ny; j++) for (let i = 0; i <= nx; i++) {
      const p = [origin[0]+i*h, origin[1]+j*h, origin[2]+k*h];
      if (obs.some(b => inBoxOpen(p, b, tol))) continue;
      gid[(k*(ny+1)+j)*(nx+1)+i] = X.length; X.push(p); ijk.push([i,j,k]);
    }
    const nGrid = X.length;
    const E = [], kind = [], mid = new Map();      // members [a,b], kind, key->member id
    const blocked = (p,q) => obs.some(b => segHitsBoxInterior(p, q, b, tol));
    const addM = (a, b, kd) => {
      const key = a < b ? a+','+b : b+','+a;
      if (mid.has(key) || blocked(X[a], X[b])) return -1;
      mid.set(key, E.length); E.push(a < b ? [a,b] : [b,a]); kind.push(kd); return E.length-1;
    };
    const mOf = (i,j,k,o) => {                     // grid member id from (i,j,k)+offset o, or -1
      const a = gidx(i,j,k), b = gidx(i+OFFS[o][0], j+OFFS[o][1], k+OFFS[o][2]);
      if (a < 0 || b < 0) return -1;
      const key = a < b ? a+','+b : b+','+a;
      return mid.has(key) ? mid.get(key) : -1;
    };
    for (let k = 0; k <= nz; k++) for (let j = 0; j <= ny; j++) for (let i = 0; i <= nx; i++) {
      const a = gidx(i,j,k); if (a < 0) continue;
      OFFS.forEach((o, oi) => { const b = gidx(i+o[0], j+o[1], k+o[2]); if (b >= 0 && kinds[OFF_KIND[oi]]) addM(a, b, OFF_KIND[oi]); });
    }
    const groups = [];                             // crossing groups: at most one member per group
    const pushG = g => { g = g.filter(m => m >= 0); if (g.length > 1) groups.push(g); };
    for (let k = 0; k <= nz; k++) for (let j = 0; j <= ny; j++) for (let i = 0; i <= nx; i++) {
      if (i < nx && j < ny) pushG([mOf(i,j,k,3), mOf(i,j+1,k,4)]);   // xy face
      if (i < nx && k < nz) pushG([mOf(i,j,k,5), mOf(i,j,k+1,6)]);   // xz face
      if (j < ny && k < nz) pushG([mOf(i,j,k,7), mOf(i,j,k+1,8)]);   // yz face
      if (i < nx && j < ny && k < nz)                                 // body diagonals of cell
        pushG([mOf(i,j,k,9), mOf(i,j,k+1,10), mOf(i,j+1,k,11), mOf(i,j+1,k+1,12)]);
    }

    // ---------- extra (off-grid) nodes ----------
    const extraInfo = [];
    const userExtra = (P.extraNodes || []).map(e => Array.isArray(e) ? e : e.p), faceExtra = faceNodes(P);
    const extra = userExtra.concat(faceExtra);
    const cellRange = (v, o, n) => {               // cells containing coordinate (tolerant)
      const t = (v-o)/h, f = Math.floor(t+1e-9), r = [];
      if (n === 0) return Math.abs(t) < 1e-9 ? [0] : [];   // 2D mode: the single node plane
      if (Math.abs(t-Math.round(t)) < 1e-9) { const c = Math.round(t); if (c-1 >= 0) r.push(c-1); if (c < n) r.push(c); }
      else if (f >= 0 && f < n) r.push(f);
      return r;
    };
    const cellOfExtra = [];
    extra.forEach((p, ei) => {
      if (obs.some(b => inBoxOpen(p, b, tol))) { extraInfo.push({ei, status:'장애물 내부 — 제외'}); return; }
      const dup = X.findIndex(x => norm(sub(x,p)) < tol);
      if (dup >= 0) { extraInfo.push({ei, status:'기존 노드와 중복', node:dup}); return; }
      const ci = cellRange(p[0],origin[0],nx), cj = cellRange(p[1],origin[1],ny), ck = cellRange(p[2],origin[2],nz);
      if (!ci.length || !cj.length || !ck.length) { extraInfo.push({ei, status:'격자 범위 밖 — 제외'}); return; }
      const id = X.length; X.push(p);
      const cells = []; ci.forEach(i => cj.forEach(j => ck.forEach(k => cells.push([i,j,k]))));
      cellOfExtra.push({id, cells}); extraInfo.push({ei, status:'추가됨', node:id, face: ei >= userExtra.length});
    });
    const cellKey = c => c.join(',');
    const extrasInCell = new Map();
    cellOfExtra.forEach(({id, cells}) => cells.forEach(c => {
      const k = cellKey(c); if (!extrasInCell.has(k)) extrasInCell.set(k, []); extrasInCell.get(k).push(id);
    }));
    const removed = new Set();
    extrasInCell.forEach((ids, key) => {
      const [i,j,k] = key.split(',').map(Number);
      const corners = [];
      for (let dk = 0; dk < 2; dk++) for (let dj = 0; dj < 2; dj++) for (let di = 0; di < 2; di++) {
        const g = gidx(i+di, j+dj, k+dk); if (g >= 0) corners.push(g);
      }
      const local = corners.concat(ids);
      const box = {min:[origin[0]+i*h, origin[1]+j*h, origin[2]+k*h], max:[origin[0]+(i+1)*h, origin[1]+(j+1)*h, origin[2]+(k+1)*h]};
      const newM = [];
      ids.forEach(p => local.forEach(q => {
        if (p === q) return;
        if (local.some(m => m !== p && m !== q && pointInSegInterior(X[m], X[p], X[q], tol))) return;
        const id = addM(p, q, 3); if (id >= 0) newM.push(id);
      }));
      const cellMembers = [];
      E.forEach((e, m) => { if (!removed.has(m) && inBoxClosed(X[e[0]], box, tol) && inBoxClosed(X[e[1]], box, tol)) cellMembers.push(m); });
      cellMembers.forEach(m => { if (ids.some(p => !E[m].includes(p) && pointInSegInterior(X[p], X[E[m][0]], X[E[m][1]], tol))) removed.add(m); });
      newM.forEach(m => cellMembers.forEach(n => {
        if (m === n || removed.has(m) || removed.has(n)) return;
        const [a,b] = E[m], [c,d] = E[n];
        if (a===c||a===d||b===c||b===d) return;
        if (segsCross(X[a],X[b],X[c],X[d],tol)) pushG([m,n]);
      }));
    });
    // compact: drop members split by extra nodes
    const remap = new Int32Array(E.length).fill(-1), E2 = [], kind2 = [];
    E.forEach((e, m) => { if (!removed.has(m)) { remap[m] = E2.length; E2.push(e); kind2.push(kind[m]); } });
    const groups2 = [], seen = new Set();
    groups.forEach(g => {
      const h2 = g.map(m => remap[m]).filter(m => m >= 0).sort((a,b) => a-b);
      const key = h2.join(','); if (h2.length > 1 && !seen.has(key)) { seen.add(key); groups2.push(h2); }
    });
    const L = new Float64Array(E2.length), D = new Float64Array(3*E2.length);
    E2.forEach(([a,b], m) => { const d = sub(X[b],X[a]), l = norm(d); L[m] = l; for (let t = 0; t < 3; t++) D[3*m+t] = d[t]/l; });
    const key2m = new Map(); E2.forEach(([a,b], m) => key2m.set(a+','+b, m));
    const memberId = (a, b) => { const k = a < b ? a+','+b : b+','+a; return key2m.has(k) ? key2m.get(k) : -1; };
    const userNodes = new Set(extraInfo.filter(e => e.status === '추가됨' && !e.face).map(e => e.node));
    const GS = { X, E: E2, kind: kind2, L, D, groups: groups2, nGrid, extraInfo, ijk, gidx, memberId, h, userNodes };
    GS.addMembers = pairs => {                     // append members (long / off-grid, kind 4); returns new ids
      const ids = [], add = [];
      pairs.forEach(([a, b]) => { const lo = Math.min(a, b), hi = Math.max(a, b); if (lo === hi || key2m.has(lo+','+hi)) return; key2m.set(lo+','+hi, GS.E.length + add.length); add.push([lo, hi]); });
      if (!add.length) return ids;
      const b0 = GS.E.length, L2 = new Float64Array(b0 + add.length), D2 = new Float64Array(3*(b0 + add.length));
      L2.set(GS.L); D2.set(GS.D);
      add.forEach(([a, b], i) => { const m = b0 + i, d = sub(X[b], X[a]), l = norm(d); GS.E.push([a, b]); GS.kind.push(4); L2[m] = l; for (let t = 0; t < 3; t++) D2[3*m+t] = d[t]/l; ids.push(m); });
      GS.L = L2; GS.D = D2; return ids;
    };
    return GS;
  }

  // ---------- new off-grid joint (crossing points, brace points, moved nodes); idempotent by coordinate ----------
  function addNode(GS, BC, p) {
    const q = 1e-4*GS.h, key = x => x.map(v => Math.round(v/q)).join(',');
    if (!GS.nodeKey) { GS.nodeKey = new Map(); GS.X.forEach((x, j) => GS.nodeKey.set(key(x), j)); }
    const k = key(p); if (GS.nodeKey.has(k)) return GS.nodeKey.get(k);
    const j = GS.X.length; GS.X.push(p.slice()); GS.nodeKey.set(k, j);
    const n = GS.X.length, grow = (a, len) => { const b = new a.constructor(len); b.set(a); return b; };
    BC.fixed = grow(BC.fixed, 3*n); if (BC.flat >= 0) BC.fixed[3*j+BC.flat] = 1;
    BC.loadedAny = grow(BC.loadedAny, n); BC.cases.forEach(c => { c.f = grow(c.f, 3*n); });
    BC.dofMap = grow(BC.dofMap, 3*n); for (let t = 0; t < 3; t++) BC.dofMap[3*j+t] = BC.fixed[3*j+t] ? -1 : BC.nf++;
    return j;
  }

  // ---------- full connectivity candidate pool: node pairs beyond one cell with no node strictly between them ----------
  function candidatePool(P, GS, log = () => {}) {
    const { X, nGrid, ijk } = GS, h = P.grid.h, tol = 1e-6*h, obs = P.obstacles || [], n = X.length;
    let maxL = (P.options?.maxLen || 0)*h;         // optional member length limit [cells → mm], 0 = unlimited
    const extras = []; for (let j = nGrid; j < n; j++) extras.push(j);
    const gcd = (a, b) => { while (b) [a, b] = [b, a % b]; return a; };
    const A = [], B = [];
    const ok = (a, b) => GS.memberId(a, b) < 0 && !(maxL > 0 && norm(sub(X[b], X[a])) > maxL + tol) && !obs.some(bx => segHitsBoxInterior(X[a], X[b], bx, tol));
    for (let a = 0; a < nGrid; a++) { const [i0, j0, k0] = ijk[a]; for (let b = a+1; b < nGrid; b++) {
      const di = Math.abs(ijk[b][0]-i0), dj = Math.abs(ijk[b][1]-j0), dk = Math.abs(ijk[b][2]-k0);
      if (Math.max(di, dj, dk) <= 1 || gcd(gcd(di, dj), dk) !== 1 || !ok(a, b)) continue;   // gcd = 1: no grid point in between
      if (extras.some(e => pointInSegInterior(X[e], X[a], X[b], tol))) continue;
      A.push(a); B.push(b);
    } }
    if (extras.length*n*n > 3e8 && !(maxL > 0 && maxL <= 3*h)) { maxL = 3*h; log(`추가 노드가 많아 추가 노드의 장거리 후보는 3셀 이내로 제한합니다`); }
    extras.forEach(e => { for (let b = 0; b < n; b++) {
      if (b === e || (b >= nGrid && b < e) || !ok(e, b)) continue;
      let hit = false; for (let c = 0; c < n && !hit; c++) if (c !== e && c !== b && pointInSegInterior(X[c], X[e], X[b], tol)) hit = true;
      if (!hit) { A.push(Math.min(e, b)); B.push(Math.max(e, b)); }
    } });
    return { a: Int32Array.from(A), b: Int32Array.from(B) };
  }

  // ---------- boundary conditions & load cases ----------
  function buildBC(P, GS) {
    const { X } = GS, n = X.length, tol = 1e-6*P.grid.h;
    const fixed = new Uint8Array(3*n), loadedAny = new Uint8Array(n), supNodes = new Set(), warn = [];
    const flat = [P.grid.nx, P.grid.ny, P.grid.nz].findIndex(v => v === 0);   // 2D mode: out-of-plane axis (−1 = 3D)
    if (flat >= 0) for (let j = 0; j < n; j++) fixed[3*j+flat] = 1;           // plane truss, braced out of plane
    (P.supports || []).forEach((s, si) => {
      let cnt = 0;
      X.forEach((p, j) => { if (inBoxClosed(p, s, tol)) { cnt++; supNodes.add(j); for (let a = 0; a < 3; a++) if (s.fix[a]) fixed[3*j+a] = 1; } });
      if (!cnt) warn.push(`지지 영역 ${si+1}: 포함 노드 없음`);
    });
    const caseIds = [...new Set((P.loads || []).map(l => l.case || 1))].sort((a,b) => a-b), deck = new Set();   // deck: forced load-surface members
    const cases = caseIds.map(c => ({ name: 'LC'+c, f: new Float64Array(3*n) }));
    (P.loads || []).forEach((l, li) => {
      const nodes = []; X.forEach((p, j) => { if (inBoxClosed(p, l, tol)) nodes.push(j); });
      if (!nodes.length) { warn.push(`하중 영역 ${li+1}: 포함 노드 없음`); return; }
      const f = cases[caseIds.indexOf(l.case || 1)].f;
      if (flat >= 0 && Math.abs(l.F[flat]) > 0) warn.push(`하중 영역 ${li+1}: 2D 모드라 평면 밖(${'xyz'[flat]}) 성분 ${l.F[flat]} N은 무시합니다`);
      const trib = t => {                          // tributary length along axis t (half spacing to each neighbouring load-node coordinate)
        const cs = [...new Set(nodes.map(j => Math.round(X[j][t]/tol)*tol))].sort((a, b) => a - b);
        return v => { if (cs.length < 2) return 1; const k = cs.findIndex(c => Math.abs(c - v) <= 2*tol); return ((k > 0 ? v - cs[k-1] : 0) + (k < cs.length - 1 ? cs[k+1] - v : 0))/2; };
      };
      const tr = [0, 1, 2].map(trib), w = nodes.map(j => tr.reduce((acc, g, t) => acc*g(X[j][t]), 1)), ws = w.reduce((a, b) => a + b, 0) || 1;
      nodes.forEach((j, k) => { loadedAny[j] = 1; for (let a = 0; a < 3; a++) if (a !== flat) f[3*j+a] += l.F[a]*w[k]/ws; });   // distributed load: tributary weights
      if (l.deck === true && nodes.length > 1) {   // load-surface members (JSON opt-in only; off in the UI): one-cell edges + one radial face diagonal per cell face (mirror-consistent)
        const ns = new Set(nodes), c = [0, 1, 2].map(t => 0.5*(l.min[t] + l.max[t]));
        GS.E.forEach(([a, b], m) => { if (GS.kind[m] === 0 && ns.has(a) && ns.has(b)) deck.add(m); });
        GS.groups.forEach(g => {
          if (g.length !== 2 || !g.every(m => GS.kind[m] === 1 && GS.E[m].every(j => ns.has(j)))) return;
          const pick = g.find(m => { const [a, b] = GS.E[m], d = sub(X[b], X[a]), mid = X[a].map((v, t) => 0.5*(v + X[b][t]) - c[t]), ax = [0, 1, 2].filter(t => Math.abs(d[t]) > tol);
            const pv = mid[ax[0]]*mid[ax[1]]; return Math.abs(pv) > tol*tol && Math.sign(pv) === Math.sign(d[ax[0]]*d[ax[1]]); });
          if (pick !== undefined) deck.add(pick);
        });
      }
    });
    const Ft = [0, 1, 2].map(t => { let v = 0; if (cases.length) for (let j = 0; j < n; j++) v += cases[0].f[3*j+t]; return v; });
    const Fn = Math.hypot(...Ft);
    const pert = P.options?.perturb || 0;          // lateral perturbation ratio (0 = off)
    if (pert > 0 && cases.length && Fn > 0) {      // load cases perpendicular to the dominant load axis
      const base = cases[0].f, mag = new Float64Array(n), dom = Ft.map(Math.abs).indexOf(Math.max(...Ft.map(Math.abs)));
      for (let j = 0; j < n; j++) mag[j] = Math.hypot(base[3*j], base[3*j+1], base[3*j+2]);
      const signs = Math.abs(P.material.sigmaT - P.material.sigmaC) > 1e-9 ? [1, -1] : [1];
      [0, 1, 2].filter(t => t !== dom).forEach(t => signs.forEach(sg => {
        const f = new Float64Array(3*n);
        for (let j = 0; j < n; j++) f[3*j+t] = sg*pert*mag[j];
        cases.push({ name: 'P' + (sg > 0 ? '+' : '−') + 'xyz'[t], f, perturb: true });
      }));
    }
    const share = [];                              // reaction share rows: Σ ĝ·(Bq − f) ≥ ρ|F| on each support region (case 1)
    const rho = P.options?.supportShare || 0;
    if (rho > 0 && Fn > 0) {
      const g = Ft.map(v => -v/Fn);
      (P.supports || []).forEach((sp, si) => {
        if ([0,1,2].reduce((acc, t) => acc + (sp.fix[t] ? g[t]*g[t] : 0), 0) < 0.25) { warn.push(`지지 영역 ${si+1}: 하중 방향 구속이 약해 반력 분담 조건을 적용하지 않음`); return; }
        const dofs = []; let rhs = rho*Fn;
        X.forEach((p, j) => { if (inBoxClosed(p, sp, tol)) for (let t = 0; t < 3; t++) if (sp.fix[t] && Math.abs(g[t]) > 1e-12) { dofs.push([3*j+t, g[t]]); rhs += g[t]*cases[0].f[3*j+t]; } });
        if (dofs.length) share.push({ dofs, rhs, region: si });
      });
    }
    const dofMap = new Int32Array(3*n).fill(-1); let nf = 0;
    for (let i = 0; i < 3*n; i++) if (!fixed[i]) dofMap[i] = nf++;
    return { fixed, dofMap, nf, cases, loadedAny, supNodes, warn, share, flat, deck };
  }

  // member is part of the design if its area or any load-case force is not negligible
  function presence(r, alpha) {
    const amax = amaxOf(r.a), qmax = r.q.map(qk => qk.reduce((s, v) => Math.max(s, Math.abs(v)), 0));
    return m => r.a[m] > alpha*amax || r.q.some((qk, k) => qk[m] !== 0 && Math.abs(qk[m]) > alpha*qmax[k]);
  }

  // ---------- mirror symmetry about mid-planes of the design domain ----------
  function detectSymmetry(P, GS, BC) {            // mirror planes under which the discretised model is invariant
    const { origin, h, nx, ny, nz } = P.grid, N = [nx, ny, nz], q = 1e-4*h, X = GS.X, n = X.length;
    const key = p => p.map(v => Math.round(v/q)).join(',');
    const idx = new Map(); X.forEach((p, j) => idx.set(key(p), j));
    let fmax = 0; BC.cases.forEach(c => c.f.forEach(v => { fmax = Math.max(fmax, Math.abs(v)); }));
    return [0, 1, 2].map(t => {
      if (N[t] === 0) return 0;                    // flat axis of a 2D model: nothing to mirror
      const c2 = 2*origin[t] + N[t]*h, nm = new Int32Array(n);
      for (let j = 0; j < n; j++) { const r = X[j].slice(); r[t] = c2 - r[t]; const k = idx.get(key(r)); if (k === undefined) return 0; nm[j] = k; }
      for (let j = 0; j < n; j++) for (let k = 0; k < 3; k++) {
        if (BC.fixed[3*j+k] !== BC.fixed[3*nm[j]+k]) return 0;
        for (const c of BC.cases) if (Math.abs(c.f[3*nm[j]+k] - (k === t ? -1 : 1)*c.f[3*j+k]) > 1e-9*Math.max(fmax, 1)) return 0;
      }
      if (GS.E.some(([a, b]) => GS.memberId(nm[a], nm[b]) < 0)) return 0;   // obstacles / missing members break it
      const sig = sr => sr.dofs.map(([i]) => i).sort((x, y) => x - y).join(',');
      const rows = new Set((BC.share || []).map(sig));
      if ((BC.share || []).some(sr => !rows.has(sr.dofs.map(([i]) => 3*nm[Math.floor(i/3)] + i % 3).sort((x, y) => x - y).join(',')))) return 0;
      return 1;
    });
  }

  function symmetryMaps(P, GS) {
    const sym = P.options?.symmetry || [0, 0, 0]; if (!sym.some(Boolean)) return null;
    const { origin, h, nx, ny, nz } = P.grid, N = [nx, ny, nz], q = 1e-4*h;
    const key = p => p.map(v => Math.round(v/q)).join(',');
    const idx = new Map(); GS.X.forEach((p, j) => idx.set(key(p), j));
    const planes = [0, 1, 2].filter(t => sym[t]);
    const nodeMir = planes.map(t => { const c2 = 2*origin[t] + N[t]*h; return GS.X.map(p => { const r = p.slice(); r[t] = c2 - p[t]; const j = idx.get(key(r)); return j === undefined ? -1 : j; }); });
    const memMir = nodeMir.map(nm => GS.E.map(([a, b]) => nm[a] < 0 || nm[b] < 0 ? -1 : GS.memberId(nm[a], nm[b])));
    const par = Int32Array.from(GS.E.keys()), find = x => { while (par[x] !== x) x = par[x] = par[par[x]]; return x; };
    let missing = 0;
    memMir.forEach(mm => mm.forEach((o, m) => { if (o < 0) { missing++; return; } const a = find(m), b = find(o); if (a !== b) par[Math.max(a, b)] = Math.min(a, b); }));
    const orb = new Map(); GS.E.forEach((_, m) => { const r = find(m); if (!orb.has(r)) orb.set(r, []); orb.get(r).push(m); });
    const orbits = [...orb.values()].filter(o => o.length > 1), orbitOf = new Int32Array(GS.E.length).fill(-1);
    orbits.forEach((o, oi) => o.forEach(m => orbitOf[m] = oi));
    const out = { planes, nodeMir, memMir, orbits, orbitOf, missing };
    if (P.options?.symExact && !missing && nodeMir.every(nm => nm.every(j => j >= 0))) addReduction(out, GS);
    return out;
  }
  function addReduction(sym, GS) {                 // exact model symmetry: member groups + representative joints for the reduced LP
    const n = GS.X.length, rep = Int32Array.from({ length: n }, (_, j) => j);
    for (let pass = 0; pass < 3; pass++) sym.nodeMir.forEach(nm => { for (let j = 0; j < n; j++) { const r = Math.min(rep[j], rep[nm[j]]); rep[j] = rep[nm[j]] = r; } });
    const size = new Int32Array(n); for (let j = 0; j < n; j++) size[rep[j]]++;
    const groups = sym.orbits.slice(); GS.E.forEach((_, m) => { if (sym.orbitOf[m] < 0) groups.push([m]); });
    sym.exact = true; sym.nodeRep = rep; sym.nodeSize = size; sym.groups = groups;
  }

  // ---------- LP / MILP ----------
  function cooToCsc(nRows, nCols, R, C, V) {
    const starts = new Int32Array(nCols+1);
    for (let k = 0; k < C.length; k++) starts[C[k]+1]++;
    for (let j = 0; j < nCols; j++) starts[j+1] += starts[j];
    const pos = starts.slice(0, nCols), indices = new Int32Array(C.length), values = new Float64Array(C.length);
    for (let k = 0; k < C.length; k++) { const p = pos[C[k]]++; indices[p] = R[k]; values[p] = V[k]; }
    const s2 = new Int32Array(nCols+1), i2 = new Int32Array(C.length), v2 = new Float64Array(C.length); let w = 0;
    for (let j = 0; j < nCols; j++) {               // merge duplicate (row, col) entries — HiGHS rejects them
      const seen = new Map();
      for (let p = starts[j]; p < starts[j+1]; p++) { const r = indices[p]; if (seen.has(r)) v2[seen.get(r)] += values[p]; else { seen.set(r, w); i2[w] = r; v2[w++] = values[p]; } }
      s2[j+1] = w;
    }
    return { format: 'csc', numRows: nRows, numCols: nCols, starts: s2, indices: i2.slice(0, w), values: v2.slice(0, w) };
  }

  function addShareRows(BC, GS, qCols, addRow, put, INF) {   // reaction share on support regions (case 1)
    if (!BC.share || !BC.share.length) return [];
    const byDof = new Map(), rows = [];
    BC.share.forEach(sr => { const r = addRow(sr.rhs, INF); rows.push(r); sr.dofs.forEach(([i, g]) => { if (!byDof.has(i)) byDof.set(i, []); byDof.get(i).push([r, g]); }); });
    GS.E.forEach(([a, b], m) => {
      for (let t = 0; t < 3; t++) {
        const d = GS.D[3*m+t]; if (Math.abs(d) < 1e-15) continue;
        if (!qCols[m]) continue;
        [[a, -d], [b, d]].forEach(([j, bv]) => (byDof.get(3*j+t) || []).forEach(([r, g]) => qCols[m].forEach(([c, sg]) => put(r, c, g*bv*sg))));
      }
    });
    return rows;
  }

  function solveLayout(highs, P, GS, BC, { yMembers, groupsAct, aBar, banned, aLo, central } = {}) {
    const { E, L, D } = GS, { dofMap, nf, cases } = BC, b = E.length, K = cases.length;
    const sT = P.material.sigmaT, sC = P.material.sigmaC, s = P.options?.jointCost || 0;
    const INF = highs.infinity, R = [], C = [], V = [], cost = [], lo = [], up = [], rL = [], rU = [], integ = [];
    let nRow = 0;
    const addCol = (c, l, u, it = 0) => { cost.push(c); lo.push(l); up.push(u); integ.push(it); return cost.length-1; };
    const addRow = (l, u) => { rL.push(l); rU.push(u); return nRow++; };
    const put = (r, c, v) => { R.push(r); C.push(c); V.push(v); };
    const eqRow = [];                              // eqRow[k][dof] = row index
    const sym = GS.sym, red = !!(sym && sym.exact && sym.nodeRep && sym.nodeRep.length === GS.X.length);   // reduced: one column set per orbit, equilibrium on representative joints
    for (let k = 0; k < K; k++) { const rows = new Int32Array(nf).fill(-1); for (let i = 0; i < 3*GS.X.length; i++) if (dofMap[i] >= 0 && (!red || sym.nodeRep[Math.floor(i/3)] === Math.floor(i/3))) rows[dofMap[i]] = addRow(cases[k].f[i], cases[k].f[i]); eqRow.push(rows); }
    const putB = (k, m, col, sign) => {
      const [a, bb] = E[m];
      for (let t = 0; t < 3; t++) {
        const da = dofMap[3*a+t], db = dofMap[3*bb+t], d = D[3*m+t];
        if (Math.abs(d) < 1e-15) continue;
        if (da >= 0 && eqRow[k][da] >= 0) put(eqRow[k][da], col, -sign*d);
        if (db >= 0 && eqRow[k][db] >= 0) put(eqRow[k][db], col, sign*d);
      }
    };
    const areaTerms = [];                          // areaTerms[m] = [[col, coef], ...] so that a_m = sum
    const ub = m => banned && banned.has(m) ? 0 : INF;  // banned members are fixed to zero area
    const qCols = new Array(b).fill(null);          // qCols[m] = [[col, sign]] giving q of case 1
    for (let m = 0; m < b; m++) areaTerms.push([]);
    const lo0 = m => aLo && aLo.has(m) ? aLo.get(m) : 0;      // per-member minimum area (Map)
    const groups = red ? sym.groups : null, nG = red ? groups.length : b, gm = gi => red ? groups[gi] : [gi];
    for (let gi = 0; gi < nG; gi++) {
      const g = gm(gi); if (banned && g.some(m => banned.has(m))) continue;   // banned: no columns at all
      const Lg = g.reduce((acc, m) => acc + L[m] + s, 0);
      if (K === 1 && !aLo) {
        const cp = addCol(Lg/sT, 0, INF), cm = addCol(Lg/sC, 0, INF);
        g.forEach(m => { putB(0, m, cp, 1); putB(0, m, cm, -1); areaTerms[m] = [[cp, 1/sT], [cm, 1/sC]]; qCols[m] = [[cp, 1], [cm, -1]]; });
      } else {
        const ca = addCol(Lg, Math.max(...g.map(lo0)), INF);
        g.forEach(m => { areaTerms[m] = [[ca, 1]]; });
        for (let k = 0; k < K; k++) {
          const cq = addCol(0, -INF, INF); g.forEach(m => { putB(k, m, cq, 1); if (k === 0) qCols[m] = [[cq, 1]]; });
          const r1 = addRow(-INF, 0); put(r1, cq, 1); put(r1, ca, -sT);
          const r2 = addRow(-INF, 0); put(r2, cq, -1); put(r2, ca, -sC);
        }
      }
    }
    const shareRows = addShareRows(BC, GS, qCols, addRow, put, INF);
    if (sym && !red) sym.orbits.forEach(o => o.slice(1).forEach(m => { if (!areaTerms[m].length && !areaTerms[o[0]].length) return; const r = addRow(0, 0); areaTerms[m].forEach(([c, v]) => put(r, c, v)); areaTerms[o[0]].forEach(([c, v]) => put(r, c, -v)); }));
    const yCol = new Map();
    (yMembers || []).forEach(m => {
      const cy = addCol(0, 0, 1, 1); yCol.set(m, cy);
      const r = addRow(-INF, 0); areaTerms[m].forEach(([c, v]) => put(r, c, v)); put(r, cy, -aBar);
    });
    (groupsAct || []).forEach(g => { const r = addRow(-INF, 1); g.forEach(m => put(r, yCol.get(m), 1)); });
    if (sym) sym.orbits.forEach(o => { const ys = o.filter(m => yCol.has(m)); ys.slice(1).forEach(m => { const r = addRow(0, 0); put(r, yCol.get(m), 1); put(r, yCol.get(ys[0]), -1); }); });
    const nCol = cost.length, mip = yCol.size > 0;
    const model = highs.createModel({
      numCols: nCol, numRows: nRow, colCost: cost, colLower: lo, colUpper: up, rowLower: rL, rowUpper: rU,
      matrix: cooToCsc(nRow, nCol, R, C, V), ...(mip ? { integrality: new Int32Array(integ) } : {}) });
    try {
      model.options.set({ output_flag: false, time_limit: P.options?.timeLimit || 60,
        ...(mip ? { mip_rel_gap: P.options?.mipGap ?? 1e-3 } : { solver: P.options?.lpSolver || 'ipm', ...(central ? { run_crossover: 'off' } : {}) }) });
      const t0 = Date.now(); model.run(); const dt = Date.now()-t0;
      const MS = highs.constants.modelStatus, st = model.getModelStatus();
      const ok = st === MS.optimal || (st === MS.timeLimit && mip && Number.isFinite(model.getObjectiveValue()));
      const sol = ok ? model.getSolution() : null;
      if (!sol || !sol.colValue || sol.colValue.length !== nCol) return { ok: false, status: st, dt, timeLimited: st === MS.timeLimit };
      const x = sol.colValue, a = new Float64Array(b), q = Array.from({length: K}, () => new Float64Array(b));
      for (let m = 0; m < b; m++) {
        if (!areaTerms[m].length) continue;
        a[m] = areaTerms[m].reduce((acc, [c, v]) => acc + v*x[c], 0);
        if (K === 1 && !aLo) q[0][m] = x[areaTerms[m][0][0]] - x[areaTerms[m][1][0]];
        else for (let k = 0; k < K; k++) q[k][m] = x[areaTerms[m][0][0]+1+k];   // columns: [a, q_1..q_K]
      }
      let U = null;                                // virtual displacements per case from the equilibrium duals (LP only)
      if (!mip && sol.rowDual && sol.rowDual.length === nRow) {
        const yd = sol.rowDual, n3 = 3*GS.X.length;
        U = eqRow.map(rows => { const u = new Float64Array(n3); for (let i = 0; i < n3; i++) if (dofMap[i] >= 0 && rows[dofMap[i]] >= 0) u[i] = yd[rows[dofMap[i]]]; return u; });
        if (red) U.forEach(u => {                  // full-model duals: u_j = R_j u_rep / |joint orbit|
          const X = GS.X, nr = sym.nodeRep, ns = sym.nodeSize, base = Float64Array.from(u);
          for (let j = 0; j < X.length; j++) { const r = nr[j]; for (let t = 0; t < 3; t++) { const flip = Math.abs(X[j][t] - X[r][t]) > 1e-9*(1 + Math.abs(X[r][t])) ? -1 : 1; u[3*j+t] = BC.fixed[3*j+t] ? 0 : flip*base[3*r+t]/ns[r]; } }
        });
        BC.share.forEach((sr, i) => { const y = yd[shareRows[i]]; if (y) sr.dofs.forEach(([dof, g]) => { U[0][dof] += y*g; }); });
      }
      return { ok: true, status: st, dt, obj: model.getObjectiveValue(), a, q, U, nCol, nRow, mip, timeLimited: st === highs.constants.modelStatus.timeLimit };
    } finally { model.dispose(); }
  }

  function optimize(highs, P, GS, BC, log = () => {}) {
    const alpha = P.options?.areaCut ?? 1e-3;      // relative area cutoff
    const noCross = P.options?.noCrossing !== false;
    const maxIt = P.options?.maxCrossIter || 20;
    let method = P.options?.crossMethod || 'auto'; // 'milp' exact | 'fix' heuristic | 'auto'
    if (method === 'auto') method = BC.cases.length === 1 && GS.E.length <= (P.options?.milpMaxMembers || 20000) ? 'milp' : 'fix';   // very large ground structures: heuristic
    let res = solveLayout(highs, P, GS, BC);
    if (!res.ok) return { ok: false, msg: `LP 실패 (status ${res.status})`, iters: [] };
    log(`LP: 변수 ${res.nCol}, 제약 ${res.nRow}, 목적함수 ${res.obj.toExponential(4)}, ${res.dt} ms`);
    const iters = [{ kind: 'LP', obj: res.obj, dt: res.dt }], lpObj = res.obj;
    const sym = GS.sym, gKey = new Map(); GS.groups.forEach((g, gi) => gKey.set([...g].sort((x, y) => x - y).join(','), gi));
    const mirrorGroups = gi => {                   // indices of the mirror images of a crossing group
      if (!sym) return [gi];
      const out = new Set([gi]);
      for (let pass = 0; pass < 3; pass++) [...out].forEach(i => sym.memMir.forEach(mm => {
        const mg = GS.groups[i].map(m => mm[m]); if (mg.some(m => m < 0)) return;
        const j = gKey.get(mg.sort((x, y) => x - y).join(',')); if (j !== undefined) out.add(j);
      }));
      return [...out];
    };
    const closeBans = bset => { if (sym) for (let pass = 0; pass < 3; pass++) [...bset].forEach(m => sym.memMir.forEach(mm => { if (mm[m] >= 0) bset.add(mm[m]); })); };
    const activeGroups = (r, skip) => {
      const on = presence(r, alpha);
      return GS.groups.map((g, gi) => [g, gi]).filter(([g, gi]) => !(skip && skip.has(gi)) && g.filter(on).length > 1);
    };
    let poolBan = null;                            // large ground structure: crossing MILP / fix on the LP support grown by one cell ring
    if (noCross && GS.E.length > (P.options?.milpPoolMax || 3000) && activeGroups(res).length) {
      const on = presence(res, alpha), nodes = new Set(); GS.E.forEach(([a, b], m) => { if (on(m)) { nodes.add(a); nodes.add(b); } });
      const ring = new Set(nodes); GS.E.forEach(([a, b], m) => { if (GS.kind[m] < 4 && (nodes.has(a) || nodes.has(b))) { ring.add(a); ring.add(b); } });
      if (sym) for (let pass = 0; pass < 3; pass++) [...ring].forEach(j => sym.nodeMir.forEach(nm => { if (nm[j] >= 0) ring.add(nm[j]); }));
      poolBan = new Set(); GS.E.forEach(([a, b], m) => { if (!(ring.has(a) && ring.has(b))) poolBan.add(m); });
      log(`교차 처리 후보: LP 사용 노드 + 1셀 이웃 사이 부재 ${GS.E.length - poolBan.size}개 (전체 ${GS.E.length})`);
    }
    if (noCross && method === 'milp') {
      const aBar = 10*amaxOf(res.a), Gset = new Set(), Gact = [];
      for (let it = 0; it < maxIt; it++) {
        const fresh = activeGroups(res, Gset); if (!fresh.length) break;
        fresh.forEach(([, gi]) => mirrorGroups(gi).forEach(j => { if (!Gset.has(j)) { Gset.add(j); Gact.push(GS.groups[j]); } }));
        const ym = [...new Set(Gact.flat())];
        log(`교차 제약(MILP) 반복 ${it+1}: 활성 교차군 +${fresh.length} (누적 ${Gact.length}), 이진변수 ${ym.length}`);
        const r2 = solveLayout(highs, P, GS, BC, { yMembers: ym, groupsAct: Gact, aBar, banned: poolBan || undefined });
        if (!r2.ok) { log(`  MILP 해 없음${r2.timeLimited ? ' (시간 제한)' : ''} → 고정 휴리스틱으로 전환`); method = 'fix'; break; }
        res = r2; iters.push({ kind: 'MILP', obj: res.obj, dt: res.dt, timeLimited: res.timeLimited });
        log(`  목적함수 ${res.obj.toExponential(4)}, ${res.dt} ms${res.timeLimited ? ' (시간 제한 도달, 최적성 미보장)' : ''}`);
      }
    }
    if (noCross && method === 'fix') {
      const banned = new Set(poolBan || []);
      for (let it = 0; it < maxIt; it++) {
        const act = activeGroups(res); if (!act.length) break;
        act.forEach(([g]) => { const keep = g.reduce((p, m) => res.a[m] > res.a[p] + 1e-9 ? m : p, g[0]); g.forEach(m => { if (m !== keep) banned.add(m); }); });
        closeBans(banned);
        const r2 = solveLayout(highs, P, GS, BC, { banned });
        if (!r2.ok) return { ok: false, msg: `교차 고정 LP 실패 (status ${r2.status})`, iters };
        res = r2; iters.push({ kind: 'FIX', obj: res.obj, dt: res.dt });
        log(`교차 고정 반복 ${it+1}: 교차군 ${act.length}개 정리, 제외 부재 누적 ${banned.size}, 목적함수 ${res.obj.toExponential(4)}, ${res.dt} ms`);
      }
    }
    return { ok: true, a: res.a, q: res.q, obj: res.obj, lpObj, iters, method };
  }



  // ---------- variable section with minimum area: a_i = 0 or a_i ≥ A_min (LP bounds, no binaries) ----------
  function solveMinArea(highs, P, GS, BC, log = () => {}) {
    const o = P.options, alpha = o.areaCut ?? 1e-3, noCross = o.noCrossing !== false, sym = GS.sym;
    const rc = o.ratioC ?? 0.2, rt = o.ratioT ?? 0.01;  // section ratio limits: general/compression, tension rods
    const base = optimize(highs, P, GS, BC, log);
    if (!base.ok) return base;
    const amax0 = amaxOf(base.a), thrC = o.minA > 0 ? o.minA : rc*amax0, thrT = o.allowRods ? Math.min(thrC, rt*amax0) : thrC;   // no rods: one threshold for all
    const lb = base.lpObj, orbit = m => sym && sym.orbitOf[m] >= 0 ? sym.orbits[sym.orbitOf[m]] : [m];
    const close = set => { if (sym) for (let k = 0; k < 3; k++) [...set].forEach(m => sym.memMir.forEach(mm => { if (mm[m] >= 0) set.add(mm[m]); })); return set; };
    const qTol = r => 1e-9*Math.max(1, ...r.q.map(qk => qk.reduce((s, v) => Math.max(s, Math.abs(v)), 0)));
    const isRod = (r, m) => { const t = qTol(r); return r.q.every(qk => qk[m] >= -t); };   // tension (or zero) in every load case
    const thrOf = (r, m) => isRod(r, m) ? thrT : thrC;
    log(`단면 하한: ${thrC.toFixed(3)} mm² (${o.minA > 0 ? 'A_min' : (100*rc).toFixed(0) + '% A_max'})${o.allowRods ? `, 텐션 로드 ${thrT.toFixed(3)} mm² (${(100*rt).toFixed(1)}% A_max)` : ' — 텐션 로드 없음: 하한 미만 부재는 제거 후 우회, 불가하면 하한으로'}`);
    const T = 1000*(o.timeLimit || 60);
    const strategy = (mode, until) => {            // mode 'support': current members only; 'nodes': any member between used nodes
      const keep = new Map();                     // member → held lower bound
      let banned = new Set(), res = base, lps = 0;
      const solveX = bn => {                      // LP with bans/lower bounds; crossings that reappear are fixed by banning the thinner member
        bn = new Set([...bn].filter(m => !keep.has(m)));
        let r = solveLayout(highs, P, GS, BC, { banned: bn, aLo: keep }); lps++;
        if (globalThis.GSTO_DEBUG && !r.ok && keep.size) console.log('minA solveX fail status', r.status, 'keep', keep.size, 'ban', bn.size);
        for (let k = 0; r.ok && noCross && k < 20; k++) {
          const on = presence(r, alpha), act = GS.groups.filter(g => g.filter(on).length > 1);
          if (!act.length) break;
          act.forEach(g => { const kp = g.find(m => keep.has(m)) ?? g.reduce((p, m) => r.a[m] > r.a[p] + 1e-9 ? m : p, g[0]); g.forEach(m => { if (m !== kp && !keep.has(m)) bn.add(m); }); });
          close(bn); [...bn].forEach(m => { if (keep.has(m)) bn.delete(m); });
          r = solveLayout(highs, P, GS, BC, { banned: bn, aLo: keep }); lps++;
        }
        return { r, bn };
      };
      const hold = (m, v) => orbit(m).forEach(x => keep.set(x, Math.max(keep.get(x) || 0, v)));
      (BC.deck || new Set()).forEach(m => hold(m, thrC));   // load-surface members always present
      GS.groups.forEach(g => { if (g.some(m => keep.has(m))) g.forEach(m => { if (!keep.has(m)) banned.add(m); }); });
      { const on = presence(res, 1e-9), used = new Set();   // keep every force-carrying member: tiny ones are pruned below, not cut here
        GS.E.forEach(([a, b], m) => { if (on(m)) { used.add(a); used.add(b); } });
        GS.E.forEach(([a, b], m) => { if (mode === 'support' ? !on(m) : !(used.has(a) && used.has(b))) banned.add(m); });
        if (noCross) GS.groups.forEach(g => { const kp = g.find(on); if (kp !== undefined) g.forEach(m => { if (m !== kp) banned.add(m); }); });
        close(banned); }
      let removed = 0;
      const smallList = r => { const on = presence(r, alpha); return GS.E.map((_, m) => m).filter(m => on(m) && !keep.has(m) && r.a[m] < thrOf(r, m)*(1 - 1e-6)).sort((u, v) => r.a[u] - r.a[v]); };
      while (Date.now() < until) {
        const small = smallList(res);
        if (!small.length) break;
        const tie = k => { let e = k; while (e < small.length && res.a[small[e]] <= res.a[small[k-1]]*(1 + 1e-6) + 1e-12) e++; return e; };   // never split equal-area members
        let bs = tie(Math.max(1, Math.ceil(0.3*small.length)));
        for (;;) {
          const trial = close(new Set([...banned, ...small.slice(0, bs).flatMap(orbit)]));
          const x = solveX(trial);
          if (x.r.ok) { res = x.r; banned = x.bn; removed += bs; break; }
          const half = tie(Math.ceil(bs/2)); if (half < bs) { bs = half; continue; }
          small.slice(0, tie(1)).forEach(m => hold(m, thrOf(res, m)));   // cannot be removed: hold it (and its equal-area peers) at the threshold
          const y = solveX(banned); if (y.r.ok) { res = y.r; banned = y.bn; }
          break;
        }
      }
      for (let pass = 0; pass < 6; pass++) {      // held rods that went into compression, and leftovers at the time limit
        let changed = false;
        keep.forEach((v, m) => { if (v < thrC && !isRod(res, m)) { hold(m, thrC); changed = true; } });
        smallList(res).forEach(m => { hold(m, thrOf(res, m)); changed = true; });
        if (!changed) break;
        const x = solveX(banned); if (!x.r.ok) break; res = x.r; banned = x.bn;
      }
      if (globalThis.GSTO_DEBUG) { const on = presence(res, alpha); let below = 0, held = 0, heldBelow = 0; GS.E.forEach((_, m) => { if (on(m) && res.a[m] < thrC*(1 - 1e-6)) below++; if (keep.has(m)) { held++; if (res.a[m] < keep.get(m)*(1 - 1e-6)) heldBelow++; } }); console.log('minA dbg below', below, 'held', held, 'heldBelow', heldBelow, 'exact', !!(GS.sym && GS.sym.exact)); }
      log(`  전략 ${mode === 'support' ? '현재 부재 안' : '사용 노드 사이'}: 제거 ${removed}개, 하한 유지 ${keep.size}개, LP ${lps}회, 체적 ${res.obj.toExponential(4)}`);
      return res;
    };
    const t0 = Date.now(), only = o.minAOnly === 'support';   // geometry stage: current members only
    const r1 = strategy('support', only ? t0 + T : t0 + T/2), r2 = only ? r1 : strategy('nodes', Math.max(Date.now() + T/4, t0 + T));
    const res = r2.obj < r1.obj ? r2 : r1;
    const gap = res.obj > 0 ? (res.obj - lb)/res.obj : 0;
    log(`단면 하한 적용: 체적 ${res.obj.toExponential(4)} (LP 하한 대비 +${(100*gap).toFixed(1)}%)`);
    return { ok: true, a: res.a, q: res.q, obj: res.obj, lpObj: lb, method: 'minA', gap, thrC, thrT, amax0, iters: [{ kind: 'MINA', obj: res.obj }] };
  }

  // ---------- minimum section by round-up: free LP → a_i ≥ A_min on the used members → drop thin members that do not pay ----------
  function solveRoundUp(highs, P, GS, BC, log = () => {}) {
    const o = P.options || {}, sym = GS.sym, s = o.jointCost || 0;
    const base = optimize(highs, P, GS, BC, log); if (!base.ok) return base;
    const cat = (o.catalog || []).filter(v => v > 0).sort((x, y) => x - y);   // section catalogue [mm²], ascending (optional)
    const amax0 = amaxOf(base.a);
    const Amin = o.minA > 0 ? o.minA : cat.length ? cat[0] : (o.ratioC ?? 0.2)*amax0;   // minimum section [mm²]
    const orbit = m => sym && sym.orbitOf[m] >= 0 ? sym.orbits[sym.orbitOf[m]] : [m];
    const on = presence(base, 1e-9); let S = new Set(); GS.E.forEach((_, m) => { if (on(m)) S.add(m); });
    const deck = BC.deck || new Set(), mates = new Set();   // load-surface members are always kept; their crossing partners never
    GS.groups.forEach(g => { if (g.some(m => deck.has(m))) g.forEach(m => { if (!deck.has(m)) mates.add(m); }); });
    deck.forEach(m => S.add(m)); mates.forEach(m => S.delete(m));
    const naive = [...S].reduce((v, m) => v + (GS.L[m] + s)*Math.max(Amin, base.a[m]), 0);
    let lps = 0;
    const lpOn = set => { const banned = new Set(), aLo = new Map(); GS.E.forEach((_, m) => { if (set.has(m)) aLo.set(m, Amin); else banned.add(m); }); lps++; return solveLayout(highs, P, GS, BC, { banned, aLo }); };
    let res = lpOn(S);
    if (!res.ok) return { ok: false, msg: `최소 단면 재최적화 LP 실패 (status ${res.status})${deck.size ? ' — 하중면 부재와 겹치는 대각재를 뺀 뒤 평형이 성립하지 않습니다' : ''}` };
    const thinCount = r => [...S].filter(m => r.a[m] <= Amin*(1 + 1e-6)).length;
    const V1 = res.obj, n1 = S.size, t1 = thinCount(res), T = Date.now() + 1000*(o.timeLimit || 30), tried = new Set(), sMin = Math.min(P.material.sigmaT, P.material.sigmaC);
    let removed = 0; const R = new Set();          // R: members taken out of the ground structure
    const near = () => {                            // members between the joints of S and their one-cell neighbours (rerouting pool)
      const nodes = new Set(); S.forEach(m => GS.E[m].forEach(j => nodes.add(j)));
      const ring = new Set(nodes); GS.E.forEach(([a, b], m) => { if (GS.kind[m] < 4 && (nodes.has(a) || nodes.has(b))) { ring.add(a); ring.add(b); } });
      const out = new Set(); GS.E.forEach(([a, b], m) => { if (!(ring.has(a) && ring.has(b))) out.add(m); }); return out;
    };
    const tolV = Math.max(0, o.simplifyTol ?? 0.02);
    const noX = o.noCrossing !== false;  // allowed volume increase for removing thin members (fraction of the re-optimised volume)
    const reroute = B => {                          // drop B and let members near B (plus the kept ones) carry its force, then round newcomers up
      const h = P.grid.h, loc = new Set(), X = GS.X;   // local region: joints within max(2h, member length) of each removed member's midpoint
      B.forEach(m => { const [a, b] = GS.E[m], c = X[a].map((v, t) => 0.5*(v + X[b][t])), rr = Math.max(2*h, GS.L[m]); X.forEach((p, j) => { if (Math.hypot(p[0]-c[0], p[1]-c[1], p[2]-c[2]) <= rr + 1e-9) loc.add(j); }); });
      const ban = new Set(), aLo = new Map();
      GS.E.forEach(([a, b], m) => { if (S.has(m) && !B.has(m)) aLo.set(m, Amin); else if (R.has(m) || B.has(m) || mates.has(m) || !(loc.has(a) && loc.has(b))) ban.add(m); });
      if (noX) GS.groups.forEach(g => { if (g.some(m => S.has(m) && !B.has(m))) g.forEach(m => { if (!S.has(m)) ban.add(m); }); });   // no new X crossings with kept members
      const tr0 = Date.now(); lps++; const r = solveLayout(highs, P, GS, BC, { banned: ban, aLo }); if (globalThis.GSTO_DEBUG) console.log('  reroute LP', Date.now() - tr0, 'ms, pool', GS.E.length - ban.size); if (!r.ok) return null;
      const keep = [...S].filter(m => !B.has(m)), fresh = GS.E.map((_, m) => m).filter(m => !ban.has(m) && !S.has(m) && r.a[m] > 1e-9*amaxOf(r.a)).sort((x, y) => r.a[y] - r.a[x]);
      const tryK = k => {                           // kept members + the k largest newcomers (crossing groups resolved by size)
        const S2 = new Set([...keep, ...fresh.slice(0, k).flatMap(orbit)]);
        if (noX) GS.groups.forEach(g => { const gm = g.filter(m => S2.has(m)); gm.sort((x, y) => r.a[y] - r.a[x]).slice(1).forEach(m => S2.delete(m)); });
        const r2 = lpOn(S2); return r2.ok ? { r: r2, S: S2 } : null;
      };
      let lo = 0, hi = fresh.length, best = tryK(hi); if (!best) return null;
      while (hi - lo > 1) { const mid = (lo + hi) >> 1, t = tryK(mid); if (t) { hi = mid; if (t.r.obj < best.r.obj) best = t; } else lo = mid; }   // fewest newcomers that still carry the load
      if (globalThis.GSTO_DEBUG) console.log('  reroute total', Date.now() - tr0, 'ms'); if (globalThis.GSTO_DEBUG) console.log('  reroute newcomers', fresh.length, '→', hi, 'r2', best.r.obj.toExponential(4), 'res', res.obj.toExponential(4));
      return best;
    };
    const accept = (B, t) => t && t.S.size <= S.size && (t.r.obj < res.obj*(1 - 1e-9) || (t.r.obj <= V1*(1 + tolV) && t.S.size < S.size));   // never more members; lighter, or fewer within the tolerance
    const use = m => deck.has(m) ? Infinity : Math.max(...res.q.map(qk => Math.abs(qk[m])))/(sMin*Amin);
    { const idle = new Set([...S].filter(m => use(m) < 1e-6).flatMap(orbit));   // members kept only by the bound: try all at once
      if (idle.size) { const t = reroute(idle); if (accept(idle, t)) { idle.forEach(m => R.add(m)); S = t.S; res = t.r; removed += idle.size; } } }
    for (let fails = 0; Date.now() < T && fails < (o.simplifyTries || 8); ) {   // then thin members one orbit at a time, least loaded first
      const thin = [...S].filter(m => !tried.has(m) && !deck.has(m) && res.a[m] <= Amin*(1 + 1e-6)).sort((a, b) => use(a) - use(b));
      if (!thin.length) break;
      const B = new Set(orbit(thin[0])), t = reroute(B);
      if (accept(B, t)) { B.forEach(m => R.add(m)); S = t.S; res = t.r; removed += B.size; fails = 0; }
      else { B.forEach(m => tried.add(m)); fails++; }
    }
    const lb = base.lpObj, gap = res.obj > 0 ? (res.obj - lb)/res.obj : 0;
    log(`최소 단면 ${Amin.toFixed(3)} mm²${o.minA > 0 ? ' (A_min)' : cat.length ? ' (카탈로그 최소)' : ` (${(100*(o.ratioC ?? 0.2)).toFixed(0)}% A_max)`}: 단순 올림 ${naive.toExponential(4)} → 재최적화 ${V1.toExponential(4)} → 가는 부재 ${removed}개 제거${tolV > 0 ? ` (체적 증가 허용 ${(100*tolV).toFixed(1)}%)` : ''} ${res.obj.toExponential(4)} (LP 하한 대비 +${(100*gap).toFixed(1)}%, LP ${lps}회) · 부재 ${n1}→${S.size}, 그중 A_min 부재 ${t1}→${thinCount(res)}`);
    return { ok: true, a: res.a, q: res.q, obj: res.obj, lpObj: lb, method: 'roundup', gap, thrC: Amin, thrT: Amin, amax0, naive, catalog: cat, iters: [{ kind: 'ROUNDUP', obj: res.obj }] };
  }

  // ---------- load capacity of the final structure: multiplier λ on each load case (Euler buckling without γ) ----------
  function loadCapacity(highs, P, GS, BC, S, fem) {
    const E = P.material.E, sT = P.material.sigmaT, sC = P.material.sigmaC, X = GS.X, mem = S.members, INF = highs.infinity;
    const capT = mem.map(e => sT*e.A), pcr = mem.map(e => Math.PI**2*E*sectionI(P, e.A)/(e.L*e.L)), capC = mem.map((e, i) => Math.min(sC*e.A, pcr[i]));
    const modeC = i => pcr[i] < sC*mem[i].A ? '좌굴' : '압축 항복';
    return BC.cases.map((cs, k) => {
      const out = { name: cs.name, elastic: null, plastic: null };
      if (fem && fem.results[k]) {                 // first member reaching its limit in the linear analysis
        let lam = Infinity, gov = -1, mode = '';
        fem.results[k].q.forEach((v, i) => { const c = v >= 0 ? capT[i] : capC[i]; if (Math.abs(v) > 1e-12*c && c/Math.abs(v) < lam) { lam = c/Math.abs(v); gov = i; mode = v >= 0 ? '인장 항복' : modeC(i); } });
        if (gov >= 0) Object.assign(out, { elastic: lam, elasticMember: gov, elasticMode: mode });
      }
      const rows = new Map(), R = [], C = [], V = [], rl = [];   // collapse: max λ s.t. B q = λ f, −capC ≤ q ≤ capT
      const row = i => { if (!rows.has(i)) { rows.set(i, rl.length); rl.push(0); } return rows.get(i); };
      mem.forEach((e, c) => { const d = sub(X[e.b], X[e.a]), l = norm(d);
        for (let t = 0; t < 3; t++) { const u = d[t]/l; if (Math.abs(u) < 1e-15) continue;
          if (BC.dofMap[3*e.a+t] >= 0) { R.push(row(3*e.a+t)); C.push(c); V.push(-u); }
          if (BC.dofMap[3*e.b+t] >= 0) { R.push(row(3*e.b+t)); C.push(c); V.push(u); } } });
      const cl = mem.length;
      for (let i = 0; i < cs.f.length; i++) if (cs.f[i] && BC.dofMap[i] >= 0) { R.push(row(i)); C.push(cl); V.push(-cs.f[i]); }
      const nr = rl.length, nc = cl + 1, cost = new Float64Array(nc); cost[cl] = -1;
      const model = highs.createModel({ numCols: nc, numRows: nr, colCost: cost, colLower: [...capC.map(v => -v), 0], colUpper: [...capT, INF], rowLower: rl, rowUpper: rl, matrix: cooToCsc(nr, nc, R, C, V) });
      try { model.options.set({ output_flag: false, time_limit: 30 }); model.run();
        if (model.getModelStatus() === highs.constants.modelStatus.optimal) {
          const x = model.getSolution().colValue, lam = x[cl], at = [];
          mem.forEach((e, i) => { if (x[i] >= capT[i]*(1 - 1e-6)) at.push([i, '인장 항복']); else if (-x[i] >= capC[i]*(1 - 1e-6)) at.push([i, modeC(i)]); });
          Object.assign(out, { plastic: lam, plasticAt: at });
        } else if (model.getModelStatus() === highs.constants.modelStatus.unbounded) out.plastic = Infinity;
      } catch {} finally { model.dispose(); }
      return out;
    });
  }

  // ---------- member adding (Gilbert & Tyas 2003): grow the ground structure from dual virtual strains ----------
  function growGroundStructure(highs, P, GS, BC, log = () => {}) {
    const t0 = Date.now(), o = P.options || {}, X = GS.X;
    const pool = candidatePool(P, GS, log), np = pool.a.length;
    const sT = P.material.sigmaT, sC = P.material.sigmaC, s = o.jointCost || 0;
    const tolV = o.addTol ?? 1e-3;                 // violation tolerance on the dual strain ratio
    const maxIt = o.addIter || 30;                 // member adding iterations
    log(`전 연결 후보: 장거리 부재 ${np}개 (겹침 제거 후), 시작 부재 ${GS.E.length}개, ${Date.now()-t0} ms`);
    if (!np) return { iters: 0, added: 0, pool: 0 };
    const inGS = new Uint8Array(np), nodeMir = GS.sym ? GS.sym.nodeMir : [];
    let added = 0, it = 0, last = null, last0 = null, flat = 0;
    for (; it < maxIt; it++) {
      const r = solveLayout(highs, P, GS, BC, { central: true }); last = r;   // interior (central) duals: fewer spurious violations than a vertex
      if (!r.ok || !r.U) { log(`부재 추가 반복 ${it+1}: LP ${r.ok ? '쌍대해 없음' : '실패'} → 중단`); break; }
      const ratio = (a, b) => {                    // Σ_k σ·|ε_k| / (L+s) with ε_k = d·(u_b − u_a)
        const d0 = X[b][0]-X[a][0], d1 = X[b][1]-X[a][1], d2 = X[b][2]-X[a][2], L = Math.hypot(d0, d1, d2); let v = 0;
        for (const u of r.U) { const w = (d0*(u[3*b]-u[3*a]) + d1*(u[3*b+1]-u[3*a+1]) + d2*(u[3*b+2]-u[3*a+2]))/L; v += w > 0 ? sT*w : -sC*w; }
        return v/(L + s);
      };
      const viol = [];
      for (let i = 0; i < np; i++) { if (inGS[i]) continue; const v = ratio(pool.a[i], pool.b[i]); if (v > 1 + tolV) { if (GS.memberId(pool.a[i], pool.b[i]) >= 0) inGS[i] = 1; else viol.push([i, v]); } }
      if (!viol.length) { log(`부재 추가 반복 ${it+1}: 위반 후보 없음 → 전 연결 최적 (체적 ${r.obj.toExponential(4)}, ${r.dt} ms)`); break; }
      viol.sort((x, y) => y[1] - x[1]);
      if (last0 !== null && r.obj > last0*(1 - 1e-5)) { if (++flat >= 3) { log(`부재 추가 반복 ${it+1}: 체적 변화가 3회 연속 0.001% 미만 → 수렴으로 판단 (남은 위반 ${viol.length}개, 최대 비 ${viol[0][1].toFixed(3)})`); break; } } else flat = 0;
      last0 = r.obj;
      const cap = Math.max(100, Math.ceil(0.1*GS.E.length)), pairs = [], seen = new Set();
      for (const [i] of viol) {
        if (pairs.length >= cap) break;
        const orb = [[pool.a[i], pool.b[i]]];      // add with all mirror images
        for (let pass = 0; pass < 3; pass++) orb.slice().forEach(([a, b]) => nodeMir.forEach(nm => { if (nm[a] >= 0 && nm[b] >= 0) orb.push([nm[a], nm[b]]); }));
        orb.forEach(([a, b]) => { const k = Math.min(a, b)+','+Math.max(a, b); if (!seen.has(k) && GS.memberId(a, b) < 0) { seen.add(k); pairs.push([a, b]); } });
        inGS[i] = 1;
      }
      const ids = GS.addMembers(pairs); added += ids.length;
      if (GS.sym) GS.sym = symmetryMaps({ ...P, options: { ...P.options, symmetry: [0, 1, 2].map(t => GS.sym.planes.includes(t) ? 1 : 0), symExact: !!GS.sym.exact } }, GS);   // new members join orbits
      log(`부재 추가 반복 ${it+1}: 위반 ${viol.length}개 (최대 비 ${viol[0][1].toFixed(3)}) 중 ${ids.length}개 추가 → 부재 ${GS.E.length}, 체적 ${r.obj.toExponential(4)}, ${r.dt} ms`);
    }
    return { iters: it + 1, added, pool: np, obj: last && last.obj };
  }

  function subSym(GS, pairs) {                     // mirror orbits for a member list on the nodes of GS
    if (!GS.sym) return null;
    const key = new Map(); pairs.forEach(([a, b], i) => key.set(Math.min(a, b)+','+Math.max(a, b), i));
    const memMir = GS.sym.nodeMir.map(nm => pairs.map(([a, b]) => { const k = key.get(Math.min(nm[a], nm[b])+','+Math.max(nm[a], nm[b])); return nm[a] < 0 || nm[b] < 0 || k === undefined ? -1 : k; }));
    const par = Int32Array.from(pairs.keys()), find = x => { while (par[x] !== x) x = par[x] = par[par[x]]; return x; };
    memMir.forEach(mm => mm.forEach((o, m) => { if (o >= 0) { const a = find(m), b = find(o); if (a !== b) par[Math.max(a, b)] = Math.min(a, b); } }));
    const orb = new Map(); pairs.forEach((_, m) => { const r = find(m); if (!orb.has(r)) orb.set(r, []); orb.get(r).push(m); });
    const orbits = [...orb.values()].filter(o => o.length > 1), orbitOf = new Int32Array(pairs.length).fill(-1);
    orbits.forEach((o, oi) => o.forEach(m => { orbitOf[m] = oi; }));
    const missing = memMir.reduce((acc, mm) => acc + mm.filter(v => v < 0).length, 0), out = { planes: GS.sym.planes, nodeMir: GS.sym.nodeMir, memMir, orbits, orbitOf, missing };
    if (GS.sym.exact && !missing) addReduction(out, { X: GS.X, E: pairs });
    return out;
  }

  // ---------- geometry optimisation (He & Gilbert 2015): SLP on free node positions with move limits ----------
  function geometryOptimize(highs, P, GS, BC, opt, log = () => {}) {
    const o = P.options || {}, h = P.grid.h, X = GS.X, tol = 1e-6*h, obs = P.obstacles || [];
    const R = (o.geomRange ?? 0.5)*h;              // cumulative move limit per axis [mm]
    let mu = (o.geomStep ?? 0.2)*h;                // initial per-iteration move limit [mm]
    const muMin = 0.005*h, maxIt = o.geomIter || 30, s = o.jointCost || 0, alpha = o.areaCut ?? 1e-3;
    const sT = P.material.sigmaT, sC = P.material.sigmaC, K = BC.cases.length, INF = highs.infinity;
    const on0 = presence(opt, alpha);
    let nodesAll = []; { const st = new Set(); GS.E.forEach(([a, b], m) => { if (on0(m)) { st.add(a); st.add(b); } }); nodesAll = [...st]; }
    const fixedNode = j => BC.supNodes.has(j) || BC.loadedAny[j] || GS.userNodes.has(j);
    const axes = [0, 1, 2].filter(t => t !== BC.flat);
    const dom = [0, 1, 2].map(t => [P.grid.origin[t], P.grid.origin[t] + [P.grid.nx, P.grid.ny, P.grid.nz][t]*h]);
    let Y = X.map(p => p.slice());
    const complete = Yc => {                       // all straight connections among the structure's joints (none passing through a joint)
      const out = [];
      for (let i = 0; i < nodesAll.length; i++) for (let k = i + 1; k < nodesAll.length; k++) {
        const a = nodesAll[i], b = nodesAll[k];
        if (o.maxLen > 0 && norm(sub(Yc[b], Yc[a])) > o.maxLen*h + tol) continue;
        if (obs.some(bx => segHitsBoxInterior(Yc[a], Yc[b], bx, tol))) continue;
        let hit = false; for (const c of nodesAll) if (c !== a && c !== b && pointInSegInterior(Yc[c], Yc[a], Yc[b], 1e-4*h)) { hit = true; break; }
        if (!hit) out.push([Math.min(a, b), Math.max(a, b)]);
      }
      return out;
    };
    const subGS = (Yc, prs) => {                   // temporary ground structure on moved joints
      const L = new Float64Array(prs.length), D = new Float64Array(3*prs.length), key = new Map();
      prs.forEach(([a, b], i) => { const d = sub(Yc[b], Yc[a]), l = norm(d); L[i] = l; for (let t = 0; t < 3; t++) D[3*i+t] = d[t]/l; key.set(a+','+b, i); });
      const g = { X: Yc, E: prs, L, D, kind: prs.map(() => 4), groups: [], memberId: (a, b) => { const k = key.get(Math.min(a, b)+','+Math.max(a, b)); return k === undefined ? -1 : k; } };
      g.sym = subSym(GS, prs); return g;
    };
    const solveAt = Yc => { const prs = complete(Yc), g = subGS(Yc, prs), r = solveLayout(highs, P, g, BC); r.pairs = prs; return r; };   // LP over all straight connections of the joints at Yc
    const tg0 = Date.now(); let cur = solveAt(Y); if (globalThis.GSTO_DEBUG) console.log('geo init', Date.now() - tg0, 'ms pairs', cur.pairs.length, 'nodes', nodesAll.length);
    if (!cur.ok) { log('기하 최적화: 초기 LP 실패 → 생략'); return null; }
    const V0 = opt.obj, Vc0 = cur.obj; let acc = 0, it = 0, small = 0;
    const okGeom = (Yc, j, prs) => !obs.some(bx => inBoxOpen(Yc[j], bx, tol)) && prs.every(([a, b]) => (a !== j && b !== j) || !obs.some(bx => segHitsBoxInterior(Yc[a], Yc[b], bx, tol)));
    for (; it < maxIt && mu > muMin; it++) {
      const onC = presence(cur, alpha), act = []; cur.pairs.forEach((p, i) => { if (onC(i)) act.push(i); });
      const used = new Set(act.flatMap(i => cur.pairs[i])), mov = nodesAll.filter(j => used.has(j) && !fixedNode(j)), movIdx = new Map(mov.map((j, i) => [j, i]));
      if (!mov.length) break;
      const cost = [], lo = [], up = [], Rr = [], Cc = [], Vv = [], rL = [], rU = []; let nRow = 0;
      const addCol = (c, l, u) => { cost.push(c); lo.push(l); up.push(u); return cost.length - 1; };
      const addRow = (l, u) => { rL.push(l); rU.push(u); return nRow++; };
      const put = (r, c, v) => { if (v !== 0) { Rr.push(r); Cc.push(c); Vv.push(v); } };
      const dcol = mov.map(j => axes.map(t => addCol(0, Math.max(-mu, dom[t][0] - Y[j][t], X[j][t] - R - Y[j][t]), Math.min(mu, dom[t][1] - Y[j][t], X[j][t] + R - Y[j][t]))));
      const dc = (j, t) => { const i = movIdx.get(j); if (i === undefined) return -1; const ai = axes.indexOf(t); return ai < 0 ? -1 : dcol[i][ai]; };
      const eq = []; for (let k = 0; k < K; k++) { const rows = new Map(); used.forEach(j => { for (let t = 0; t < 3; t++) if (BC.dofMap[3*j+t] >= 0) rows.set(3*j+t, addRow(BC.cases[k].f[3*j+t], BC.cases[k].f[3*j+t])); }); eq.push(rows); }
      act.forEach(i => {
        const [a, b] = cur.pairs[i], d = sub(Y[b], Y[a]), L = norm(d), u = d.map(v => v/L), a0 = cur.a[i];
        const J = [0, 1, 2].map(r => [0, 1, 2].map(c => ((r === c ? 1 : 0) - u[r]*u[c])/L));   // ∂u/∂(x_b − x_a)
        const ca = addCol(L + s, 0, INF);
        [0, 1, 2].forEach(t => { const cb = dc(b, t), cA = dc(a, t); if (cb >= 0) cost[cb] += a0*u[t]; if (cA >= 0) cost[cA] -= a0*u[t]; });   // dV ≈ a0·u·(δb − δa)
        for (let k = 0; k < K; k++) {
          const cq = addCol(0, -INF, INF), q0 = cur.q[k][i];
          for (let t = 0; t < 3; t++) {
            const ra = eq[k].get(3*a+t), rb = eq[k].get(3*b+t);
            if (ra !== undefined) put(ra, cq, -u[t]);
            if (rb !== undefined) put(rb, cq, u[t]);
            if (q0) for (let c = 0; c < 3; c++) { const v = q0*J[t][c], cb = dc(b, c), cA = dc(a, c);
              if (rb !== undefined) { if (cb >= 0) put(rb, cb, v); if (cA >= 0) put(rb, cA, -v); }
              if (ra !== undefined) { if (cb >= 0) put(ra, cb, -v); if (cA >= 0) put(ra, cA, v); } }
          }
          let r = addRow(-INF, 0); put(r, cq, 1); put(r, ca, -sT);
          r = addRow(-INF, 0); put(r, cq, -1); put(r, ca, -sC);
        }
      });
      const model = highs.createModel({ numCols: cost.length, numRows: nRow, colCost: cost, colLower: lo, colUpper: up, rowLower: rL, rowUpper: rU, matrix: cooToCsc(nRow, cost.length, Rr, Cc, Vv) });
      let x = null;
      try { model.options.set({ output_flag: false, solver: 'ipm', time_limit: o.timeLimit || 60 }); model.run();
        if (model.getModelStatus() === highs.constants.modelStatus.optimal) x = model.getSolution().colValue; } catch { x = null; } finally { model.dispose(); }
      if (!x) { mu *= 0.5; continue; }
      let dl = mov.map((j, i) => [0, 1, 2].map(t => { const ai = axes.indexOf(t); return ai < 0 ? 0 : x[dcol[i][ai]]; }));
      if (GS.sym) GS.sym.planes.forEach((t, p) => {   // keep the step mirror-symmetric
        const nm = GS.sym.nodeMir[p];
        dl = dl.map((v, i) => { const im = movIdx.get(nm[mov[i]]); if (im === undefined) return v; const w = dl[im]; return [0, 1, 2].map(c => 0.5*(v[c] + (c === t ? -w[c] : w[c]))); });
      });
      const Y2 = Y.map(p => p.slice()); mov.forEach((j, i) => { for (let t = 0; t < 3; t++) Y2[j][t] += dl[i][t]; });
      const actPairs = act.map(i => cur.pairs[i]);
      mov.forEach(j => { if (!okGeom(Y2, j, actPairs)) Y2[j] = Y[j].slice(); });
      const r2 = solveAt(Y2);
      if (globalThis.GSTO_DEBUG) console.log('slp', (Date.now()/1000 % 1000).toFixed(1), it, 'mu', mu.toFixed(3), 'cur', cur.obj, 'try', r2.ok && r2.obj, r2.status);
      if (r2.ok && r2.obj < cur.obj*(1 - 1e-7)) { small = r2.obj > cur.obj*(1 - 1e-4) ? small + 1 : 0; Y = Y2; cur = r2; acc++; if (small >= 3) break; } else mu *= 0.5;   // stop after 3 accepted steps gaining < 0.01 %
    }
    if (globalThis.GSTO_DEBUG) console.log('geo loop end', Date.now() - tg0);
    // merge joints that came together
    const mergeTol = (o.geomMerge ?? 0.1)*h; let merged = 0;
    for (let changed = true; changed;) {
      changed = false;
      const onC = presence(cur, alpha), usedN = [...new Set(cur.pairs.filter((_, i) => onC(i)).flat())];
      outer: for (let i = 0; i < usedN.length; i++) for (let k = i + 1; k < usedN.length; k++) {
        const a = usedN[i], b = usedN[k]; if (norm(sub(Y[a], Y[b])) > mergeTol) continue;
        const fa = fixedNode(a), fb = fixedNode(b); if (fa && fb) continue;
        const keep = fb ? b : a, drop = keep === a ? b : a;
        if (!fa && !fb) Y[keep] = Y[a].map((v, t) => 0.5*(v + Y[b][t]));
        nodesAll = nodesAll.filter(j => j !== drop); merged++; changed = true; break outer;
      }
      if (changed) { const r3 = solveAt(Y); if (!r3.ok) break; cur = r3; }
    }
    // final sizing on the moved joints (minimum-area rules re-applied when active)
    let fin = cur;
    if (opt.method === 'roundup') { const g = subGS(Y, cur.pairs), r = solveRoundUp(highs, { ...P, options: { ...P.options, timeLimit: Math.min(P.options?.timeLimit || 60, 20), minA: opt.thrC } }, g, BC); if (r.ok) fin = { ...r, pairs: cur.pairs }; }
    if (opt.method === 'minA') { const g = subGS(Y, cur.pairs), r = solveMinArea(highs, { ...P, options: { ...P.options, minAOnly: 'support', timeLimit: Math.min(P.options?.timeLimit || 60, 20) } }, g, BC, globalThis.GSTO_DEBUG ? (m => console.log('  [geo minA] ' + m)) : () => {}); if (r.ok) { fin = { ...r, pairs: cur.pairs }; } else if (globalThis.GSTO_DEBUG) console.log('geo minA fail', r.msg); }
    if (globalThis.GSTO_DEBUG) console.log('geo minA end', Date.now() - tg0, fin === cur ? 'minA not applied' : 'minA ok');
    if (fin.obj >= V0*(1 - 1e-9)) { log(`기하 최적화: 체적이 줄지 않아 (${V0.toExponential(4)} → ${fin.obj.toExponential(4)}) 격자 결과를 유지합니다`); return null; }
    const onF = presence(fin, alpha), keepI = fin.pairs.map((_, i) => i).filter(i => onF(i));
    const nodeOf = new Map(); GS.dead = GS.dead || new Set();
    [...new Set(keepI.flatMap(i => fin.pairs[i]))].forEach(j => { if (norm(sub(Y[j], X[j])) > 1e-6*h) { nodeOf.set(j, addNode(GS, BC, Y[j])); GS.dead.add(j); } else nodeOf.set(j, j); });
    const prsF = keepI.map(i => fin.pairs[i].map(j => nodeOf.get(j)));
    GS.addMembers(prsF);
    const b = GS.E.length, A = new Float64Array(b), Q = Array.from({ length: K }, () => new Float64Array(b));
    keepI.forEach((i, n) => { const m = GS.memberId(prsF[n][0], prsF[n][1]); A[m] = fin.a[i]; for (let k = 0; k < K; k++) Q[k][m] = fin.q[k][i]; });
    const nMoved = [...nodeOf].filter(([j, k]) => j !== k).length;
    log(`기하 최적화: SLP ${it}회 (채택 ${acc}), 이동 노드 ${nMoved}, 노드 병합 ${merged}, 체적 ${V0.toExponential(4)} → ${fin.obj.toExponential(4)} (${(100*(fin.obj - V0)/V0).toFixed(2)}%; 같은 노드 전 연결 ${Vc0.toExponential(4)})`);
    return { ...opt, ...(fin.thrC !== undefined ? { thrC: fin.thrC, thrT: fin.thrT, amax0: fin.amax0 } : {}), a: A, q: Q, obj: fin.obj, geom: { V0, V: fin.obj, iters: it, accepted: acc, merged } };
  }

  // ---------- crossings of long members: joint at the crossing point (also T-contacts on member interiors) ----------
  function splitCrossings(P, GS, BC, S) {
    const X = GS.X, tol = 1e-6*P.grid.h, mem = S.members, cuts = mem.map(() => []);
    const hit = (p, q, c, d) => {                  // interior intersection point of segments pq and cd (or null)
      const u = sub(q, p), v = sub(d, c), w = sub(p, c), a = dot(u, u), b = dot(u, v), cc = dot(v, v), dd = dot(u, w), e = dot(v, w), den = a*cc - b*b;
      if (den < 1e-12*a*cc) return null;
      const sP = (b*e - cc*dd)/den, tP = (a*e - b*dd)/den, Lu = Math.sqrt(a), Lv = Math.sqrt(cc);
      if (sP < -tol/Lu || sP > 1 + tol/Lu || tP < -tol/Lv || tP > 1 + tol/Lv) return null;
      const P1 = p.map((x, i) => x + sP*u[i]), P2 = c.map((x, i) => x + tP*v[i]);
      if (norm(sub(P1, P2)) > 10*tol) return null;
      return { s: sP, t: tP, p: P1.map((x, i) => 0.5*(x + P2[i])), inS: sP > tol/Lu && sP < 1 - tol/Lu, inT: tP > tol/Lv && tP < 1 - tol/Lv };
    };
    let nCross = 0;
    for (let i = 0; i < mem.length; i++) for (let j = i + 1; j < mem.length; j++) {
      const A = mem[i], B = mem[j];
      if (!A.long && !B.long) continue;              // one-cell members are handled by the crossing groups
      if (A.a === B.a || A.a === B.b || A.b === B.a || A.b === B.b) continue;
      const r = hit(X[A.a], X[A.b], X[B.a], X[B.b]); if (!r || (!r.inS && !r.inT)) continue;
      if (r.inS && r.inT) { const c = addNode(GS, BC, r.p); cuts[i].push(c); cuts[j].push(c); nCross++; }
      else if (r.inS) cuts[i].push(r.t < 0.5 ? B.a : B.b);   // endpoint of B lies on A's interior
      else cuts[j].push(r.s < 0.5 ? A.a : A.b);
    }
    if (!cuts.some(c => c.length)) return 0;
    const out = [];
    mem.forEach((e, i) => {
      if (!cuts[i].length) { out.push(e); return; }
      const u = sub(X[e.b], X[e.a]).map(v => v/e.L), t = j => dot(sub(X[j], X[e.a]), u);
      const seq = [e.a, ...[...new Set(cuts[i])].sort((p, q) => t(p) - t(q)), e.b];
      for (let k = 0; k + 1 < seq.length; k++) { const a = seq[k], b = seq[k+1]; if (a === b) continue;
        out.push({ ...e, a: Math.min(a, b), b: Math.max(a, b), L: norm(sub(X[b], X[a])), crossSplit: true }); }
    });
    S.members = out; S.volume = out.reduce((v, e) => v + e.L*e.A, 0);
    return nCross;
  }

  // ---------- post-processing ----------
  function extractStructure(P, GS, BC, opt) {
    const alpha = P.options?.areaCut ?? 1e-3, K = BC.cases.length, on = presence(opt, alpha);
    let mem = [];
    GS.E.forEach(([a,b], m) => { if (on(m)) mem.push({ a, b, A: opt.a[m], q: opt.q.map(qk => qk[m]), gs: [m], long: GS.kind[m] === 4 }); });
    let merged = 0;
    if (P.options?.mergeCollinear !== false) {
      const X = GS.X, isFixed = j => BC.supNodes.has(j);   // support nodes keep their joints
      let changed = true;
      while (changed) {
        changed = false;
        const inc = new Map();
        mem.forEach((e, i) => { [e.a, e.b].forEach(j => { if (!inc.has(j)) inc.set(j, []); inc.get(j).push(i); }); });
        for (const [j, list] of inc) {
          if (list.length !== 2 || BC.loadedAny[j] || isFixed(j)) continue;
          const [e1, e2] = list.map(i => mem[i]), o1 = e1.a === j ? e1.b : e1.a, o2 = e2.a === j ? e2.b : e2.a;
          const u = sub(X[o1], X[j]), v = sub(X[o2], X[j]);
          if (dot(u,v)/(norm(u)*norm(v)) > -1+1e-9) continue;
          mem = mem.filter((_, i) => i !== list[0] && i !== list[1]);
          mem.push({ a: Math.min(o1,o2), b: Math.max(o1,o2), A: Math.max(e1.A, e2.A), q: e1.q.map((x, k) => 0.5*(x + e2.q[k])), gs: e1.gs.concat(e2.gs), mergedFrom: j, long: e1.long || e2.long });
          merged++; changed = true; break;
        }
      }
    }
    let V = 0; mem.forEach(e => { e.L = norm(sub(GS.X[e.b], GS.X[e.a])); V += e.L*e.A; });
    if (opt.thrT !== undefined) {                   // tension rods: tension (or zero) in every load case, thinner than the general limit
      const qm = Math.max(1, ...mem.flatMap(e => e.q.map(Math.abs)));
      mem.forEach(e => { e.rod = e.q.every(v => v >= -1e-9*qm) && e.q.some(v => v > 1e-6*qm) && e.A < opt.thrC*(1 - 1e-6); });
    }
    return { members: mem, volume: V, merged, K };
  }

  // ---------- symmetric eigen (Householder tridiagonalisation + implicit QL, JAMA) ----------
  function symEig(Ain, n) {
    const V = Ain, d = new Float64Array(n), e = new Float64Array(n);  // V: Float64Array n*n row-major, overwritten by eigenvectors (columns)
    const I = (i,j) => i*n+j;
    for (let j = 0; j < n; j++) d[j] = V[I(n-1,j)];
    for (let i = n-1; i > 0; i--) {
      let scale = 0, hh = 0;
      for (let k = 0; k < i; k++) scale += Math.abs(d[k]);
      if (scale === 0) { e[i] = d[i-1]; for (let j = 0; j < i; j++) { d[j] = V[I(i-1,j)]; V[I(i,j)] = 0; V[I(j,i)] = 0; } }
      else {
        for (let k = 0; k < i; k++) { d[k] /= scale; hh += d[k]*d[k]; }
        let f = d[i-1], g = Math.sqrt(hh); if (f > 0) g = -g;
        e[i] = scale*g; hh -= f*g; d[i-1] = f-g;
        for (let j = 0; j < i; j++) e[j] = 0;
        for (let j = 0; j < i; j++) {
          f = d[j]; V[I(j,i)] = f; g = e[j] + V[I(j,j)]*f;
          for (let k = j+1; k <= i-1; k++) { g += V[I(k,j)]*d[k]; e[k] += V[I(k,j)]*f; }
          e[j] = g;
        }
        f = 0; for (let j = 0; j < i; j++) { e[j] /= hh; f += e[j]*d[j]; }
        const hh2 = f/(hh+hh); for (let j = 0; j < i; j++) e[j] -= hh2*d[j];
        for (let j = 0; j < i; j++) { f = d[j]; g = e[j]; for (let k = j; k <= i-1; k++) V[I(k,j)] -= (f*e[k] + g*d[k]); d[j] = V[I(i-1,j)]; V[I(i,j)] = 0; }
      }
      d[i] = hh;
    }
    for (let i = 0; i < n-1; i++) {
      V[I(n-1,i)] = V[I(i,i)]; V[I(i,i)] = 1; const hh = d[i+1];
      if (hh !== 0) {
        for (let k = 0; k <= i; k++) d[k] = V[I(k,i+1)]/hh;
        for (let j = 0; j <= i; j++) { let g = 0; for (let k = 0; k <= i; k++) g += V[I(k,i+1)]*V[I(k,j)]; for (let k = 0; k <= i; k++) V[I(k,j)] -= g*d[k]; }
      }
      for (let k = 0; k <= i; k++) V[I(k,i+1)] = 0;
    }
    for (let j = 0; j < n; j++) { d[j] = V[I(n-1,j)]; V[I(n-1,j)] = 0; }
    V[I(n-1,n-1)] = 1; e[0] = 0;
    for (let i = 1; i < n; i++) e[i-1] = e[i];
    e[n-1] = 0;
    let f = 0, tst1 = 0; const eps = Math.pow(2,-52);
    for (let l = 0; l < n; l++) {
      tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
      let m = l; while (m < n) { if (Math.abs(e[m]) <= eps*tst1) break; m++; }
      if (m > l) {
        do {
          let g = d[l], p = (d[l+1]-g)/(2*e[l]), r = Math.hypot(p,1); if (p < 0) r = -r;
          d[l] = e[l]/(p+r); d[l+1] = e[l]*(p+r); const dl1 = d[l+1]; let hh = g-d[l];
          for (let i = l+2; i < n; i++) d[i] -= hh;
          f += hh; p = d[m]; let c = 1, c2 = c, c3 = c; const el1 = e[l+1]; let s = 0, s2 = 0;
          for (let i = m-1; i >= l; i--) {
            c3 = c2; c2 = c; s2 = s; g = c*e[i]; hh = c*p; r = Math.hypot(p, e[i]); e[i+1] = s*r; s = e[i]/r; c = p/r; p = c*d[i] - s*g; d[i+1] = hh + s*(c*g + s*d[i]);
            for (let k = 0; k < n; k++) { hh = V[I(k,i+1)]; V[I(k,i+1)] = s*V[I(k,i)] + c*hh; V[I(k,i)] = c*V[I(k,i)] - s*hh; }
          }
          p = -s*s2*c3*el1*e[l]/dl1; e[l] = s*p; d[l] = c*p;
        } while (Math.abs(e[l]) > eps*tst1);
      }
      d[l] += f; e[l] = 0;
    }
    return { values: d, vectors: V };              // vectors: column j = eigenvector j
  }

  // ---------- kinematic (DOF) check: SVD of equilibrium matrix via A·Aᵀ ----------
  function stabilityCheck(GS, BC, S, maxDof = 3000) {
    const X = GS.X, mem = S.members, nodes = [...new Set(mem.flatMap(e => [e.a, e.b]))].sort((a,b) => a-b);
    const loc = new Map(); let nf = 0; const dofOf = [];
    nodes.forEach(j => { const r = []; for (let t = 0; t < 3; t++) r.push(BC.fixed[3*j+t] ? -1 : nf++); loc.set(j, r); dofOf.push(r); });
    const b = mem.length, out = { nNodes: nodes.length, nf, b };
    // local rank check per node
    const local = [];
    nodes.forEach(j => {
      const M = new Float64Array(9); let cnt = 0;
      mem.forEach(e => { if (e.a !== j && e.b !== j) return; const d = sub(X[e.b], X[e.a]), l = norm(d); for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) M[3*r+c] += d[r]*d[c]/(l*l); cnt++; });
      for (let t = 0; t < 3; t++) if (BC.fixed[3*j+t]) M[4*t] += 1;
      const ev = symEig(M, 3).values, rk = [...ev].filter(v => v > 1e-8).length;
      if (rk < 3) local.push({ node: j, rank: rk });
    });
    out.localDeficient = local;
    if (nf === 0) return { ...out, r: 0, m: 0, s: b, mechNodes: [], skipped: false };
    if (nf > maxDof) return { ...out, skipped: true, msg: `자유도 ${nf} > ${maxDof}: 전역 SVD 생략` };
    const Mx = new Float64Array(nf*nf);
    mem.forEach(e => {
      const d = sub(X[e.b], X[e.a]), l = norm(d), u = [d[0]/l, d[1]/l, d[2]/l], ia = loc.get(e.a), ib = loc.get(e.b);
      const idx = [ia[0],ia[1],ia[2],ib[0],ib[1],ib[2]], val = [-u[0],-u[1],-u[2],u[0],u[1],u[2]];
      for (let p = 0; p < 6; p++) { if (idx[p] < 0) continue; for (let q = 0; q < 6; q++) if (idx[q] >= 0) Mx[idx[p]*nf+idx[q]] += val[p]*val[q]; }
    });
    const { values, vectors } = symEig(Mx, nf), lmax = Math.max(...values), thr = 1e-10*lmax;
    const nullIdx = []; values.forEach((v, i) => { if (v <= thr) nullIdx.push(i); });
    const r = nf - nullIdx.length, m = nullIdx.length, s = b - r;
    const phi = new Float64Array(nodes.length);
    nodes.forEach((j, ni) => { let acc = 0; dofOf[ni].forEach(dd => { if (dd < 0) return; nullIdx.forEach(c => { const v = vectors[dd*nf+c]; acc += v*v; }); }); phi[ni] = Math.sqrt(acc); });
    const mechNodes = []; nodes.forEach((j, ni) => { if (phi[ni] > 0.1) mechNodes.push({ node: j, phi: phi[ni] }); });
    const modes = nullIdx.slice(0, 300).map(c => {   // mechanism modes in global DOF numbering
      const v = new Float64Array(3*X.length);
      nodes.forEach((j, ni) => dofOf[ni].forEach((dd, t) => { if (dd >= 0) v[3*j+t] = vectors[dd*nf+c]; }));
      return v;
    });
    return { ...out, r, m, s, mechNodes, modes, lambdaMax: lmax, skipped: false };
  }


  // ---------- mechanism removal: greedy stabilising members from the ground structure (method B) ----------
  function stabilizeGreedy(P, GS, BC, S, stab) {
    const vnorm = w => { let q = 0; for (let i = 0; i < w.length; i++) q += w[i]*w[i]; return Math.sqrt(q); };   // no argument spreading
    const X = GS.X, ratio = P.options?.stabArea ?? 0.01;  // stabilising member area / max area
    const h = P.grid.h, Abase = S.aStruct ?? ratio*Math.max(...S.members.map(e => e.A));   // structural minimum for added members
    const Pc = new Map();                          // largest compression force meeting each node
    S.members.forEach(e => { const c = Math.max(0, ...e.q.map(v => -v)); [e.a, e.b].forEach(j => Pc.set(j, Math.max(Pc.get(j) || 0, c))); });
    const brace = P.options?.braceRatio ?? 0.02;   // notional lateral force / compression force
    const sizeOf = (a, b) => Math.max(Abase, brace*Math.max(Pc.get(a) || 0, Pc.get(b) || 0)/P.material.sigmaC);
    const active = new Set(S.members.flatMap(e => [e.a, e.b])), present = new Set(S.members.flatMap(e => e.gs || []));
    for (let j = 0; j < X.length; j++) if (BC.fixed[3*j] && BC.fixed[3*j+1] && BC.fixed[3*j+2]) active.add(j);   // fully fixed support nodes are free anchors
    const groupsOf = new Map(); GS.groups.forEach((g, gi) => g.forEach(m => { if (!groupsOf.has(m)) groupsOf.set(m, []); groupsOf.get(m).push(gi); }));
    const blocked = m => (groupsOf.get(m) || []).some(gi => GS.groups[gi].some(o => o !== m && present.has(o)));
    const used = new Set(); present.forEach(m => GS.E[m].forEach(j => used.add(j)));  // grid points occupied by members
    const dead = GS.dead || new Set();             // grid joints vacated by geometry optimisation
    const cand = [], okPart = m => m >= 0 && GS.kind[m] < 4 && !present.has(m) && !blocked(m) && !dead.has(GS.E[m][0]) && !dead.has(GS.E[m][1]);
    const incM = new Map(); GS.E.forEach(([a, b], m) => { if (GS.kind[m] < 4) [a, b].forEach(j => { if (!incM.has(j)) incM.set(j, []); incM.get(j).push(m); }); });
    const Rf = (P.options?.braceRadius ?? 6)*h, Rn = 1.75*h, obsF = P.obstacles || [], tolF = 1e-6*h, offGrid = j => j >= GS.nGrid;
    const fKey = c => c.a+','+c.b, freeUsed = new Set();
    const freeCands = (c, Rmax = Rf) => {         // straight members (not one-cell grid members) to active joints within Rmax
      const list = [];
      active.forEach(d => {
        if (d === c) return;
        const m = GS.memberId(c, d); if (m >= 0 && GS.kind[m] < 4) return;   // one-cell members are grid candidates
        const L = norm(sub(X[d], X[c])); if (L > Rmax || L < 1e-9) return;
        const a = Math.min(c, d), b = Math.max(c, d);
        if (obsF.some(bx => segHitsBoxInterior(X[a], X[b], bx, tolF))) return;
        const lo = [0, 1, 2].map(t => Math.min(X[a][t], X[b][t]) - tolF), hi = [0, 1, 2].map(t => Math.max(X[a][t], X[b][t]) + tolF);
        for (const j of active) { const p = X[j]; if (j === a || j === b || p[0] < lo[0] || p[0] > hi[0] || p[1] < lo[1] || p[1] > hi[1] || p[2] < lo[2] || p[2] > hi[2]) continue; if (pointInSegInterior(p, X[a], X[b], tolF)) return; }
        const dd = sub(X[b], X[a]); list.push({ a, b, parts: [], free: true, L, d: dd.map(v => v/L) });
      });
      return list;
    };
    const crossesAny = c => [...S.members, ...out].some(e => (c.free || e.long) && e.a !== c.a && e.a !== c.b && e.b !== c.a && e.b !== c.b && segsCross(X[c.a], X[c.b], X[e.a], X[e.b], tolF));
    GS.E.forEach(([a, b], m) => { if (active.has(a) && active.has(b) && okPart(m)) cand.push({ a, b, parts: [m] }); });
    const midOf = new Map();                       // interior grid point of a merged member → member index
    S.members.forEach((e, i) => { if (e.gs && e.gs.length > 1) e.gs.forEach(m => GS.E[m].forEach(j => { if (j !== e.a && j !== e.b) midOf.set(j, i); })); });
    const splits = [];                              // [node, member index]: interior points turned into joints (lateral bracing)
    const halfOf = c => { const e = S.members[midOf.get(c)], d = sub(X[e.a], X[c]), L = norm(d); return { a: c, b: e.a, L, d: [-d[0]/L, -d[1]/L, -d[2]/L] }; };
    active.forEach(a => {                          // straight chains along grid lines (merged one-cell members)
      if (a >= GS.nGrid) return;
      const [i0, j0, k0] = GS.ijk[a];
      for (const o of OFFS) for (const sg of [1, -1]) {
        const parts = [], mids = []; let prev = a;
        for (let st = 1; st < 64; st++) {
          const nb = GS.gidx(i0+sg*o[0]*st, j0+sg*o[1]*st, k0+sg*o[2]*st); if (nb < 0) break;
          const m = GS.memberId(prev, nb); if (!okPart(m)) break;
          parts.push(m); prev = nb;
          if (active.has(nb)) { if (st > 1 && a < nb) cand.push({ a, b: nb, parts: parts.slice(), mids: mids.slice() }); break; }
          if (used.has(nb)) break;
          mids.push(nb);
        }
      }
    });
    cand.forEach(c => { const d = sub(X[c.b], X[c.a]); c.L = norm(d); c.d = d.map(x => x/c.L); });
    let N = stab.modes.map(v => Float64Array.from(v)), out = [];
    const anchor = j => BC.fixed[3*j] && BC.fixed[3*j+1] && BC.fixed[3*j+2];
    { const fs = new Set(); active.forEach(c => { freeCands(c).forEach(x => {   // long straight candidates: from off-grid joints or support anchors; short ones otherwise
      if (!(offGrid(x.a) || offGrid(x.b) || anchor(x.a) || anchor(x.b)) && x.L > Rn) return;
      if (!fs.has(fKey(x))) { fs.add(fKey(x)); cand.push(x); } }); }); }
    const proj = c => {                            // w = Nᵀ b_c (b_c: unit bar vector)
      const w = new Float64Array(N.length);
      N.forEach((v, i) => { let s = 0; for (let t = 0; t < 3; t++) s += c.d[t]*(v[3*c.b+t] - v[3*c.a+t]); w[i] = s; });
      return w;
    };
    const reduce = bw => {                         // N ← basis of {x ∈ span N : bᵀx = 0} (Householder in mode space)
      const wn = vnorm(bw), u = Array.from(bw, x => x/wn); u[0] -= 1;
      const un = vnorm(u), k = N.length;
      if (un < 1e-14) { N = N.slice(1); return; }
      for (let i = 0; i < k; i++) u[i] /= un;
      const len = N[0].length, w = new Float64Array(len);   // H = I − 2uuᵀ:  (N H)_j = N_j − 2 u_j (Σ_i u_i N_i)
      for (let i = 0; i < k; i++) if (u[i] !== 0) { const Ni = N[i], ui = u[i]; for (let p = 0; p < len; p++) w[p] += ui*Ni[p]; }
      const N2 = [];
      for (let j = 1; j < k; j++) { const v = Float64Array.from(N[j]), c = 2*u[j]; if (c !== 0) for (let p = 0; p < len; p++) v[p] -= c*w[p]; N2.push(v); }
      N = N2;
    };
    const nodeMembers = c => {                     // one-cell members from inactive node c to active nodes
      const list = [];
      (incM.get(c) || []).forEach(m => { const [a, b] = GS.E[m]; if ((a === c && active.has(b)) || (b === c && active.has(a))) if (okPart(m)) list.push(m); });
      return list.map(m => { const [a, b] = GS.E[m], d = sub(X[b], X[a]), L = norm(d); return { a, b, parts: [m], L, d: d.map(x => x/L) }; }).concat(freeCands(c, Rn));
    };
    const unit = (j, t) => { const v = new Float64Array(3*X.length); v[3*j+t] = 1; return v; };
    const pushUnit = (j, t) => { N.push(unit(j, t)); };
    // lazy greedy: a candidate's score can only fall while the mechanism space shrinks, so stale scores are upper bounds
    const heap = [], hUp = i => { while (i > 0) { const p = (i - 1) >> 1; if (heap[p].ub >= heap[i].ub) break; [heap[p], heap[i]] = [heap[i], heap[p]]; i = p; } };
    const hDown = i => { for (;;) { const l = 2*i + 1, r = l + 1; let m = i; if (l < heap.length && heap[l].ub > heap[m].ub) m = l; if (r < heap.length && heap[r].ub > heap[m].ub) m = r; if (m === i) break; [heap[m], heap[i]] = [heap[i], heap[m]]; i = m; } };
    const hPush = c => { c.inHeap = true; heap.push(c); hUp(heap.length - 1); };
    const hPop = () => { const top = heap[0], last = heap.pop(); if (heap.length) { heap[0] = last; hDown(0); } top.inHeap = false; return top; };
    const gone = c => c.skip || c.parts.some(m => present.has(m) || blocked(m)) || (c.mids && c.mids.some(j => used.has(j))) || (c.free && freeUsed.has(fKey(c)));
    const revive = c => { c.ub = Infinity; if (!c.inHeap) hPush(c); };
    const sym = GS.sym, gOf = m => groupsOf.get(m) || [];
    const orbitOf = c0 => {                        // mirror orbit of a candidate; null if any image is missing or conflicts
      const ck = x => x.free ? 'n' + Math.min(x.a, x.b) + ',' + Math.max(x.a, x.b) : x.parts.slice().sort((u, v) => u - v).join(',');
      const out = [c0], seen = new Set([ck(c0)]);
      if (sym) for (let pass = 0; pass < 3; pass++) for (const x of out.slice()) for (let p = 0; p < sym.planes.length; p++) {
        const nm = sym.nodeMir[p], mm = sym.memMir[p], a = nm[x.a], b = nm[x.b], parts = x.parts.map(m => mm[m]), mids = (x.mids || []).map(j => nm[j]);
        if (a < 0 || b < 0 || parts.some(m => m < 0) || mids.some(j => j < 0) || !active.has(a) || !active.has(b)) return null;
        const img = { a: Math.min(a, b), b: Math.max(a, b), parts, mids, free: x.free }, k = ck(img); if (seen.has(k)) continue; seen.add(k);
        const d = sub(X[img.b], X[img.a]), L = norm(d); img.L = L; img.d = d.map(v => v/L); out.push(img);
      }
      const allParts = out.flatMap(c => c.parts), allMids = out.flatMap(c => c.mids || []);
      if (new Set(allParts).size !== allParts.length || allParts.some(m => present.has(m) || blocked(m)) || allMids.some(j => used.has(j)) || new Set(allMids).size !== allMids.length) return null;   // mirror chains sharing a mid node would cross there
      if (out.some(c => (c.free && freeUsed.has(fKey(c))) || crossesAny(c))) return null;
      for (let i = 0; i < out.length; i++) for (let j = i + 1; j < out.length; j++) {   // mirror images must not cross each other
        const A = out[i], B = out[j]; if ((!A.free && !B.free) || A.a === B.a || A.a === B.b || A.b === B.a || A.b === B.b) continue;
        if (segsCross(X[A.a], X[A.b], X[B.a], X[B.b], tolF)) return null;
      }
      const gs = allParts.flatMap(gOf); if (new Set(gs).size !== gs.length) return null;   // two images in one crossing group
      return out;
    };
    for (let guard = 0; N.length && guard < 200; guard++) {
      for (;;) {                                   // greedy member addition
        if (!N.length) break;
        if (guard === 0 && !heap.length) cand.forEach(revive);
        let best = null;
        while (heap.length) {
          const c = hPop(); if (gone(c)) continue;
          const sc = vnorm(proj(c))*Math.sqrt(h/c.L); if (!(sc > 1e-6)) { c.ub = 0; continue; }   // dormant until a touching joint is activated
          if (!heap.length || sc >= heap[0].ub - 1e-12) { best = c; break; }
          c.ub = sc; hPush(c);
        }
        if (!best) { if (globalThis.GSTO_DEBUG) { let mx = 0, why = ''; for (const c of cand) { const w = proj(c), sc = vnorm(w); if (sc > mx) { mx = sc; why = (c.skip ? 'skip(orbit)' : '') + (c.parts.some(m => present.has(m)) ? 'present' : '') + (c.parts.some(m => blocked(m)) ? 'blocked' : '') + ((c.mids && c.mids.some(j => used.has(j))) ? 'mid-used' : ''); } } console.log('stab: no member; best raw proj', mx.toExponential(2), why); } break; }
        const orb = orbitOf(best);                 // add the candidate together with its mirror images
        if (!orb) { best.skip = true; continue; }
        if (!best.inHeap) { best.ub = 0; }
        orb.forEach(c => {
          c.parts.forEach(m => { present.add(m); GS.E[m].forEach(j => used.add(j)); });
          if (c.free) freeUsed.add(fKey(c));
          out.push({ a: c.a, b: c.b, A: sizeOf(c.a, c.b), q: S.members[0].q.map(() => 0), gs: c.parts, stab: true, L: c.L, long: !!c.free });
          if (N.length) { const w = proj(c); if (vnorm(w) > 1e-9) reduce(w); }
        });
      }
      if (!N.length) break;
      // node activation: pick an inactive grid node whose members reduce the net nullity the most
      const near = new Set();
      active.forEach(a => (incM.get(a) || []).forEach(m => { const [p, q] = GS.E[m], o = p === a ? q : p; if (!active.has(o) && !dead.has(o)) near.add(o); }));
      active.forEach(a => { if (offGrid(a)) for (let j = 0; j < GS.nGrid; j++) if (!active.has(j) && !dead.has(j) && norm(sub(X[j], X[a])) <= Rn) near.add(j); });
      let bestC = -1, bestGain = 0, bestCost = Infinity;
      const evalNode = (c, allowMid) => {
        const isMid = midOf.has(c) && !active.has(c);
        if (isMid && !allowMid) return;
        if (used.has(c) && !isMid) return;
        const nFree = [0, 1, 2].filter(t => !BC.fixed[3*c+t]).length, mem = nodeMembers(c);
        if (isMid) { const hm = halfOf(c); mem.push({ a: hm.b, b: c, L: hm.L, d: hm.d.map(v => -v), half: true }); }
        if (mem.length < Math.max(1, nFree)) return;   // support nodes need no own bracing
        const ext = N.concat([0,1,2].filter(t => !BC.fixed[3*c+t]).map(t => unit(c, t)));
        const W = mem.map(e => ext.map(v => { let s = 0; for (let t = 0; t < 3; t++) s += e.d[t]*(v[3*e.b+t] - v[3*e.a+t]); return s; }));
        const k = W.length, Gm = new Float64Array(k*k);
        for (let i = 0; i < k; i++) for (let j = 0; j < k; j++) Gm[i*k+j] = W[i].reduce((s, x, p) => s + x*W[j][p], 0);
        const ev = symEig(Gm, k).values, emax = Math.max(...ev), rk = Array.from(ev).filter(x => x > 1e-9*emax).length;
        const gain = rk - (ext.length - N.length), cost = mem.reduce((s, e) => s + (e.half ? 0 : e.L), 0);
        if (gain > bestGain || (gain === bestGain && gain > 0 && cost < bestCost)) { bestC = c; bestGain = gain; bestCost = cost; }
      };
      near.forEach(c => evalNode(c, false));        // plain grid nodes first
      if (bestC < 0) midOf.forEach((i, c) => { if (!active.has(c)) evalNode(c, true); });   // then interior points of long members
      let pick = bestC >= 0 ? [bestC] : null;
      if (!pick) {                                  // two-node lookahead: some mechanisms need a new node pair (e.g. sway towards a support line)
        const nbr = new Map(); GS.E.forEach(([a, b], m) => { if (GS.kind[m] >= 4) return; [[a, b], [b, a]].forEach(([u, v]) => { if (!nbr.has(u)) nbr.set(u, []); nbr.get(u).push(v); }); });
        const evalSet = set => {
          const inSet = new Set(set), mem = [], seenM = new Set();
          set.forEach(c => { if (midOf.has(c) && !active.has(c)) { const hm = halfOf(c); mem.push({ a: hm.b, b: c, L: 0, d: hm.d.map(v => -v) }); } });
          set.forEach(c => (incM.get(c) || []).forEach(m => { const [a, b] = GS.E[m];
            if (seenM.has(m) || !okPart(m)) return;
            const o = a === c ? b : b === c ? a : -1; if (o < 0 || !(active.has(o) || inSet.has(o))) return;
            seenM.add(m); const d = sub(X[b], X[a]), L = norm(d); mem.push({ a, b, L, d: d.map(v => v/L) });
          }));
          set.forEach(c => freeCands(c, Rn).forEach(x => mem.push(x)));
          const ext = N.concat(set.flatMap(c => [0, 1, 2].filter(t => !BC.fixed[3*c+t]).map(t => unit(c, t))));
          if (!mem.length) return { gain: 0 };
          const W = mem.map(e => ext.map(v => { let s = 0; for (let t = 0; t < 3; t++) s += e.d[t]*(v[3*e.b+t] - v[3*e.a+t]); return s; }));
          const k = W.length, Gm = new Float64Array(k*k);
          for (let i = 0; i < k; i++) for (let j = 0; j < k; j++) Gm[i*k+j] = W[i].reduce((s, x, p) => s + x*W[j][p], 0);
          const ev = symEig(Gm, k).values, emax = Math.max(...ev);
          return { gain: Array.from(ev).filter(x => x > 1e-9*emax).length - (ext.length - N.length), cost: mem.reduce((s, e) => s + e.L, 0) };
        };
        let bg = 0, bc = Infinity;
        const isFreeMid = c => midOf.has(c) && !active.has(c);
        near.forEach(c => { if (used.has(c) && !isFreeMid(c)) return; (nbr.get(c) || []).forEach(d => {
          if (d <= c && near.has(d)) return; if (active.has(d) || (used.has(d) && !isFreeMid(d))) return;
          const r = evalSet([c, d]); if (r.gain > bg || (r.gain === bg && r.gain > 0 && r.cost < bc)) { bg = r.gain; bc = r.cost; pick = [c, d]; } }); });
        if (globalThis.GSTO_DEBUG) console.log('stab: pair lookahead gain', bg, pick);
      }
      if (!pick) break;
      const newNodes = new Set(pick);
      if (sym) for (let pass = 0; pass < 3; pass++) [...newNodes].forEach(c => sym.nodeMir.forEach(nm => { if (nm[c] >= 0 && !active.has(nm[c]) && !used.has(nm[c])) newNodes.add(nm[c]); }));
      if (sym) for (let pass = 0; pass < 3; pass++) [...newNodes].forEach(c => sym.nodeMir.forEach(nm => { const d = nm[c]; if (d >= 0 && !active.has(d) && midOf.has(d)) newNodes.add(d); }));
      newNodes.forEach(c => {
        if (midOf.has(c)) {                         // interior point: split the member there and constrain the point along the bar
          splits.push([c, midOf.get(c)]);
          const hm = halfOf(c), w = new Float64Array(N.length);
          N.forEach((v, i) => { let s = 0; for (let t = 0; t < 3; t++) s += hm.d[t]*(v[3*hm.b+t] - v[3*c+t]); w[i] = s; });
          [0,1,2].filter(t => !BC.fixed[3*c+t]).forEach(t => pushUnit(c, t));
          active.add(c); used.add(c);
          const w2 = new Float64Array(N.length); N.forEach((v, i) => { let s = 0; for (let t = 0; t < 3; t++) s += hm.d[t]*(v[3*hm.b+t] - v[3*c+t]); w2[i] = s; });
          if (vnorm(w2) > 1e-9) reduce(w2);
          return;
        }
        active.add(c); used.add(c); [0,1,2].filter(t => !BC.fixed[3*c+t]).forEach(t => pushUnit(c, t));
      });
      cand.forEach(c => { if (newNodes.has(c.a) || newNodes.has(c.b)) revive(c); });
      newNodes.forEach(c => { nodeMembers(c).forEach(x => { if (!x.free) { cand.push(x); revive(x); } }); freeCands(c).forEach(x => { cand.push(x); revive(x); }); });
    }
    out.splits = splits;
    return out;
  }


  // ---------- folding analysis: rigid triangular faces, hinges between them ----------
  function foldAnalysis(GS, S, stab) {
    const X = GS.X, adj = new Map(), key = (a, b) => a < b ? a+','+b : b+','+a, edge = new Set();
    S.members.forEach(e => { [[e.a, e.b], [e.b, e.a]].forEach(([u, v]) => { if (!adj.has(u)) adj.set(u, new Set()); adj.get(u).add(v); }); edge.add(key(e.a, e.b)); });
    const tris = [];
    adj.forEach((nb, a) => nb.forEach(b => { if (b <= a) return; nb.forEach(c => { if (c <= b || !adj.get(b).has(c)) return;
      if (norm(cross(sub(X[b], X[a]), sub(X[c], X[a]))) > 1e-9) tris.push([a, b, c]); }); }));
    const out = { triangles: tris.length, hinges: [], bodies: [], triBody: [] };
    if (!stab || stab.skipped || !stab.m || !stab.modes || !tris.length) { out.bodyCount = tris.length ? 1 : 0; return out; }
    const modes = stab.modes, scale = amaxOf(S.members.map(e => e.L));
    const omega = (t, v) => {                       // rigid rotation rate of a triangle in mode v (least squares)
      const [a, b, c] = t, M = new Float64Array(9), r = [0, 0, 0];
      for (const j of [b, c]) {
        const d = sub(X[j], X[a]), du = [0, 1, 2].map(k => v[3*j+k] - v[3*a+k]);
        const A = [[0, d[2], -d[1]], [-d[2], 0, d[0]], [d[1], -d[0], 0]];   // du = ω × d = A ω
        for (let i = 0; i < 3; i++) for (let k = 0; k < 3; k++) { for (let l = 0; l < 3; l++) M[3*k+l] += A[i][k]*A[i][l]; r[k] += A[i][k]*du[i]; }
      }
      const det = M[0]*(M[4]*M[8]-M[5]*M[7]) - M[1]*(M[3]*M[8]-M[5]*M[6]) + M[2]*(M[3]*M[7]-M[4]*M[6]);
      const inv = [(M[4]*M[8]-M[5]*M[7]), -(M[1]*M[8]-M[2]*M[7]), (M[1]*M[5]-M[2]*M[4]), -(M[3]*M[8]-M[5]*M[6]), (M[0]*M[8]-M[2]*M[6]), -(M[0]*M[5]-M[2]*M[3]), (M[3]*M[7]-M[4]*M[6]), -(M[0]*M[7]-M[1]*M[6]), (M[0]*M[4]-M[1]*M[3])];
      return [0, 1, 2].map(k => (inv[3*k]*r[0] + inv[3*k+1]*r[1] + inv[3*k+2]*r[2])/det);
    };
    const W = tris.map(t => modes.map(v => omega(t, v)));
    const par = tris.map((_, i) => i), find = x => { while (par[x] !== x) x = par[x] = par[par[x]]; return x; };
    const byEdge = new Map();
    tris.forEach((t, i) => [[t[0], t[1]], [t[1], t[2]], [t[0], t[2]]].forEach(([a, b]) => { const k = key(a, b); if (!byEdge.has(k)) byEdge.set(k, []); byEdge.get(k).push(i); }));
    const same = (i, j) => W[i].every((w, mi) => Math.hypot(w[0]-W[j][mi][0], w[1]-W[j][mi][1], w[2]-W[j][mi][2])*scale < 1e-6);
    const pairs = [];
    byEdge.forEach((list, k) => { for (let x = 0; x < list.length; x++) for (let y = x+1; y < list.length; y++) { const i = list[x], j = list[y]; if (same(i, j)) par[find(i)] = find(j); else pairs.push([i, j, k]); } });
    const seen = new Set();
    pairs.forEach(([i, j, k]) => { if (find(i) !== find(j) && !seen.has(k)) { seen.add(k); out.hinges.push(k.split(',').map(Number)); } });
    const ids = new Map(); out.triBody = tris.map((t, i) => { const r = find(i); if (!ids.has(r)) ids.set(r, ids.size); return ids.get(r); });
    out.tris = tris; out.bodyCount = ids.size;
    return out;
  }


  function splitMembers(GS, S, splits) {         // split merged members at interior grid points
    if (!splits || !splits.length) return 0;
    const X = GS.X, by = new Map(); splits.forEach(([c, i]) => { if (!by.has(i)) by.set(i, new Set()); by.get(i).add(c); });
    const mem = [];
    S.members.forEach((e, i) => {
      if (!by.has(i)) { mem.push(e); return; }
      const u = sub(X[e.b], X[e.a]).map(v => v/e.L), t = j => dot(sub(X[j], X[e.a]), u);
      const cuts = [e.a, ...[...by.get(i)].sort((p, q) => t(p) - t(q)), e.b];
      for (let k = 0; k + 1 < cuts.length; k++) {
        const lo = t(cuts[k]) - 1e-9, hi = t(cuts[k+1]) + 1e-9;
        const parts = (e.gs || []).filter(m => GS.E[m].every(j => t(j) >= lo && t(j) <= hi));
        mem.push({ ...e, a: Math.min(cuts[k], cuts[k+1]), b: Math.max(cuts[k], cuts[k+1]), gs: parts, L: hi - lo - 2e-9, braced: true });
      }
    });
    S.members = mem; return by.size;
  }

  // ---------- member buckling (Euler, pinned ends) ----------
  function sectionI(P, A) {                        // second moment of area from area for the assumed section
    const o = P.options || {};
    if (o.section === 'solid') return A*A/(4*Math.PI);
    const k = o.DoverT || 20, D2 = A*k*k/(Math.PI*(k - 1)), D = Math.sqrt(D2), di = D*(1 - 2/k);
    return Math.PI/64*(D**4 - di**4);
  }
  function bucklingFix(P, GS, BC, S, noBrace = new Set()) {
    const o = P.options || {}, E = P.material.E, g = o.gammaB ?? 2, X = GS.X;
    const c1 = sectionI(P, 1);                      // I = c1·A² for both section models
    const out = { split: 0, braceNodes: [], enlarged: 0, over: 0, origin: new Map() };
    const mem = [];
    S.members.forEach((e, ei) => {
      const Pc = Math.max(0, ...e.q.map(v => -v));
      if (!(Pc > 0) || e.stab) { mem.push(e); return; }
      const Pcr = Math.PI**2*E*sectionI(P, e.A)/(e.L*e.L);
      if (g*Pc <= Pcr) { mem.push(e); return; }
      out.over++;
      const Lcr = Math.PI*Math.sqrt(E*sectionI(P, e.A)/(g*Pc));
      let pieces = [e];
      const u = sub(X[e.b], X[e.a]).map(v => v/e.L), tt = j => dot(sub(X[j], X[e.a]), u), tl = 1e-6*P.grid.h;
      const chain = e.gs && !e.long && e.gs.length > 1;   // straight chain of one-cell members: brace at its grid points
      if (!chain && e.long && e.L > 1.5*P.grid.h && !noBrace.has(ei) && BC) {   // long member: brace at new joints at equal spacing
        const k = Math.ceil(e.L/Lcr - 1e-9);
        if (k > 1) {
          const ids = [e.a]; for (let i = 1; i < k; i++) ids.push(addNode(GS, BC, X[e.a].map((v, t) => v + u[t]*i*e.L/k))); ids.push(e.b);
          pieces = [];
          for (let i = 0; i + 1 < ids.length; i++) { pieces.push({ ...e, a: Math.min(ids[i], ids[i+1]), b: Math.max(ids[i], ids[i+1]), L: e.L/k, braced: true }); if (i > 0) { out.braceNodes.push(ids[i]); out.origin.set(ids[i], ei); } }
          out.split++;
        }
      }
      if (chain && !noBrace.has(ei)) {
        const nodes = [...new Set([e.a, e.b, ...e.gs.flatMap(m => GS.E[m])])].filter(j => tt(j) > -tl && tt(j) < e.L + tl).sort((p, q) => tt(p) - tt(q));
        const pos = nodes.map(j => dot(sub(X[j], X[e.a]), u)), k = Math.ceil(e.L/Lcr - 1e-9), cut = [0];
        for (let i = 1; i < k; i++) {               // grid points nearest to equal spacing (mirror-consistent)
          const target = i*e.L/k; let best = -1, bd = Infinity;
          pos.forEach((p, idx) => { const d = Math.abs(p - target); if (idx > cut[cut.length-1] && idx < nodes.length - 1 && d < bd - 1e-9) { bd = d; best = idx; } });
          if (best > 0) cut.push(best);
        }
        cut.push(nodes.length - 1);
        if (cut.length > 2) {
          pieces = [];
          for (let i = 0; i + 1 < cut.length; i++) {
            const a = nodes[cut[i]], b = nodes[cut[i+1]], parts = e.gs.filter(m => { const [p, q] = GS.E[m]; const lo = pos[cut[i]] - 1e-9, hi = pos[cut[i+1]] + 1e-9; return [p, q].every(j => { const t = dot(sub(X[j], X[e.a]), u); return t >= lo && t <= hi; }); });
            pieces.push({ ...e, a: Math.min(a, b), b: Math.max(a, b), gs: parts, L: pos[cut[i+1]] - pos[cut[i]], braced: true });
            if (i > 0) { out.braceNodes.push(a); out.origin.set(a, ei); }
          }
          out.split++;
        }
      }
      pieces.forEach(pc => {                         // pieces that still buckle (one-cell or no grid point): enlarge the section
        const Pcr2 = Math.PI**2*E*sectionI(P, pc.A)/(pc.L*pc.L);
        if (g*Pc > Pcr2) {
          const Areq = Math.sqrt(g*Pc*pc.L*pc.L/(Math.PI**2*E*c1)); pc.A = Areq; pc.enlarged = true; out.enlarged++;
        }
        mem.push(pc);
      });
    });
    S.members = mem; S.volume = mem.reduce((v, e) => v + e.L*e.A, 0);
    out.braceNodes = [...new Set(out.braceNodes)];
    return out;
  }

  // ---------- linear stiffness analysis (sparse PCG) ----------
  function stiffnessAnalysis(P, GS, BC, S, softNodes = []) {
    const X = GS.X, mem = S.members, E = P.material.E, n = X.length;
    const nodes = [...new Set(mem.flatMap(e => [e.a, e.b]))], dof = new Int32Array(3*n).fill(-1); let nf = 0;
    nodes.forEach(j => { for (let t = 0; t < 3; t++) if (!BC.fixed[3*j+t]) dof[3*j+t] = nf++; });
    const rows = Array.from({length: nf}, () => new Map());
    mem.forEach(e => {
      const d = sub(X[e.b], X[e.a]), l = norm(d), u = d.map(x => x/l), k = E*e.A/l;
      const idx = [3*e.a, 3*e.a+1, 3*e.a+2, 3*e.b, 3*e.b+1, 3*e.b+2].map(i => dof[i]), v = [-u[0],-u[1],-u[2],u[0],u[1],u[2]];
      for (let p = 0; p < 6; p++) { if (idx[p] < 0) continue; const row = rows[idx[p]]; for (let q = 0; q < 6; q++) if (idx[q] >= 0) row.set(idx[q], (row.get(idx[q]) || 0) + k*v[p]*v[q]); }
    });
    if (softNodes.length) {                        // weak springs at mechanism nodes (display only)
      let dm = 0, cnt = 0; rows.forEach((r, i) => { const v = r.get(i); if (v) { dm += v; cnt++; } });
      const ks = 1e-3*dm/Math.max(cnt, 1);
      softNodes.forEach(j => { for (let t = 0; t < 3; t++) { const d = dof[3*j+t]; if (d >= 0) rows[d].set(d, (rows[d].get(d) || 0) + ks); } });
    }
    const rp = new Int32Array(nf+1); rows.forEach((r, i) => rp[i+1] = rp[i] + r.size);
    const ci = new Int32Array(rp[nf]), cv = new Float64Array(rp[nf]), diag = new Float64Array(nf);
    rows.forEach((r, i) => { let p = rp[i]; for (const [c, v] of r) { ci[p] = c; cv[p++] = v; if (c === i) diag[i] = v; } });
    const mv = (x, y) => { for (let i = 0; i < nf; i++) { let s = 0; for (let p = rp[i]; p < rp[i+1]; p++) s += cv[p]*x[ci[p]]; y[i] = s; } };
    const results = BC.cases.map(cs => {
      const f = new Float64Array(nf); let fn = 0;
      for (let i = 0; i < 3*n; i++) if (dof[i] >= 0) { f[dof[i]] = cs.f[i]; fn += cs.f[i]*cs.f[i]; }
      fn = Math.sqrt(fn);
      const x = new Float64Array(nf), r = f.slice(), z = r.map((v, i) => v/diag[i]), p = z.slice(), Ap = new Float64Array(nf);
      let rz = r.reduce((s, v, i) => s + v*z[i], 0), it = 0, rn = fn;
      for (; it < 20*nf + 100 && fn > 0; it++) {
        mv(p, Ap); const al = rz / p.reduce((s, v, i) => s + v*Ap[i], 0);
        for (let i = 0; i < nf; i++) { x[i] += al*p[i]; r[i] -= al*Ap[i]; }
        rn = Math.sqrt(r.reduce((s, v) => s + v*v, 0)); if (rn < 1e-10*fn) break;
        for (let i = 0; i < nf; i++) z[i] = r[i]/diag[i];
        const rz2 = r.reduce((s, v, i) => s + v*z[i], 0); for (let i = 0; i < nf; i++) p[i] = z[i] + (rz2/rz)*p[i]; rz = rz2;
      }
      const U = new Float64Array(3*n); for (let i = 0; i < 3*n; i++) if (dof[i] >= 0) U[i] = x[dof[i]];
      const q = mem.map(e => { const d = sub(X[e.b], X[e.a]), l = norm(d); let s = 0; for (let t = 0; t < 3; t++) s += d[t]/l*(U[3*e.b+t]-U[3*e.a+t]); return E*e.A/l*s; });
      let umax = 0; for (let j = 0; j < n; j++) umax = Math.max(umax, Math.hypot(U[3*j], U[3*j+1], U[3*j+2]));
      const comp = cs.f.reduce((s, v, i) => s + v*U[i], 0);
      const util = q.map((qq, i) => qq >= 0 ? qq/(mem[i].A*P.material.sigmaT) : -qq/(mem[i].A*P.material.sigmaC));
      const gB = P.options?.gammaB ?? 2;
      const buckUtil = q.reduce((mx, qq, i) => qq < 0 ? Math.max(mx, -qq*gB*mem[i].L**2/(Math.PI**2*P.material.E*sectionI(P, mem[i].A))) : mx, 0);
      return { name: cs.name, U, q, umax, compliance: comp, converged: rn < 1e-8*Math.max(fn,1e-30), iters: it, utilMax: Math.max(...util), buckUtil };
    });
    return { nf, results, regularized: softNodes.length > 0 };
  }

  // ---------- Abaqus .inp export ----------
  function toInp(P, GS, BC, S) {
    const X = GS.X, used = [...new Set(S.members.flatMap(e => [e.a, e.b]))].sort((a,b) => a-b), nid = new Map();
    used.forEach((j, i) => nid.set(j, i+1));
    const L = ['*HEADING', 'GSTO truss layout optimization result (units: mm, N, MPa)', '*NODE'];
    used.forEach(j => L.push(`${nid.get(j)}, ${X[j].map(v => v.toFixed(6)).join(', ')}`));
    L.push('*ELEMENT, TYPE=T3D2, ELSET=TRUSS');
    S.members.forEach((e, i) => L.push(`${i+1}, ${nid.get(e.a)}, ${nid.get(e.b)}`));
    S.members.forEach((e, i) => { L.push(`*ELSET, ELSET=M${i+1}`, `${i+1}`, `*SOLID SECTION, ELSET=M${i+1}, MATERIAL=${e.rod ? 'MAT_ROD' : 'MAT'}`, `${e.A.toExponential(6)}`); });
    const grp = (name, f) => { const ids = S.members.map((e, i) => f(e) ? i+1 : 0).filter(Boolean); if (!ids.length) return; L.push(`*ELSET, ELSET=${name}`); for (let i = 0; i < ids.length; i += 16) L.push(ids.slice(i, i+16).join(', ')); };
    grp('RODS', e => e.rod); grp('BRACING', e => e.stab); grp('DECK', e => e.deck); grp('MAIN', e => !e.stab && !e.rod && !e.deck);
    L.push('*MATERIAL, NAME=MAT', '*ELASTIC', `${P.material.E}, ${P.material.nu ?? 0.3}`);
    if (S.members.some(e => e.rod)) L.push('** MAT_ROD: tension-only rods (all load cases in tension in the optimisation); add *NO COMPRESSION if desired', '*MATERIAL, NAME=MAT_ROD', '*ELASTIC', `${P.material.E}, ${P.material.nu ?? 0.3}`);
    L.push('*BOUNDARY');
    used.forEach(j => { for (let t = 0; t < 3; t++) if (BC.fixed[3*j+t]) L.push(`${nid.get(j)}, ${t+1}, ${t+1}`); });
    BC.cases.forEach(cs => {
      L.push(`*STEP, NAME=${cs.name}`, '*STATIC', '*CLOAD, OP=NEW');
      used.forEach(j => { for (let t = 0; t < 3; t++) { const v = cs.f[3*j+t]; if (Math.abs(v) > 0) L.push(`${nid.get(j)}, ${t+1}, ${v.toExponential(6)}`); } });
      L.push('*OUTPUT, FIELD', '*NODE OUTPUT', 'U, RF', '*ELEMENT OUTPUT', 'S, SF', '*END STEP');
    });
    return L.join('\n') + '\n';
  }

  async function run(highs, P, log = () => {}) {
    const t0 = Date.now();
    const GS = buildGroundStructure(P);
    log(`Ground structure: 노드 ${GS.X.length} (격자 ${GS.nGrid}, 추가 ${GS.X.length-GS.nGrid}), 후보 부재 ${GS.E.length}, 교차군 ${GS.groups.length}`);
    const BC = buildBC(P, GS); BC.warn.forEach(w => log('경고: ' + w));
    const symPlanes = P.options?.symAuto !== false ? detectSymmetry(P, GS, BC) : (P.options?.symmetry || [0, 0, 0]);
    const rebuildSym = () => { GS.sym = symmetryMaps({ ...P, options: { ...P.options, symmetry: symPlanes, symExact: P.options?.symAuto !== false } }, GS); };   // detected symmetry is exact → reduced LP   // after nodes / members are added
    rebuildSym();
    if (P.options?.symAuto !== false) log(symPlanes.some(Boolean) ? `대칭 감지 (모델 기준): ${symPlanes.map((v, t) => v ? 'xyz'[t] : '').join('')}` : '대칭 감지 (모델 기준): 없음');
    if (GS.sym) log(`대칭 강제 (${GS.sym.planes.map(t => 'xyz'[t]).join(', ')} 중앙면): 대칭 궤도 ${GS.sym.orbits.length}개${GS.sym.missing ? `, 거울상이 없는 부재 ${GS.sym.missing}개 (장애물·추가 노드가 비대칭)` : ''}`);
    if (BC.share.length) log(`반력 분담: 지지 영역 ${BC.share.length}곳 각각 하중의 ${(100*P.options.supportShare).toFixed(0)}% 이상`);
    if (!BC.cases.length) return { ok: false, msg: '하중 케이스가 없습니다', GS };
    if (BC.flat >= 0) log(`2D 모드: ${['YZ', 'XZ', 'XY'][BC.flat]} 평면 트러스 (평면 밖 ${'xyz'[BC.flat]} 방향은 횡지지된 것으로 가정)`);
    log(`자유도 ${BC.nf}, 하중 케이스 ${BC.cases.map(c => c.name).join(', ')}`);
    let ob = P.options?.objective || 'roundup';
    if (ob === 'uniform') { log('균일 단면 모드는 없어졌습니다 → 최소 단면(올림 + 재최적화)으로 실행합니다'); ob = 'roundup'; }
    let grow = null, geom = null;
    const tm = {}, tick = (k, t) => { tm[k] = (tm[k] || 0) + Date.now() - t; };   // stage timing [ms]
    let tS = Date.now();
    if ((P.options?.connect || 'full') === 'full') { grow = growGroundStructure(highs, P, GS, BC, log); rebuildSym(); }
    tick('grow', tS); tS = Date.now();
    let opt = ob === 'roundup' ? solveRoundUp(highs, P, GS, BC, log) : ob === 'minA' ? solveMinArea(highs, P, GS, BC, log) : optimize(highs, P, GS, BC, log);
    if (!opt.ok) return { ok: false, msg: opt.msg, GS, BC };
    tick('opt', tS); tS = Date.now();
    if (P.options?.geomOpt) {
      const g = geometryOptimize(highs, P, GS, BC, opt, log); if (g) { opt = g; geom = g.geom; rebuildSym(); }
    }
    tick('geom', tS); tS = Date.now();
    const probe = P.options?._probe || (() => {});   // debug hook: probe(stage, GS, members)
    probe('opt', GS, GS.E.map(([a, b], m) => ({ a, b, A: opt.a[m] })).filter(e => e.A > 1e-9));
    const S = extractStructure(P, GS, BC, opt);
    if (BC.deck && BC.deck.size) {                  // load-surface members: flag, and add any the optimiser left out
      const have = new Set(S.members.flatMap(e => e.gs || [])), dA = opt.thrC ?? (P.options?.stabArea ?? 0.01)*amaxOf(opt.a);
      S.members.forEach(e => { if ((e.gs || []).some(m => BC.deck.has(m))) e.deck = true; });
      let addD = 0; BC.deck.forEach(m => { if (!have.has(m)) { const [a, b] = GS.E[m]; S.members.push({ a, b, A: dA, q: BC.cases.map(() => 0), gs: [m], L: GS.L[m], deck: true }); addD++; } });
      S.volume = S.members.reduce((v, e) => v + e.L*e.A, 0);
      log(`하중면 부재 ${BC.deck.size}개 (최적화에 포함 ${BC.deck.size - addD}, 후처리로 추가 ${addD}, 단면 ≥ ${dA.toFixed(3)} mm²)`);
    }
    const nCrossSplit = splitCrossings(P, GS, BC, S);
    if (nCrossSplit) { log(`장거리 부재 교차 ${nCrossSplit}곳에 절점 생성 (부재 분할)`); rebuildSym(); }
    probe('extract', GS, S.members);
    S.aStruct = ob === 'minA' || ob === 'roundup' ? (P.options.braceMin > 0 ? P.options.braceMin*opt.amax0 : opt.thrC) : undefined;   // locking / bracing member minimum
    const nRod = S.members.filter(e => e.rod).length;
    log(`결과 부재 ${S.members.length} (직선 병합 ${S.merged}${nRod ? `, 텐션 로드 ${nRod}` : ''}), 체적 ${S.volume.toExponential(4)} mm³`);
    const S0 = S.members.map(e => ({ ...e })), V0 = S.volume;
    let buck = null, stab, added = [], fold0 = null, braceSplit = 0;
    const stage = noBrace => {                      // buckling fix → DOF check → locking / lateral bracing
      S.members = S0.map(e => ({ ...e })); S.volume = V0; buck = null; added = []; braceSplit = 0;
      if (P.options?.buckling !== false) { buck = bucklingFix(P, GS, BC, S, noBrace); rebuildSym(); }
      stab = stabilityCheck(GS, BC, S);
      fold0 = stab.skipped ? null : foldAnalysis(GS, S, stab);
      if (!stab.skipped && stab.m > 0 && P.options?.stabilize !== false) {
        added = stabilizeGreedy(P, GS, BC, S, stab);
        const nSplit = splitMembers(GS, S, added.splits);
        if (nSplit) braceSplit += nSplit;
        if (added.length || nSplit) { added.forEach(e => S.members.push(e)); S.volume = S.members.reduce((v, e) => v + e.L*e.A, 0); stab = stabilityCheck(GS, BC, S); }
      }
    };
    stage(new Set());
    probe('lock', GS, S.members);
    if (buck && buck.braceNodes.length && stab.m > 0) {   // brace points that could not be tied back: enlarge those members instead
      const bad = new Set(stab.mechNodes.map(o => o.node).filter(j => buck.origin.has(j)).map(j => buck.origin.get(j)));
      if (bad.size) { log(`횡보강 불가 (주변에 붙잡을 구조 없음): 긴 압축재 ${bad.size}개는 단면 증대로 전환`); stage(bad); }
    }
    if (buck) log(buck.over ? `좌굴 검토: 초과 ${buck.over}개 → 긴 봉 ${buck.split}개를 중간 격자점 ${buck.braceNodes.length}곳에서 횡보강, 단면 증대 ${buck.enlarged}개${S.members.some(e => e.buckleOver) ? `, 균일 단면이라 증대하지 못한 부재 ${S.members.filter(e => e.buckleOver).length}개` : ''}` : '좌굴 검토: 초과 부재 없음');
    if (fold0 && fold0.hinges.length) log(`접힘 분석(잠금 전): 강체 삼각형 ${fold0.triangles}개, 강체 덩어리 ${fold0.bodyCount}개, 자유 힌지 ${fold0.hinges.length}개`);
    if (added.length) { const As = added.map(e => e.A); log(`잠금·횡보강 부재 ${added.length}개 추가 (단면 ${Math.min(...As).toFixed(2)}–${Math.max(...As).toFixed(2)} mm²)${braceSplit ? `, 긴 부재 ${braceSplit}개의 중간점을 절점으로 사용` : ''}`); }
    log(stab.skipped ? stab.msg : `자유도 판정: n_f=${stab.nf}, r=${stab.r}, m=${stab.m}, s=${stab.s}`);
    if (!stab.skipped && stab.m > 0) log(`남은 메커니즘 ${stab.m}개: 참여 노드 ${stab.mechNodes.length}개를 강접합 후보로 표시합니다`);
    const fold = stab.skipped ? null : foldAnalysis(GS, S, stab);
    if (fold) log(`접힘 분석: 강체 삼각형 ${fold.triangles}개, 자유 힌지 ${fold.hinges.length}개${fold.hinges.length ? ` (강체 덩어리 ${fold.bodyCount}개)` : ' — 접히는 삼각형 면 없음'}`);
    const onF = presence(opt, P.options?.areaCut ?? 1e-3);
    const crossLeft = GS.groups.filter(g => g.filter(onF).length > 1).length;
    tick('stab', tS); tS = Date.now();
    let fem = null, fsd = 0;
    if (!stab.skipped) {
      const soft = stab.m > 0 ? stab.mechNodes.map(o => o.node) : [];
      const sT = P.material.sigmaT, sC = P.material.sigmaC, E = P.material.E, gB = P.options?.gammaB ?? 2, c1 = sectionI(P, 1);
      const fsdOn = P.options?.fsd !== false;
      for (let it = 0; ; it++) {                    // fully stressed design: grow sections until the linear check passes
        fem = stiffnessAnalysis(P, GS, BC, S, soft);
        if (!fsdOn || it >= 30) break;
        let grown = 0;
        S.members.forEach((e, i) => {
          let req = 0;
          fem.results.forEach(r => { const q = r.q[i];
            if (q >= 0) req = Math.max(req, q/sT);
            else { req = Math.max(req, -q/sC); if (P.options?.buckling !== false) req = Math.max(req, Math.sqrt(-q*gB*e.L*e.L/(Math.PI**2*E*c1))); if (e.rod && -q > 1e-6*sC*e.A) e.rod = false; } });
          if (req > e.A*(1 + 1e-4)) { e.A = req*1.0001; grown++; }
        });
        if (!grown) break;
        fsd++;
      }
      S.volume = S.members.reduce((v, e) => v + e.L*e.A, 0);
      probe('fsd', GS, S.members);
      if (fsd) log(`선형 검증 단면 보정(FSD): ${fsd}회 반복, 체적 ${S.volume.toExponential(4)} mm³`);
      const cat = opt.catalog || [];
      if (cat.length) {                             // round every member up to the next catalogue section
        let over = 0; S.members.forEach(e => { const c = cat.find(v => v >= e.A*(1 - 1e-9)); if (c === undefined) { e.overCat = true; over++; } else e.A = c; });
        fem = stiffnessAnalysis(P, GS, BC, S, soft);
        for (let it = 0; it < 30; it++) {             // stiffer members attract force: bump any member still over its limit to the next size
          let bumped = 0;
          S.members.forEach((e, i) => { let req = 0; fem.results.forEach(r => { const q = r.q[i]; req = Math.max(req, q >= 0 ? q/sT : -q/sC); if (q < 0 && P.options?.buckling !== false) req = Math.max(req, Math.sqrt(-q*gB*e.L*e.L/(Math.PI**2*E*c1))); });
            if (req > e.A*(1 + 1e-6)) { const c = cat.find(v => v >= req*(1 - 1e-9)) ?? cat.find(v => v > e.A); if (c !== undefined && c > e.A) { e.A = c; bumped++; } } });
          if (!bumped) break;
          fem = stiffnessAnalysis(P, GS, BC, S, soft);
        }
        S.volume = S.members.reduce((v, e) => v + e.L*e.A, 0);
        log(`카탈로그 단면 적용: ${[...new Set(S.members.map(e => e.A))].length}종 사용, 체적 ${S.volume.toExponential(4)} mm³${over ? `, 카탈로그 최대보다 큰 부재 ${over}개 (연속값 유지)` : ''}`);
      }
      fem.results.forEach(r => log(`강성 해석 ${r.name}${soft.length ? ' (참고용: 메커니즘 노드에 약한 스프링)' : ''}: |u|max=${r.umax.toExponential(3)} mm, 응력 이용률 ${r.utilMax.toFixed(3)}, 좌굴 이용률 ${r.buckUtil.toFixed(3)}, PCG ${r.iters}회${r.converged ? '' : ' (미수렴)'}`));
    } else log('자유도가 커서 강성 해석을 생략했습니다');
    let cap = null;
    if (!stab.skipped) {
      cap = loadCapacity(highs, P, GS, BC, S, fem);
      cap.forEach(c => log(`최대 하중 배수 ${c.name} (좌굴 안전율 γ 미포함): 첫 부재 한계 λ=${c.elastic ? c.elastic.toFixed(3) : '–'}${c.elastic ? ` (${c.elasticMode}, 부재 ${c.elasticMember + 1})` : ''} · 붕괴(재분배 허용) λ=${c.plastic === null ? '–' : Number.isFinite(c.plastic) ? c.plastic.toFixed(3) : '∞'}${c.plasticAt && c.plasticAt.length ? ` (한계 도달 ${c.plasticAt.length}개 중 좌굴 ${c.plasticAt.filter(x => x[1] === '좌굴').length})` : ''}${stab.m > 0 ? ' — 메커니즘이 남아 있어 참고용' : ''}`));
    }
    tick('fem', tS);
    const sec = k => ((tm[k] || 0)/1000).toFixed(1);
    log(`단계별 시간 [s]: 부재 추가 ${sec('grow')} · 최적화(교차·단면) ${sec('opt')} · 기하 ${sec('geom')} · 좌굴·잠금 ${sec('stab')} · 검증(FSD·강성) ${sec('fem')}`);
    log(`완료: ${Date.now()-t0} ms`);
    return { ok: true, GS, BC, opt, S, stab, fem, crossLeft, added: added.length, buck, fold, fold0, fsd, symPlanes, grow, geom, nCrossSplit, cap, timing: tm, time: Date.now()-t0 };
  }

  return { solveRoundUp, loadCapacity, solveMinArea, symmetryMaps, solveLayout, candidatePool, growGroundStructure, geometryOptimize, splitCrossings, addNode, buildGroundStructure, buildBC, optimize, extractStructure, stabilityCheck, stabilizeGreedy, stiffnessAnalysis, toInp, run, symEig, segsCross, foldAnalysis, sectionI };
})();
if (typeof module !== 'undefined') module.exports = GSTO;
