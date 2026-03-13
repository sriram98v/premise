// ── Tab switching ──────────────────────────────────────────────────────────
document.querySelectorAll('.nav-tab').forEach(btn => {
  btn.addEventListener('click', () => {
    document.querySelectorAll('.nav-tab').forEach(b => {
      b.classList.remove('active');
      b.setAttribute('aria-selected', 'false');
    });
    document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
    btn.classList.add('active');
    btn.setAttribute('aria-selected', 'true');
    document.getElementById('panel-' + btn.dataset.tab).classList.add('active');
  });
});

// ── Shared helpers ─────────────────────────────────────────────────────────
function fmtSize(bytes) {
  if (bytes < 1024) return bytes + ' B';
  if (bytes < 1048576) return (bytes / 1024).toFixed(1) + ' KB';
  return (bytes / 1048576).toFixed(1) + ' MB';
}

function setStatus(tab, state /* 'running' | 'done' | 'error' | null */) {
  ['running', 'done', 'error'].forEach(s => {
    const el = document.getElementById(tab + '-status-' + s);
    if (el) el.classList.toggle('visible', s === state);
  });
}

// ── File validation helpers ─────────────────────────────────────────────────
function readFileHead(file, bytes) {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload  = e => resolve(e.target.result);
    reader.onerror = () => reject(new Error('Could not read file'));
    reader.readAsText(file.slice(0, bytes));
  });
}

function readFileBytesHead(file, bytes) {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload  = e => resolve(new Uint8Array(e.target.result));
    reader.onerror = () => reject(new Error('Could not read file'));
    reader.readAsArrayBuffer(file.slice(0, bytes));
  });
}

// Returns null if valid, or an error string if not.
async function validateFile(file, type) {
  const name = file.name.toLowerCase();

  if (type === 'fasta') {
    if (!name.endsWith('.fasta') && !name.endsWith('.fa'))
      return 'Expected a .fasta or .fa file.';
    const text = await readFileHead(file, 2048);
    if (!text.trimStart().startsWith('>'))
      return 'File does not look like a FASTA file (first non-whitespace character must be ">").';
    return null;
  }

  if (type === 'fastq') {
    const isGz = name.endsWith('.gz');
    if (!name.endsWith('.fastq') && !name.endsWith('.fq') && !isGz)
      return 'Expected a .fastq, .fq, or .fastq.gz file.';
    if (isGz) {
      const bytes = await readFileBytesHead(file, 2);
      if (bytes[0] !== 0x1f || bytes[1] !== 0x8b)
        return 'File has a .gz extension but does not appear to be a valid gzip file.';
    } else {
      const text = await readFileHead(file, 2048);
      if (!text.trimStart().startsWith('@'))
        return 'File does not look like a FASTQ file (first non-whitespace character must be "@").';
    }
    return null;
  }

  if (type === 'fmidx') {
    if (!name.endsWith('.fmidx'))
      return 'Expected a .fmidx index file (build one with the Build tab).';
    return null;
  }

  return null;
}

function setupDropzone(dropzoneId, inputId, chosenId, nameId, sizeId, onFile, errorId, validate) {
  const dz    = document.getElementById(dropzoneId);
  const input = document.getElementById(inputId);

  function setError(msg) {
    if (!errorId) return;
    const el = document.getElementById(errorId);
    if (!el) return;
    el.textContent = msg || '';
    el.classList.toggle('visible', !!msg);
  }

  async function showFile(file) {
    setError(null);
    if (validate) {
      const err = await validate(file);
      if (err) {
        setError(err);
        document.getElementById(chosenId).classList.remove('visible');
        return;
      }
    }
    document.getElementById(nameId).textContent = file.name;
    document.getElementById(sizeId).textContent = fmtSize(file.size);
    document.getElementById(chosenId).classList.add('visible');
    onFile(file);
  }

  input.addEventListener('change', () => { if (input.files[0]) showFile(input.files[0]); });
  dz.addEventListener('dragover',  e => { e.preventDefault(); dz.classList.add('dragover'); });
  dz.addEventListener('dragleave', () => dz.classList.remove('dragover'));
  dz.addEventListener('drop', e => {
    e.preventDefault();
    dz.classList.remove('dragover');
    const file = e.dataTransfer.files[0];
    if (file) showFile(file);
  });
}

function escHtml(s) {
  return String(s)
    .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

function ncbiLink(refId) {
  const url = 'https://www.ncbi.nlm.nih.gov/search/all/?term=' + encodeURIComponent(refId);
  return `<a href="${url}" target="_blank" rel="noopener noreferrer">${escHtml(refId)}</a>`;
}

// ── Build tab ──────────────────────────────────────────────────────────────
const buildBtn = document.getElementById('build-btn');

setupDropzone(
  'build-dropzone', 'build-fasta-input',
  'build-file-chosen', 'build-file-name', 'build-file-size',
  () => { buildBtn.disabled = false; },
  'build-fasta-error', f => validateFile(f, 'fasta')
);

buildBtn.addEventListener('click', async () => {
  const input = document.getElementById('build-fasta-input');
  const file  = input.files[0];
  if (!file) return;

  setStatus('build', 'running');
  buildBtn.disabled = true;

  try {
    const res = await fetch('/api/build', {
      method:  'POST',
      headers: { 'Content-Type': 'application/octet-stream' },
      body:    file,
    });

    if (!res.ok) {
      const msg = await res.text();
      throw new Error(msg || `Server error ${res.status}`);
    }

    const logRaw  = res.headers.get('X-Premise-Log') || '';
    const blob    = await res.blob();
    const stem    = file.name.replace(/\.[^.]+$/, '');
    const outName = stem + '.fmidx';
    const url     = URL.createObjectURL(blob);

    document.getElementById('build-result-filename').textContent = outName;
    const link    = document.getElementById('build-download-link');
    link.href     = url;
    link.download = outName;

    const buildLog = document.getElementById('build-log');
    if (logRaw) {
      buildLog.textContent = logRaw.replace(/ \| /g, '\n');
      buildLog.style.display = '';
    }

    setStatus('build', 'done');
  } catch (err) {
    document.getElementById('build-error-message').textContent = err.message;
    setStatus('build', 'error');
  } finally {
    buildBtn.disabled = false;
  }
});

// ── Align tab ──────────────────────────────────────────────────────────────
const alignBtn   = document.getElementById('align-btn');
const alignFiles = { index: null, r1: null, r2: null };
let   alignSession = null;

function checkAlignReady() {
  if (alignBtn) alignBtn.disabled = !(alignFiles.index && alignFiles.r1 && alignFiles.r2);
}

if (alignBtn) {
setupDropzone('align-index-dropzone', 'align-fmidx-input',
  'align-index-chosen', 'align-index-name', 'align-index-size',
  f => { alignFiles.index = f; checkAlignReady(); },
  'align-index-error', f => validateFile(f, 'fmidx'));

setupDropzone('align-r1-dropzone', 'align-r1-input',
  'align-r1-chosen', 'align-r1-name', 'align-r1-size',
  f => { alignFiles.r1 = f; checkAlignReady(); },
  'align-r1-error', f => validateFile(f, 'fastq'));

setupDropzone('align-r2-dropzone', 'align-r2-input',
  'align-r2-chosen', 'align-r2-name', 'align-r2-size',
  f => { alignFiles.r2 = f; checkAlignReady(); },
  'align-r2-error', f => validateFile(f, 'fastq'));
} // end if (alignBtn)

// Progress bar helpers
let _pRaf = null, _pStart = null, _pFrom = 0;

function setProgress(pct, label) {
  if (_pRaf) { cancelAnimationFrame(_pRaf); _pRaf = null; }
  const bar = document.getElementById('align-progress-bar');
  const lbl = document.getElementById('align-progress-label');
  if (bar) bar.style.width = pct + '%';
  if (lbl && label) lbl.textContent = label;
  _pFrom = pct;
}

function animateProgress(label, fromPct, targetPct = 92, halfTimeSec = 18) {
  if (_pRaf) { cancelAnimationFrame(_pRaf); _pRaf = null; }
  const bar = document.getElementById('align-progress-bar');
  const lbl = document.getElementById('align-progress-label');
  if (lbl && label) lbl.textContent = label;
  _pStart = Date.now();
  const range = targetPct - fromPct;
  function tick() {
    const t = (Date.now() - _pStart) / 1000;
    const p = fromPct + range * (1 - Math.exp(-t / halfTimeSec));
    if (bar) bar.style.width = p.toFixed(1) + '%';
    _pRaf = requestAnimationFrame(tick);
  }
  _pRaf = requestAnimationFrame(tick);
}

// Upload a single file to the server; returns the session ID
async function uploadAlignFile(part, file, sessionId) {
  const ext  = file.name.includes('.') ? file.name.split('.').pop() : 'fastq';
  const sess = sessionId || 'new';
  const res  = await fetch(
    `/api/align/upload?part=${part}&session=${sess}&ext=${encodeURIComponent(ext)}`,
    { method: 'POST', headers: { 'Content-Type': 'application/octet-stream' }, body: file }
  );
  if (!res.ok) throw new Error(`Upload failed (${part}): ${await res.text()}`);
  const data = await res.json();
  return data.session;
}

// ── Dark mode ──────────────────────────────────────────────────────────────
(function () {
  const btn  = document.getElementById('dark-toggle');
  const moon = document.getElementById('icon-moon');
  const sun  = document.getElementById('icon-sun');

  function applyDark(on) {
    document.documentElement.setAttribute('data-theme', on ? 'dark' : 'light');
    moon.style.display = on ? 'none' : '';
    sun.style.display  = on ? ''     : 'none';
  }

  applyDark(localStorage.getItem('dark') === '1');

  btn.addEventListener('click', () => {
    const next = document.documentElement.getAttribute('data-theme') !== 'dark';
    localStorage.setItem('dark', next ? '1' : '0');
    applyDark(next);
  });
})();

// ── Re-render charts on theme change ───────────────────────────────────────
new MutationObserver(() => {
  if (_convergenceData.length) renderConvergenceChart(_convergenceData);
  if (_lastPropsTsv) renderQueryPie(_lastPropsTsv);
}).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });

// ── Align table state ──────────────────────────────────────────────────────
let _tsvHeaders  = [];
let _tsvRows     = [];   // raw rows (never mutated)
let _filteredRows = [];  // after regex filter
let _sortedRows   = [];  // after sort  (source for pagination)
let _tPage    = 0;       // current page (0-based)
let _sortCol  = null;    // null | 0 | 1 | 2
let _sortDir  = 'asc';   // 'asc' | 'desc'
let _filterRx = [null, null]; // compiled RegExp or null for col 0 & 1

function tPerPage() {
  return parseInt(document.getElementById('align-rows-per-page').value, 10) || 50;
}

// ── Filter ─────────────────────────────────────────────────────────────────
function compileFilter(inputId, colIdx) {
  const input = document.getElementById(inputId);
  const val   = input.value.trim();
  if (!val) {
    input.classList.remove('invalid');
    _filterRx[colIdx] = null;
    return true;
  }
  try {
    _filterRx[colIdx] = new RegExp(val, 'i');
    input.classList.remove('invalid');
    return true;
  } catch {
    input.classList.add('invalid');
    _filterRx[colIdx] = null;
    return false;
  }
}

function applyFilters() {
  _filteredRows = _tsvRows.filter(row => {
    const cells = row.split('\t');
    for (let i = 0; i < 2; i++) {
      if (_filterRx[i] && !_filterRx[i].test(cells[i] || '')) return false;
    }
    return true;
  });
}

// ── Sort ───────────────────────────────────────────────────────────────────
function applySort() {
  if (_sortCol === null) { _sortedRows = _filteredRows.slice(); return; }
  const col = _sortCol, dir = _sortDir;
  _sortedRows = _filteredRows.slice().sort((a, b) => {
    const ca = a.split('\t')[col] ?? '';
    const cb = b.split('\t')[col] ?? '';
    let cmp;
    if (col === 2) {
      // numeric; treat '-' (unclassified) as always-last
      if (ca === '-' && cb === '-') return 0;
      if (ca === '-') return 1;
      if (cb === '-') return -1;
      cmp = parseFloat(ca) - parseFloat(cb);
    } else {
      cmp = ca.localeCompare(cb, undefined, { sensitivity: 'base' });
    }
    return dir === 'asc' ? cmp : -cmp;
  });
}

// ── Full recompute (filter → sort → page 0 → render) ──────────────────────
function recompute() {
  applyFilters();
  applySort();
  _tPage = 0;
  renderTablePage();
}

// ── Init (called once after TSV arrives) ───────────────────────────────────
function initAlignTable(tsv) {
  const lines = tsv.trim().split('\n');
  if (lines.length < 2) { _tsvHeaders = []; _tsvRows = []; return; }
  _tsvHeaders   = lines[0].split('\t');
  _tsvRows      = lines.slice(1);
  _filteredRows = _tsvRows.slice();
  _sortedRows   = _tsvRows.slice();
  _sortCol = null;
  _sortDir = 'asc';
  _filterRx = [null, null];
  document.getElementById('filter-read-id').value = '';
  document.getElementById('filter-ref-id').value  = '';
  document.getElementById('filter-read-id').classList.remove('invalid');
  document.getElementById('filter-ref-id').classList.remove('invalid');
  _tPage = 0;
  renderTablePage();
}

// ── Render the current page ────────────────────────────────────────────────
function renderTablePage() {
  const perPage   = tPerPage();
  const total     = _sortedRows.length;
  const rawTotal  = _tsvRows.length;
  const pageCount = Math.max(1, Math.ceil(total / perPage));
  _tPage = Math.min(_tPage, pageCount - 1);

  const start = _tPage * perPage;
  const end   = Math.min(start + perPage, total);
  const slice = _sortedRows.slice(start, end);

  // Info line
  const filtered = total < rawTotal
    ? ` (${total.toLocaleString()} matching)`
    : '';
  document.getElementById('align-table-info').innerHTML =
    `<strong>${total > 0 ? (start + 1).toLocaleString() : 0}–${end.toLocaleString()}</strong>`
    + ` of <strong>${rawTotal.toLocaleString()}</strong> rows${filtered}`;

  // Sort-icon helper
  const ICON_NONE = '<span class="sort-icon">↕</span>';
  const ICON_ASC  = '<span class="sort-icon">↑</span>';
  const ICON_DESC = '<span class="sort-icon">↓</span>';

  // Table
  let html = '<thead><tr>';
  _tsvHeaders.forEach((h, i) => {
    const active = _sortCol === i;
    const icon   = active ? (_sortDir === 'asc' ? ICON_ASC : ICON_DESC) : ICON_NONE;
    html += `<th class="th-sortable${active ? ' th-sort-active' : ''}" data-col="${i}">`
          + escHtml(h) + icon + '</th>';
  });
  html += '</tr></thead><tbody>';

  for (const row of slice) {
    const cells = row.split('\t');
    const isUnclassified = cells[1] === 'unclassified';
    html += `<tr${isUnclassified ? ' class="unclassified"' : ''}>`
      + cells.map((c, i) => `<td>${i === 1 && c !== 'unclassified' ? ncbiLink(c) : escHtml(c)}</td>`).join('') + '</tr>';
  }
  html += '</tbody>';
  document.getElementById('align-table').innerHTML = html;

  // Attach sort click handlers
  document.getElementById('align-table').querySelectorAll('th[data-col]').forEach(th => {
    th.addEventListener('click', () => {
      const col = parseInt(th.dataset.col, 10);
      if (_sortCol === col) {
        _sortDir = _sortDir === 'asc' ? 'desc' : 'asc';
      } else {
        _sortCol = col;
        _sortDir = 'asc';
      }
      applySort();
      _tPage = 0;
      renderTablePage();
    });
  });

  renderPagination(pageCount);
}

// ── Pagination ─────────────────────────────────────────────────────────────
function renderPagination(pageCount) {
  const pg  = _tPage;
  const bar = document.getElementById('align-pagination');
  if (pageCount <= 1) { bar.innerHTML = ''; return; }

  const pages = new Set([0, pageCount - 1]);
  for (let i = Math.max(0, pg - 2); i <= Math.min(pageCount - 1, pg + 2); i++) pages.add(i);
  const sorted = [...pages].sort((a, b) => a - b);

  let html = `<button class="pg-btn" id="pg-prev" ${pg === 0 ? 'disabled' : ''}
                title="Previous page">&#8249;</button>`;
  let prev = -1;
  for (const p of sorted) {
    if (prev !== -1 && p - prev > 1) html += `<span class="pg-ellipsis">…</span>`;
    html += `<button class="pg-btn${p === pg ? ' active' : ''}" data-page="${p}">${(p + 1).toLocaleString()}</button>`;
    prev = p;
  }
  html += `<button class="pg-btn" id="pg-next" ${pg === pageCount - 1 ? 'disabled' : ''}
             title="Next page">&#8250;</button>`;

  bar.innerHTML = html;
  bar.querySelector('#pg-prev')?.addEventListener('click', () => { _tPage--; renderTablePage(); });
  bar.querySelector('#pg-next')?.addEventListener('click', () => { _tPage++; renderTablePage(); });
  bar.querySelectorAll('[data-page]').forEach(btn =>
    btn.addEventListener('click', () => { _tPage = parseInt(btn.dataset.page, 10); renderTablePage(); })
  );
}

// ── Wire up filter inputs + rows-per-page ──────────────────────────────────
let _filterTimer = null;
function onFilterInput() {
  clearTimeout(_filterTimer);
  _filterTimer = setTimeout(() => {
    compileFilter('filter-read-id', 0);
    compileFilter('filter-ref-id',  1);
    recompute();
  }, 250);  // debounce 250 ms
}

if (alignBtn) {
  document.getElementById('filter-read-id').addEventListener('input', onFilterInput);
  document.getElementById('filter-ref-id').addEventListener('input', onFilterInput);
  document.getElementById('align-rows-per-page')
    .addEventListener('change', () => { _tPage = 0; renderTablePage(); });
}

// ── Query tab ──────────────────────────────────────────────────────────────
const queryFiles   = { index: null, r1: null, r2: null };
let   querySession = null;

function checkQueryReady() {
  document.getElementById('query-btn').disabled =
    !(queryFiles.index && queryFiles.r1 && queryFiles.r2);
}

setupDropzone('query-index-dropzone', 'query-fmidx-input',
  'query-index-chosen', 'query-index-name', 'query-index-size',
  f => { queryFiles.index = f; checkQueryReady(); },
  'query-index-error', f => validateFile(f, 'fmidx'));

setupDropzone('query-r1-dropzone', 'query-r1-input',
  'query-r1-chosen', 'query-r1-name', 'query-r1-size',
  f => { queryFiles.r1 = f; checkQueryReady(); },
  'query-r1-error', f => validateFile(f, 'fastq'));

setupDropzone('query-r2-dropzone', 'query-r2-input',
  'query-r2-chosen', 'query-r2-name', 'query-r2-size',
  f => { queryFiles.r2 = f; checkQueryReady(); },
  'query-r2-error', f => validateFile(f, 'fastq'));

// Query progress helpers
let _qRaf = null;

function setQueryProgress(pct, label) {
  if (_qRaf) { cancelAnimationFrame(_qRaf); _qRaf = null; }
  const bar = document.getElementById('query-progress-bar');
  const lbl = document.getElementById('query-progress-label');
  if (bar) bar.style.width = pct + '%';
  if (lbl && label) lbl.textContent = label;
}

function animateQueryProgress(label, fromPct, targetPct = 92, halfTimeSec = 18) {
  if (_qRaf) { cancelAnimationFrame(_qRaf); _qRaf = null; }
  const bar   = document.getElementById('query-progress-bar');
  const lbl   = document.getElementById('query-progress-label');
  if (lbl && label) lbl.textContent = label;
  const start = Date.now();
  const range = targetPct - fromPct;
  function tick() {
    const t = (Date.now() - start) / 1000;
    const p = fromPct + range * (1 - Math.exp(-t / halfTimeSec));
    if (bar) bar.style.width = p.toFixed(1) + '%';
    _qRaf = requestAnimationFrame(tick);
  }
  _qRaf = requestAnimationFrame(tick);
}

// Upload a file to /api/query/upload
async function uploadQueryFile(part, file, sessionId) {
  const ext  = file.name.includes('.') ? file.name.split('.').pop() : 'fastq';
  const sess = sessionId || 'new';
  const res  = await fetch(
    `/api/query/upload?part=${part}&session=${sess}&ext=${encodeURIComponent(ext)}`,
    { method: 'POST', headers: { 'Content-Type': 'application/octet-stream' }, body: file }
  );
  if (!res.ok) throw new Error(`Upload failed (${part}): ${await res.text()}`);
  const data = await res.json();
  return data.session;
}

// ── Query table state ───────────────────────────────────────────────────────
let _qTsvHeaders  = [];
let _qTsvRows     = [];
let _qFilteredRows = [];
let _qSortedRows   = [];
let _qPage    = 0;
let _qSortCol = null;
let _qSortDir = 'asc';
let _qFilterRx = [null, null];
let _qMatchesUrl = null, _qPosteriorsUrl = null, _qPropsUrl = null, _qAlignsUrl = null;
let _lastPropsTsv = null;

function qTPerPage() {
  return parseInt(document.getElementById('query-rows-per-page').value, 10) || 50;
}

function qCompileFilter(inputId, colIdx) {
  const input = document.getElementById(inputId);
  const val   = input.value.trim();
  if (!val) {
    input.classList.remove('invalid');
    _qFilterRx[colIdx] = null;
    return true;
  }
  try {
    _qFilterRx[colIdx] = new RegExp(val, 'i');
    input.classList.remove('invalid');
    return true;
  } catch {
    input.classList.add('invalid');
    _qFilterRx[colIdx] = null;
    return false;
  }
}

function qApplyFilters() {
  _qFilteredRows = _qTsvRows.filter(row => {
    const cells = row.split('\t');
    for (let i = 0; i < 2; i++) {
      if (_qFilterRx[i] && !_qFilterRx[i].test(cells[i] || '')) return false;
    }
    return true;
  });
}

function qApplySort() {
  if (_qSortCol === null) { _qSortedRows = _qFilteredRows.slice(); return; }
  const col = _qSortCol, dir = _qSortDir;
  _qSortedRows = _qFilteredRows.slice().sort((a, b) => {
    const ca = a.split('\t')[col] ?? '';
    const cb = b.split('\t')[col] ?? '';
    let cmp;
    if (col === 2) {
      if (ca === '-' && cb === '-') return 0;
      if (ca === '-') return 1;
      if (cb === '-') return -1;
      cmp = parseFloat(ca) - parseFloat(cb);
    } else {
      cmp = ca.localeCompare(cb, undefined, { sensitivity: 'base' });
    }
    return dir === 'asc' ? cmp : -cmp;
  });
}

function qRecompute() {
  qApplyFilters();
  qApplySort();
  _qPage = 0;
  qRenderTablePage();
}

function initQueryTable(tsv) {
  const lines = tsv.trim().split('\n');
  if (lines.length < 2) { _qTsvRows = []; _qFilteredRows = []; _qSortedRows = []; return; }
  _qTsvHeaders  = lines[0].split('\t');
  _qTsvRows     = lines.slice(1);
  _qFilteredRows = _qTsvRows.slice();
  _qSortedRows   = _qTsvRows.slice();
  _qSortCol = null;
  _qSortDir = 'asc';
  _qFilterRx = [null, null];
  document.getElementById('query-filter-read-id').value = '';
  document.getElementById('query-filter-ref-id').value  = '';
  document.getElementById('query-filter-read-id').classList.remove('invalid');
  document.getElementById('query-filter-ref-id').classList.remove('invalid');
  _qPage = 0;
  qRenderTablePage();
}

function qRenderTablePage() {
  const perPage   = qTPerPage();
  const total     = _qSortedRows.length;
  const rawTotal  = _qTsvRows.length;
  const pageCount = Math.max(1, Math.ceil(total / perPage));
  _qPage = Math.min(_qPage, pageCount - 1);

  const start = _qPage * perPage;
  const end   = Math.min(start + perPage, total);
  const slice = _qSortedRows.slice(start, end);

  const filtered = total < rawTotal
    ? ` (${total.toLocaleString()} matching)` : '';
  document.getElementById('query-table-info').innerHTML =
    `<strong>${total > 0 ? (start + 1).toLocaleString() : 0}–${end.toLocaleString()}</strong>`
    + ` of <strong>${rawTotal.toLocaleString()}</strong> rows${filtered}`;

  const ICON_NONE = '<span class="sort-icon">↕</span>';
  const ICON_ASC  = '<span class="sort-icon">↑</span>';
  const ICON_DESC = '<span class="sort-icon">↓</span>';

  let html = '<thead><tr>';
  _qTsvHeaders.forEach((h, i) => {
    const active = _qSortCol === i;
    const icon   = active ? (_qSortDir === 'asc' ? ICON_ASC : ICON_DESC) : ICON_NONE;
    html += `<th class="th-sortable${active ? ' th-sort-active' : ''}" data-col="${i}">`
          + escHtml(h) + icon + '</th>';
  });
  html += '</tr></thead><tbody>';

  for (const row of slice) {
    const cells = row.split('\t');
    const isUnclassified = cells[1] === 'unclassified';
    html += `<tr${isUnclassified ? ' class="unclassified"' : ''}>`
      + cells.map((c, i) => `<td>${i === 1 && c !== 'unclassified' ? ncbiLink(c) : escHtml(c)}</td>`).join('') + '</tr>';
  }
  html += '</tbody>';
  document.getElementById('query-table').innerHTML = html;

  document.getElementById('query-table').querySelectorAll('th[data-col]').forEach(th => {
    th.addEventListener('click', () => {
      const col = parseInt(th.dataset.col, 10);
      if (_qSortCol === col) {
        _qSortDir = _qSortDir === 'asc' ? 'desc' : 'asc';
      } else {
        _qSortCol = col;
        _qSortDir = 'asc';
      }
      qApplySort();
      _qPage = 0;
      qRenderTablePage();
    });
  });

  qRenderPagination(pageCount);
}

function qRenderPagination(pageCount) {
  const pg   = _qPage;
  const bars = [
    document.getElementById('query-pagination-top'),
    document.getElementById('query-pagination'),
  ].filter(Boolean);

  if (pageCount <= 1) { bars.forEach(b => b.innerHTML = ''); return; }

  const pages = new Set([0, pageCount - 1]);
  for (let i = Math.max(0, pg - 2); i <= Math.min(pageCount - 1, pg + 2); i++) pages.add(i);
  const sorted = [...pages].sort((a, b) => a - b);

  let html = `<button class="pg-btn" data-qpg-dir="prev" ${pg === 0 ? 'disabled' : ''}
                title="Previous page">&#8249;</button>`;
  let prev = -1;
  for (const p of sorted) {
    if (prev !== -1 && p - prev > 1) html += `<span class="pg-ellipsis">…</span>`;
    html += `<button class="pg-btn${p === pg ? ' active' : ''}" data-qpage="${p}">${(p + 1).toLocaleString()}</button>`;
    prev = p;
  }
  html += `<button class="pg-btn" data-qpg-dir="next" ${pg === pageCount - 1 ? 'disabled' : ''}
             title="Next page">&#8250;</button>`;

  bars.forEach(bar => {
    bar.innerHTML = html;
    bar.querySelector('[data-qpg-dir="prev"]')?.addEventListener('click', () => { _qPage--; qRenderTablePage(); });
    bar.querySelector('[data-qpg-dir="next"]')?.addEventListener('click', () => { _qPage++; qRenderTablePage(); });
    bar.querySelectorAll('[data-qpage]').forEach(btn =>
      btn.addEventListener('click', () => { _qPage = parseInt(btn.dataset.qpage, 10); qRenderTablePage(); })
    );
  });
}

// Filter + rows-per-page wiring
let _qFilterTimer = null;
function onQueryFilterInput() {
  clearTimeout(_qFilterTimer);
  _qFilterTimer = setTimeout(() => {
    qCompileFilter('query-filter-read-id', 0);
    qCompileFilter('query-filter-ref-id',  1);
    qRecompute();
  }, 250);
}
document.getElementById('query-filter-read-id').addEventListener('input', onQueryFilterInput);
document.getElementById('query-filter-ref-id').addEventListener('input', onQueryFilterInput);
document.getElementById('query-rows-per-page')
  .addEventListener('change', () => { _qPage = 0; qRenderTablePage(); });

// ── D3 Pie chart ───────────────────────────────────────────────────────────
const PIE_COLORS = [
  '#0d6e6e','#e97b2a','#2a7be9','#e92a6c','#7be92a',
  '#9b59b6','#e9c02a','#2ae9c0','#e94a2a','#2ac0e9',
  '#8e8e00','#e92ab6'
];

function renderQueryPie(tsv) {
  const wrap   = document.getElementById('query-pie-wrap');
  const svgEl  = document.getElementById('query-pie');
  const legend = document.getElementById('query-pie-legend');
  svgEl.innerHTML    = '';
  legend.innerHTML   = '';

  const lines = tsv.trim().split('\n').filter(l => l.trim());
  if (!lines.length) return;

  let data = lines.map(l => {
    const p = l.split('\t');
    return { id: p[0] || '', proportion: parseFloat(p[1]) || 0 };
  }).filter(d => d.proportion > 0);

  if (!data.length) return;

  // Normalize and sort
  const total = data.reduce((s, d) => s + d.proportion, 0);
  data.forEach(d => { d.proportion /= total; });
  data.sort((a, b) => b.proportion - a.proportion);

  const W = 260, H = 260, R = 110;

  // Resolve CSS vars for PNG-safe colours
  const cs    = getComputedStyle(document.documentElement);
  const cSurf = cs.getPropertyValue('--pico-card-background-color').trim() || '#ffffff';
  const cTxt  = cs.getPropertyValue('--pico-color').trim()                 || '#1a1d23';

  svgEl.setAttribute('width',   W);
  svgEl.setAttribute('height',  H);
  svgEl.setAttribute('viewBox', `0 0 ${W} ${H}`);

  const svg = d3.select(svgEl);

  // Background rect (needed for PNG export)
  svg.append('rect').attr('width', W).attr('height', H).attr('fill', cSurf);

  const g = svg.append('g').attr('transform', `translate(${W/2},${H/2})`);

  const pie  = d3.pie().value(d => d.proportion).sort(null);
  const arc  = d3.arc().innerRadius(0).outerRadius(R);
  const arcH = d3.arc().innerRadius(0).outerRadius(R + 8); // hover expand

  // Tooltip div
  const tip = document.getElementById('pie-tooltip');

  const slices = g.selectAll('path')
    .data(pie(data))
    .enter()
    .append('path')
      .attr('d', arc)
      .attr('fill', (d, i) => PIE_COLORS[i % PIE_COLORS.length])
      .attr('stroke', cSurf)
      .attr('stroke-width', 1.5)
      .style('cursor', 'pointer')
      .style('transition', 'd 0.15s');

  slices
    .on('mouseover', function(event, d) {
      d3.select(this).attr('d', arcH);
      const pct = (d.data.proportion * 100).toFixed(2);
      tip.innerHTML = `<strong>${escHtml(d.data.id)}</strong><br>${pct}%`;
      tip.style.display = 'block';
    })
    .on('mousemove', function(event) {
      const wrapRect = wrap.getBoundingClientRect();
      let tx = event.clientX - wrapRect.left + 14;
      let ty = event.clientY - wrapRect.top  - 48;
      if (tx + 200 > wrapRect.width) tx = event.clientX - wrapRect.left - 206;
      if (ty < 0) ty = 4;
      tip.style.left = tx + 'px';
      tip.style.top  = ty + 'px';
    })
    .on('mouseleave', function() {
      d3.select(this).attr('d', arc);
      tip.style.display = 'none';
    });

  // Legend
  data.forEach((d, i) => {
    const color = PIE_COLORS[i % PIE_COLORS.length];
    const pct   = (d.proportion * 100).toFixed(2);
    const li    = document.createElement('li');
    li.innerHTML =
      `<span class="pie-legend-dot" style="background:${color}"></span>` +
      `<span>${ncbiLink(d.id)} — ${pct}%</span>`;
    legend.appendChild(li);
  });
}

// ── Pie chart PNG download ─────────────────────────────────────────────────
function downloadPiePng() {
  const svgEl = document.getElementById('query-pie');
  if (!svgEl || !svgEl.getAttribute('width')) return;

  // Hide tooltip before serialising
  const tip = document.getElementById('pie-tooltip');
  const tipWas = tip.style.display;
  tip.style.display = 'none';

  const W = parseInt(svgEl.getAttribute('width'));
  const H = parseInt(svgEl.getAttribute('height'));
  const svgData = new XMLSerializer().serializeToString(svgEl);

  tip.style.display = tipWas;

  const blob = new Blob([svgData], { type: 'image/svg+xml;charset=utf-8' });
  const url  = URL.createObjectURL(blob);
  const img  = new Image();
  img.onerror = () => URL.revokeObjectURL(url);
  img.onload  = () => {
    const canvas = document.createElement('canvas');
    const dpr = window.devicePixelRatio || 2;
    canvas.width  = W * dpr;
    canvas.height = H * dpr;
    const ctx = canvas.getContext('2d');
    ctx.scale(dpr, dpr);
    ctx.drawImage(img, 0, 0);
    URL.revokeObjectURL(url);
    const a = document.createElement('a');
    a.download = 'abundances.png';
    a.href = canvas.toDataURL('image/png');
    a.click();
  };
  img.src = url;
}

// ── Sub-tab switching (query results) ──────────────────────────────────────
document.querySelectorAll('.sub-tab-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    document.querySelectorAll('.sub-tab-btn').forEach(b => b.classList.remove('active'));
    btn.classList.add('active');
    const tab = btn.dataset.subtab;
    document.getElementById('query-assignments-panel').style.display =
      tab === 'assignments' ? '' : 'none';
    document.getElementById('query-abundances-panel').style.display =
      tab === 'abundances' ? '' : 'none';
    document.getElementById('query-convergence-panel').style.display =
      tab === 'convergence' ? '' : 'none';
  });
});

// ── D3 EM Convergence chart ────────────────────────────────────────────────
let _convergenceData = [];

function renderConvergenceChart(values) {
  _convergenceData = values || [];
  const svgEl = document.getElementById('query-convergence-chart');
  svgEl.innerHTML = '';
  if (!_convergenceData.length) return;

  const W = 560, H = 310;
  const M = { top: 24, right: 28, bottom: 54, left: 84 };
  const cW = W - M.left - M.right;
  const cH = H - M.top  - M.bottom;

  svgEl.setAttribute('width',   W);
  svgEl.setAttribute('height',  H);
  svgEl.setAttribute('viewBox', `0 0 ${W} ${H}`);

  // Resolve CSS vars to concrete colours for PNG export
  const cs    = getComputedStyle(document.documentElement);
  const cAcc  = cs.getPropertyValue('--pico-primary').trim()                  || '#0d6e6e';
  const cSurf = cs.getPropertyValue('--pico-card-background-color').trim()    || '#ffffff';
  const cBord = cs.getPropertyValue('--pico-form-element-border-color').trim()|| '#dde1e8';
  const cMut  = cs.getPropertyValue('--pico-muted-color').trim()              || '#6b7280';
  const cTxt  = cs.getPropertyValue('--pico-color').trim()                    || '#1a1d23';

  const n    = _convergenceData.length;
  const minV = d3.min(_convergenceData);
  const maxV = d3.max(_convergenceData);

  const xScale = d3.scaleLinear().domain([0, n - 1]).range([0, cW]);
  const yScale = d3.scaleLinear().domain([minV, maxV]).nice().range([cH, 0]);

  const svg = d3.select(svgEl);

  // Background
  svg.append('rect').attr('width', W).attr('height', H).attr('fill', cSurf);

  const g = svg.append('g').attr('transform', `translate(${M.left},${M.top})`);

  // Grid lines (5 ticks on Y)
  const yTicks = yScale.ticks(5);
  g.selectAll('.grid-line')
    .data(yTicks)
    .enter().append('line')
      .attr('x1', 0).attr('x2', cW)
      .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
      .attr('stroke', cBord)
      .attr('stroke-width', d => d === yScale.domain()[0] ? 1.5 : 1)
      .attr('stroke-dasharray', d => d === yScale.domain()[0] ? null : '3,4');

  // Y axis
  const yAxis = d3.axisLeft(yScale)
    .ticks(5)
    .tickFormat(d => d.toExponential(2));

  g.append('g')
    .call(yAxis)
    .call(ax => ax.select('.domain').attr('stroke', cBord))
    .call(ax => ax.selectAll('.tick line').attr('stroke', cBord))
    .call(ax => ax.selectAll('.tick text')
      .attr('fill', cMut)
      .attr('font-family', 'ui-monospace,monospace')
      .attr('font-size', '10'));

  // X axis
  const xAxis = d3.axisBottom(xScale)
    .ticks(Math.min(n, 9))
    .tickFormat(d3.format('d'));

  g.append('g')
    .attr('transform', `translate(0,${cH})`)
    .call(xAxis)
    .call(ax => ax.select('.domain').attr('stroke', cBord))
    .call(ax => ax.selectAll('.tick line').attr('stroke', cBord))
    .call(ax => ax.selectAll('.tick text')
      .attr('fill', cMut)
      .attr('font-family', 'ui-monospace,monospace')
      .attr('font-size', '10'));

  // Axis labels
  g.append('text')
    .attr('x', cW / 2).attr('y', cH + 42)
    .attr('text-anchor', 'middle')
    .attr('font-size', '12').attr('fill', cTxt).attr('font-family', 'sans-serif')
    .text('EM Iteration');

  g.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', -(cH / 2)).attr('y', -M.left + 16)
    .attr('text-anchor', 'middle')
    .attr('font-size', '12').attr('fill', cTxt).attr('font-family', 'sans-serif')
    .text('Data Likelihood');

  // Area fill
  if (n > 1) {
    const area = d3.area()
      .x((d, i) => xScale(i))
      .y0(cH)
      .y1(d => yScale(d))
      .curve(d3.curveLinear);

    g.append('path')
      .datum(_convergenceData)
      .attr('fill', cAcc)
      .attr('opacity', 0.1)
      .attr('d', area);
  }

  // Line
  const line = d3.line()
    .x((d, i) => xScale(i))
    .y(d => yScale(d))
    .curve(d3.curveLinear);

  g.append('path')
    .datum(_convergenceData)
    .attr('fill', 'none')
    .attr('stroke', cAcc)
    .attr('stroke-width', 2)
    .attr('stroke-linejoin', 'round')
    .attr('stroke-linecap', 'round')
    .attr('d', line);

  // Interactive overlay elements
  const crosshair = g.append('line')
    .attr('id', 'convergence-crosshair')
    .attr('y1', 0).attr('y2', cH)
    .attr('stroke', cMut)
    .attr('stroke-width', 1)
    .attr('stroke-dasharray', '4,3')
    .style('display', 'none');

  const dot = g.append('circle')
    .attr('id', 'convergence-dot')
    .attr('r', 4.5)
    .attr('fill', cAcc)
    .attr('stroke', cSurf)
    .attr('stroke-width', 2)
    .style('display', 'none');

  // Transparent overlay rect for mouse events
  const overlay = g.append('rect')
    .attr('width', cW).attr('height', cH)
    .attr('fill', 'transparent')
    .style('cursor', 'crosshair');

  const tip = document.getElementById('convergence-tooltip');

  overlay.on('mousemove', function(event) {
    const [mx] = d3.pointer(event, this);
    let idx = n <= 1 ? 0 : Math.round((mx / cW) * (n - 1));
    idx = Math.max(0, Math.min(n - 1, idx));

    const cx = xScale(idx), cy = yScale(_convergenceData[idx]);
    crosshair.attr('x1', cx).attr('x2', cx).style('display', null);
    dot.attr('cx', cx).attr('cy', cy).style('display', null);

    tip.style.display = 'block';
    tip.innerHTML = `<strong>Iteration\u00a0${idx}</strong><br>${_convergenceData[idx].toExponential(6)}`;

    const wrap = document.getElementById('convergence-chart-wrap').getBoundingClientRect();
    let tx = event.clientX - wrap.left + 14;
    let ty = event.clientY - wrap.top  - 48;
    if (tx + 190 > wrap.width) tx = event.clientX - wrap.left - 196;
    if (ty < 0) ty = 4;
    tip.style.left = tx + 'px';
    tip.style.top  = ty + 'px';
  });

  overlay.on('mouseleave', () => {
    crosshair.style('display', 'none');
    dot.style('display', 'none');
    tip.style.display = 'none';
  });
}

function downloadConvergencePng() {
  const svgEl = document.getElementById('query-convergence-chart');
  if (!svgEl || !svgEl.getAttribute('width')) return;

  // Hide interactive overlays before serializing
  const cross = document.getElementById('convergence-crosshair');
  const d     = document.getElementById('convergence-dot');
  if (cross) cross.style.display = 'none';
  if (d)     d.style.display     = 'none';

  const W = parseInt(svgEl.getAttribute('width'));
  const H = parseInt(svgEl.getAttribute('height'));
  const svgData = new XMLSerializer().serializeToString(svgEl);

  // Restore immediately (serialization is sync)
  if (cross) cross.style.display = '';
  if (d)     d.style.display     = '';

  const blob = new Blob([svgData], { type: 'image/svg+xml;charset=utf-8' });
  const url  = URL.createObjectURL(blob);
  const img  = new Image();
  img.onerror = () => URL.revokeObjectURL(url);
  img.onload  = () => {
    const canvas = document.createElement('canvas');
    const dpr = window.devicePixelRatio || 2;
    canvas.width  = W * dpr;
    canvas.height = H * dpr;
    const ctx = canvas.getContext('2d');
    ctx.scale(dpr, dpr);
    ctx.drawImage(img, 0, 0);
    URL.revokeObjectURL(url);
    const a = document.createElement('a');
    a.download = 'em_convergence.png';
    a.href = canvas.toDataURL('image/png');
    a.click();
  };
  img.src = url;
}

// ── Run Query ──────────────────────────────────────────────────────────────
document.getElementById('query-btn').addEventListener('click', async () => {
  window.scrollTo({ top: 0, behavior: 'smooth' });
  setStatus('query', 'running');
  document.getElementById('query-results-section').style.display = 'none';
  document.getElementById('query-btn').disabled = true;

  try {
    setQueryProgress(0, 'Uploading index…');
    querySession = await uploadQueryFile('index', queryFiles.index, 'new');

    setQueryProgress(10, 'Uploading R1 reads…');
    querySession = await uploadQueryFile('r1', queryFiles.r1, querySession);

    setQueryProgress(20, 'Uploading R2 reads…');
    querySession = await uploadQueryFile('r2', queryFiles.r2, querySession);

    animateQueryProgress('Running EM classification…', 20);

    const mismatch  = document.getElementById('query-mismatch').value  || '5';
    const eps1      = document.getElementById('query-eps1').value      || '1e-32';
    const iter      = document.getElementById('query-iter').value      || '100';
    const threads   = document.getElementById('query-threads').value   || '0';
    const rho       = document.getElementById('query-rho').value       || '20';
    const omega     = document.getElementById('query-omega').value     || '1e-20';
    const eps2        = document.getElementById('query-eps2').value          || '1e-18';
    const emThreshold = document.getElementById('query-em-threshold').value  || '1e-6';
    const noPenalty   = document.getElementById('query-no-penalty').checked;

    const res = await fetch(
      `/api/query/run?session=${querySession}&mismatch=${mismatch}&eps_1=${eps1}` +
      `&iter=${iter}&threads=${threads}&rho=${rho}&omega=${omega}&eps_2=${eps2}` +
      `&em_threshold=${emThreshold}&no_penalty=${noPenalty}`,
      { method: 'POST' }
    );
    if (!res.ok) throw new Error(await res.text());
    const data = await res.json();

    if (_qRaf) { cancelAnimationFrame(_qRaf); _qRaf = null; }
    setQueryProgress(100, 'Done!');

    // Blob URLs
    if (_qMatchesUrl)    URL.revokeObjectURL(_qMatchesUrl);
    if (_qPosteriorsUrl) URL.revokeObjectURL(_qPosteriorsUrl);
    if (_qPropsUrl)      URL.revokeObjectURL(_qPropsUrl);
    if (_qAlignsUrl)     URL.revokeObjectURL(_qAlignsUrl);
    _qMatchesUrl    = URL.createObjectURL(new Blob([data.matches],    { type: 'text/tab-separated-values' }));
    _qPosteriorsUrl = URL.createObjectURL(new Blob([data.posteriors], { type: 'text/tab-separated-values' }));
    _qPropsUrl      = URL.createObjectURL(new Blob([data.props],      { type: 'text/tab-separated-values' }));
    _qAlignsUrl     = URL.createObjectURL(new Blob([data.aligns],     { type: 'text/tab-separated-values' }));

    document.getElementById('query-dl-matches').href    = _qMatchesUrl;
    document.getElementById('query-dl-posteriors').href = _qPosteriorsUrl;
    document.getElementById('query-dl-props').href      = _qPropsUrl;
    document.getElementById('query-dl-aligns').href     = _qAlignsUrl;

    const rowCount = data.matches.trim().split('\n').length - 1;
    document.getElementById('query-stats').textContent =
      `${rowCount.toLocaleString()} assignment row${rowCount !== 1 ? 's' : ''}`;

    initQueryTable(data.matches);
    _lastPropsTsv = data.props;
    renderQueryPie(data.props);
    renderConvergenceChart(data.convergence || []);

    const queryLog = document.getElementById('query-log');
    if (data.log) {
      queryLog.textContent = data.log;
      queryLog.style.display = '';
    }

    document.getElementById('query-results-section').style.display = '';
    const ph = document.getElementById('query-results-placeholder');
    if (ph) ph.style.display = 'none';
    // Reset to Assignments sub-tab
    document.querySelectorAll('.sub-tab-btn').forEach(b => b.classList.remove('active'));
    document.querySelector('.sub-tab-btn[data-subtab="assignments"]').classList.add('active');
    document.getElementById('query-assignments-panel').style.display  = '';
    document.getElementById('query-abundances-panel').style.display   = 'none';
    document.getElementById('query-convergence-panel').style.display  = 'none';

    setStatus('query', 'done');
  } catch (err) {
    document.getElementById('query-error-message').textContent = err.message;
    setStatus('query', 'error');
  } finally {
    document.getElementById('query-btn').disabled = false;
    if (_qRaf) { cancelAnimationFrame(_qRaf); _qRaf = null; }
  }
});

if (alignBtn) alignBtn.addEventListener('click', async () => {
  setStatus('align', 'running');
  document.getElementById('align-table-section').style.display = 'none';
  alignBtn.disabled = true;

  try {
    // ── Upload phase ──────────────────────────────────────────────
    setProgress(0, 'Uploading index…');
    alignSession = await uploadAlignFile('index', alignFiles.index, 'new');

    setProgress(15, 'Uploading R1 reads…');
    alignSession = await uploadAlignFile('r1', alignFiles.r1, alignSession);

    setProgress(30, 'Uploading R2 reads…');
    alignSession = await uploadAlignFile('r2', alignFiles.r2, alignSession);

    // ── Alignment phase ───────────────────────────────────────────
    animateProgress('Running alignment…', 30);

    const mismatch = document.getElementById('align-mismatch').value || '5';
    const threads  = document.getElementById('align-threads').value  || '0';
    const res = await fetch(
      `/api/align/run?session=${alignSession}&mismatch=${mismatch}&threads=${threads}`,
      { method: 'POST' }
    );
    if (!res.ok) throw new Error(await res.text());
    const alignData = await res.json();
    const tsv = alignData.tsv || '';

    setProgress(100, 'Done!');

    // ── Result ────────────────────────────────────────────────────
    const blob    = new Blob([tsv], { type: 'text/tab-separated-values' });
    const url     = URL.createObjectURL(blob);
    const link    = document.getElementById('align-download-link');
    link.href     = url;
    link.download = 'alignment.aligns';
    document.getElementById('align-result-filename').textContent = 'alignment.aligns';

    const rowCount = tsv.trim().split('\n').length - 1;
    document.getElementById('align-stats').textContent =
      `${rowCount.toLocaleString()} alignment row${rowCount !== 1 ? 's' : ''}`;

    const alignLog = document.getElementById('align-log');
    if (alignData.log) {
      alignLog.textContent = alignData.log;
      alignLog.style.display = '';
    }

    initAlignTable(tsv);
    document.getElementById('align-table-section').style.display = '';

    setStatus('align', 'done');
  } catch (err) {
    document.getElementById('align-error-message').textContent = err.message;
    setStatus('align', 'error');
  } finally {
    alignBtn.disabled = false;
    if (_pRaf) { cancelAnimationFrame(_pRaf); _pRaf = null; }
  }
});
