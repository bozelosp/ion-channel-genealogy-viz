import { diffLines } from 'diff';

const MAX_SOURCE_CHARACTERS = 64 * 1024;
const MAX_SOURCE_LINES = 2_000;
const MAX_OUTPUT_CHARACTERS = 4 * 1024 * 1024;

type DiffRequest = {
  id: number;
  source: string;
  target: string;
  theme: 'light' | 'dark';
};

type DiffResponse = {
  id: number;
  html?: string;
  error?: string;
};

const workerScope = self as unknown as {
  onmessage: ((event: MessageEvent<DiffRequest>) => void) | null;
  postMessage: (response: DiffResponse) => void;
};

function escapeHtml(text: string): string {
  return text.replace(/[&<>"']/g, (character) => ({
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#039;',
  })[character] as string);
}

function validateSource(source: string, label: string): void {
  if (source.length > MAX_SOURCE_CHARACTERS) {
    throw new Error(`${label} exceeds the 64 KiB comparison limit`);
  }
  if (source.split('\n').length > MAX_SOURCE_LINES) {
    throw new Error(`${label} exceeds the 2,000-line comparison limit`);
  }
}

function generateHtmlDiff(source: string, target: string, theme: 'light' | 'dark'): string {
  validateSource(source, 'Source');
  validateSource(target, 'Target');

  const rows: string[] = [];
  let outputCharacters = 0;
  let sourceLine = 1;
  let targetLine = 1;

  const pushRow = (row: string) => {
    outputCharacters += row.length;
    if (outputCharacters > MAX_OUTPUT_CHARACTERS) {
      throw new Error('Generated diff exceeds the 4 MiB display limit');
    }
    rows.push(row);
  };

  for (const part of diffLines(source, target)) {
    const lines = part.value.split('\n');
    if (lines.at(-1) === '') lines.pop();

    for (const line of lines) {
      if (part.added) {
        pushRow(`<tr><td class="diff_header"></td><td></td><td class="diff_header">${targetLine}</td><td class="diff_add">+&nbsp;${escapeHtml(line)}</td></tr>`);
        targetLine += 1;
      } else if (part.removed) {
        pushRow(`<tr><td class="diff_header">${sourceLine}</td><td class="diff_sub">-&nbsp;${escapeHtml(line)}</td><td class="diff_header"></td><td></td></tr>`);
        sourceLine += 1;
      } else {
        const escaped = escapeHtml(line);
        pushRow(`<tr><td class="diff_header">${sourceLine}</td><td>&nbsp;${escaped}</td><td class="diff_header">${targetLine}</td><td>&nbsp;${escaped}</td></tr>`);
        sourceLine += 1;
        targetLine += 1;
      }
    }
  }

  return `<!doctype html>
<html data-theme="${theme}">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Code Diff</title>
  <style>
    :root{--bg:#fff;--text:#0f172a;--border:#e2e8f0;--header:#e5e7eb;--add:#aaffaa;--sub:#ffaaaa;--hi:#111827}
    :root[data-theme='dark']{--bg:#0f172a;--text:#e2e8f0;--border:#334155;--header:#1f2937;--add:#14532d;--sub:#7f1d1d;--hi:#f8fafc}
    html,body{background:var(--bg);color:var(--text);margin:0;padding:12px;min-height:100%;overflow:auto}
    table{font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace;border:1px solid var(--border);border-collapse:collapse;width:100%;font-size:13px}
    td,th{border:1px solid var(--border);padding:4px 8px;overflow-wrap:anywhere}
    td:nth-child(1),td:nth-child(3){width:1%;white-space:nowrap;text-align:right;background:var(--header);font-size:12px}
    td:nth-child(2),td:nth-child(4){width:49%;white-space:pre-wrap;word-break:break-all}
    .diff_header{background:var(--header)} .diff_add{background:var(--add);color:var(--hi)} .diff_sub{background:var(--sub);color:var(--hi)}
  </style>
</head>
<body>
  <table class="diff">
    <thead><tr><th colspan="2">Source Code 1</th><th colspan="2">Source Code 2</th></tr></thead>
    <tbody>${rows.join('')}</tbody>
  </table>
</body>
</html>`;
}

workerScope.onmessage = (event: MessageEvent<DiffRequest>) => {
  const { id, source, target, theme } = event.data;
  const response: DiffResponse = { id };

  try {
    if (typeof source !== 'string' || typeof target !== 'string') {
      throw new Error('Comparison inputs must be strings');
    }
    response.html = generateHtmlDiff(source, target, theme === 'dark' ? 'dark' : 'light');
  } catch (error) {
    response.error = error instanceof Error ? error.message : 'Diff generation failed';
  }

  workerScope.postMessage(response);
};

export {};
