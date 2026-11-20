type DiffWorkerReply = { id: number; html?: string; error?: string };
type PendingDiff = {
  resolve: (html: string) => void;
  timeout: ReturnType<typeof setTimeout>;
};

const DIFF_WORKER_TIMEOUT_MS = 8_000;
const MAX_DIFF_SOURCE_CHARACTERS = 64 * 1024;
let diffWorker: Worker | null = null;
let nextDiffRequestId = 1;
const pendingDiffs = new Map<number, PendingDiff>();

function diffErrorDocument(message: string): string {
  const escaped = message.replace(/[&<>"']/g, (character) => ({
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#039;',
  })[character] as string);
  return `<!doctype html><html><body><div style="padding:16px;border:1px solid #fecaca;border-radius:8px;background:#fef2f2;color:#b91c1c;font-family:system-ui,sans-serif">Diff unavailable: ${escaped}</div></body></html>`;
}

function resetDiffWorker(message: string): void {
  const errorHtml = diffErrorDocument(message);
  for (const { resolve, timeout } of pendingDiffs.values()) {
    clearTimeout(timeout);
    resolve(errorHtml);
  }
  pendingDiffs.clear();
  diffWorker?.terminate();
  diffWorker = null;
}

function getDiffWorker(): Worker {
  if (diffWorker) return diffWorker;

  const worker = new Worker(new URL('./diff.worker.ts', import.meta.url), { type: 'module' });
  worker.onmessage = (event: MessageEvent<DiffWorkerReply>) => {
    const pending = pendingDiffs.get(event.data.id);
    if (!pending) return;

    clearTimeout(pending.timeout);
    pendingDiffs.delete(event.data.id);
    pending.resolve(
      event.data.html ?? diffErrorDocument(event.data.error ?? 'Diff generation failed'),
    );
  };
  worker.onerror = () => resetDiffWorker('Diff worker failed');
  diffWorker = worker;
  return worker;
}

export async function generateDiffHtml(
  sourceCode: string,
  targetCode: string,
  theme: 'light' | 'dark' = 'light',
  changesOnly = false,
): Promise<string> {
  if (
    sourceCode.length > MAX_DIFF_SOURCE_CHARACTERS ||
    targetCode.length > MAX_DIFF_SOURCE_CHARACTERS
  ) {
    return diffErrorDocument('A source file exceeds the 64 KiB comparison limit');
  }

  try {
    const worker = getDiffWorker();
    const id = nextDiffRequestId;
    nextDiffRequestId += 1;

    return await new Promise<string>((resolve) => {
      const timeout = setTimeout(() => {
        if (pendingDiffs.has(id)) {
          resetDiffWorker('Diff generation exceeded the 8-second limit');
        }
      }, DIFF_WORKER_TIMEOUT_MS);

      pendingDiffs.set(id, { resolve, timeout });
      worker.postMessage({ id, source: sourceCode, target: targetCode, theme, changesOnly });
    });
  } catch (error) {
    return diffErrorDocument(error instanceof Error ? error.message : 'Diff generation failed');
  }
}
