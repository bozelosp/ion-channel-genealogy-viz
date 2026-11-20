const modelIdPattern = /^\d{1,9}_[0-9A-Za-z_.+-]{1,48}-ID-\d{1,6}$/;

export function isModelSourceId(id: string): boolean {
  return id.length <= 64 && modelIdPattern.test(id);
}

export function comparisonPath(left: string, right: string): string {
  return '/compare?' + new URLSearchParams({ left, right });
}

export async function fetchModelSource(id: string, signal?: AbortSignal): Promise<string> {
  const response = await fetch('/api/source-code/' + encodeURIComponent(id), { signal });
  if (!response.ok) throw new Error(`Source unavailable for ${id} (HTTP ${response.status}).`);
  return response.text();
}
