import { NextRequest, NextResponse } from 'next/server';

export const runtime = 'nodejs';

const SOURCE_CODE_BASE = 'https://ion-channels.s3.eu-west-1.amazonaws.com/static/modelDB';
const MAX_UNIQUE_ID_LENGTH = 64;
const MAX_SOURCE_BYTES = 256 * 1024;
const UPSTREAM_TIMEOUT_MS = 8_000;
const UNIQUE_ID_PATTERN = /^\d{1,9}_[0-9A-Za-z_.+-]{1,48}-ID-\d{1,6}$/;

class UpstreamBodyTooLargeError extends Error {}

async function readTextWithinLimit(response: Response): Promise<string> {
  const declaredLength = Number(response.headers.get('content-length'));
  if (Number.isFinite(declaredLength) && declaredLength > MAX_SOURCE_BYTES) {
    await response.body?.cancel();
    throw new UpstreamBodyTooLargeError('Source-code response exceeded the byte limit');
  }

  if (!response.body) return '';

  const reader = response.body.getReader();
  const chunks: Uint8Array[] = [];
  let totalBytes = 0;

  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    totalBytes += value.byteLength;
    if (totalBytes > MAX_SOURCE_BYTES) {
      await reader.cancel();
      throw new UpstreamBodyTooLargeError('Source-code response exceeded the byte limit');
    }
    chunks.push(value);
  }

  const bytes = new Uint8Array(totalBytes);
  let offset = 0;
  for (const chunk of chunks) {
    bytes.set(chunk, offset);
    offset += chunk.byteLength;
  }
  return new TextDecoder('utf-8', { fatal: false }).decode(bytes);
}

export async function GET(_request: NextRequest, context: { params: Promise<{ uniqueId: string }> }) {
  const { uniqueId } = await context.params;

  if (!uniqueId) {
    return NextResponse.json({ error: 'uniqueId is required' }, { status: 400 });
  }
  if (uniqueId.length > MAX_UNIQUE_ID_LENGTH || !UNIQUE_ID_PATTERN.test(uniqueId)) {
    return NextResponse.json({ error: 'Invalid uniqueId' }, { status: 400 });
  }

  const url = `${SOURCE_CODE_BASE}/source_code/${encodeURIComponent(uniqueId)}.mod`;

  try {
    const upstream = await fetch(url, {
      redirect: 'error',
      signal: AbortSignal.timeout(UPSTREAM_TIMEOUT_MS),
      next: { revalidate: 31_536_000 },
    });

    if (!upstream.ok) {
      await upstream.body?.cancel();
      if (upstream.status === 404) {
        return NextResponse.json(
          { error: 'Source code not found' },
          { status: 404, headers: { 'Cache-Control': 'public, s-maxage=300' } },
        );
      }
      return NextResponse.json(
        { error: 'Source-code upstream is unavailable' },
        { status: 502, headers: { 'Cache-Control': 'no-store' } },
      );
    }

    const text = await readTextWithinLimit(upstream);
    return new NextResponse(text, {
      status: 200,
      headers: {
        'Content-Type': 'text/plain; charset=utf-8',
        'Cache-Control': 'public, max-age=3600, s-maxage=31536000, stale-while-revalidate=86400',
        'X-Content-Type-Options': 'nosniff',
      },
    });
  } catch (error) {
    if (error instanceof UpstreamBodyTooLargeError) {
      console.warn('source-code upstream response rejected: size limit exceeded');
      return NextResponse.json(
        { error: 'Source-code upstream response is too large' },
        { status: 502, headers: { 'Cache-Control': 'no-store' } },
      );
    }
    if (error instanceof Error && (error.name === 'TimeoutError' || error.name === 'AbortError')) {
      console.warn('source-code upstream request timed out');
      return NextResponse.json(
        { error: 'Source-code upstream timed out' },
        { status: 504, headers: { 'Cache-Control': 'no-store' } },
      );
    }
    console.error('source-code upstream request failed');
    return NextResponse.json(
      { error: 'Unable to load source code' },
      { status: 502, headers: { 'Cache-Control': 'no-store' } },
    );
  }
}
