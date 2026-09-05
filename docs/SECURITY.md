# Security model

This is a public, read-only scientific explorer with no accounts, sessions,
or tenant data. The assets that matter are **scientific integrity** (what
renders must be the canonical data or an explicit error), **availability**,
and **the supply chain** (what ships must be exactly the reviewed source).

## Controls

**Delivery.**
- Production is static files behind a CDN; there is no server-side request
  handling. Delivery is GET/HEAD only.
- Every miss is a true 404; object misses use the exported not-found page.
  Unknown model URLs are refused at the edge before any content is touched; invalid
  source identifiers are refused by grammar. Junk never resolves to an
  application shell.
- Browser protections are configured in `next.config.ts` for local serving
  and separately at the CDN for static delivery:
  a restrictive `Content-Security-Policy` (`default-src 'self'`, no
  objects, no framing, no external connections),
  `Strict-Transport-Security`, `X-Content-Type-Options: nosniff`,
  `X-Frame-Options: DENY`, `Referrer-Policy: strict-origin-when-cross-origin`,
  and a `Permissions-Policy` denying the listed device capabilities.
  Development additionally permits the framework's evaluation-based tooling;
  production CSP does not permit `unsafe-eval`. Inline bootstrap scripts
  and styles still require `unsafe-inline`; CSP is defense in depth.

**Application.**
- Scientific rendering is fail-closed: the visualizer validates its dataset
  and renders a blocking error rather than fabricated nodes; the report
  shell verifies the fetched report's identity fields against the requested
  model and shows "no substitute data is being shown" on any failure.
- Source-code diffs run in a Web Worker with bounded input (64 KiB / 2,000
  lines per side), bounded output (4 MiB), and full HTML escaping; the
  result renders inside a `sandbox=""` (unique-origin, script-less) iframe.
- The source-code route accepts only a strict identifier grammar; in
  development its proxy fetches a fixed origin with no redirects, a request
  timeout, and a response size cap.
- No `dangerouslySetInnerHTML`, no `eval`, no server-trusting client state.
- Both local server commands bind to loopback. Images are pre-generated;
  the unused server-side image optimizer is disabled in every build mode.

**Build and supply chain.**
- Installs use `npm ci --ignore-scripts` against the committed lockfile;
  the Node version is pinned.
- The static export checks required files and report/index counts against
  the canonical data. Its restoration and concurrency behavior is described
  in [Architecture](ARCHITECTURE.md#static-export-constraints).
- Data-generation scripts that read binary research artifacts use a
  restricted, NumPy-only unpickler with conversion budgets, and a shared
  path-safety module (symlink-refusing, root-contained writes). Operators
  must establish input-artifact provenance before running them. Conversion
  budgets apply after unpickling; this is not a sandbox for hostile pickles.

## Review references

Retrieved 2026-09-05: [OWASP ASVS 5.0.0](https://owasp.org/www-project-application-security-verification-standard/)
(application verification; released 2025-05-30),
[NIST SSDF 1.1 / SP 800-218](https://csrc.nist.gov/pubs/sp/800/218/final)
(final, 2022-02-03), [Next.js security release](https://nextjs.org/blog/august-2026-security-release)
(16.3.3, 2026-08-25), and [Node.js security release](https://nodejs.org/en/blog/vulnerability/july-2026-security-releases)
(22.23.2, 2026-07-29). Review official advisories as well as registry audit
results: a clean dependency audit alone does not establish security.

## Reporting

If you believe you have found a security issue, please open a GitHub issue
with minimal detail and a way to reach you, or contact the maintainer
directly; do not include exploit specifics in public.
