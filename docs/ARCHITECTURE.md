# Architecture

The application is a Next.js app that **deploys as a fully static site**:
every byte served in production is a pre-built file behind a CDN, with a
small amount of edge routing logic. `next dev` runs the complete app for
development; `npm run build:static` produces the verified production export.

## The Omnimodel report shell

The Omnimodel library contains 3,524 models. Prerendering one HTML page per
model was deliberately rejected (build cost, deploy size, a generated-code
monolith), and runtime rendering requires a server. Instead, every model URL
(`/omnimodels/<channel>/<model>` and `/omnimodels/static/<channel>/<slug>`)
is served by **one exported client shell** (`/omnimodels/report`):

- An edge function validates the requested model against a published key set
  derived from `public/data/omnimodel-shell-index.json` (both slug forms,
  case-insensitive). Valid URLs are rewritten to the shell object — the
  browser keeps the model URL; unknown models receive a **true HTTP 404**,
  never a soft client-side apology.
- The shell resolves `window.location.pathname` against the compact shell
  index (one cached fetch, ~488 KiB for all models), fetches that model's
  report object from `/data/omnimodel-reports/…`, **verifies the report's
  identity fields match the requested model**, then renders the shared
  report component with its plot images.
- Resolution is **fail-closed**: any fetch failure or identity mismatch
  renders an explicit "no substitute data is being shown" error. Nothing
  synthetic is ever rendered — the same rule the visualizer applies to its
  network dataset.

The listing pages are cheap and real: the index and the five channel pages
are genuinely prerendered, and every model grid paginates 60 cards at a time.

## The source-code route

`/api/source-code/<uniqueId>` accepts exactly the identifier grammar
`^\d{1,9}_[0-9A-Za-z_.+-]{1,48}-ID-\d{1,6}$` (bounded length, decoded before
validation). In development a bounded server proxy fetches the fixed
upstream origin (no redirects, request timeout, response size cap); in
production an edge function performs the same validation and rewrites onto
immutable stored objects. Anything else — other prefixes, traversal,
other extensions — is a true 404. Responses are plain text with
`X-Content-Type-Options: nosniff`.

## Static-export constraints

`output: 'export'` cannot contain middleware, non-static route handlers, or
runtime filesystem reads. `scripts/build-static.mjs` therefore moves the
dev-only surfaces aside for the duration of the export — the dev source-code
proxy route and the two runtime-rendered report routes replaced by the
shell — and restores them before finishing. Only the build that owns the
lock may restore paths or release it. Interrupts stop the child before
restoration; a surviving orphan child prevents conflicting recovery. The
move manifest is atomically replaced before each rename. Recovery refuses
invalid paths, unknown files, and source conflicts without discarding them.
Existing route edits are preserved; a dirty tree is not a build error.

`node --test scripts/build-static.test.mjs` exercises concurrency, failure,
interruption, crash recovery, and source restoration with tiny temporary
fixtures, without disrupting the checkout or running Next.js. The export
checks required files and compares report and shell-index counts against
the canonical index, rather than hard-coding a library size; plot counts
are reported. The default visualizer payload contains all models meeting
its initial two-copy threshold; lowering that threshold loads the
complete graph. Grouped layouts fetch only the table for their actual group
count. Canonical scientific downloads, including the aggregate layout and
dated network snapshot, remain available without adding to initial page load.

`node --test scripts/generate-visualizer-assets.test.mjs` checks the compact
graph's metadata, edges, ancestor candidates, and exact layout coordinates
against the complete source assets.

## Caching classes

Assets are grouped by volatility: HTML and route payloads are short-lived
(re-validated within minutes of a deploy), content-hashed build assets are
immutable for a year, and the scientific data (`/data/**`) is long-cached
with stale-while-revalidate — it changes only when the library is
regenerated. The initial visualizer graph is 1.70 MB uncompressed (~130 KB
gzip); the complete 9.98 MB graph is fetched only for one-copy models.
The two-group layout is 74 KB, instead of fetching the 12.69 MB aggregate.

## Design decisions

| Decision | Why |
|---|---|
| Fully static production | A read-only scientific explorer needs no server; removing the runtime removes its entire operational and attack surface. |
| One shell + edge key validation | Preserves every public model URL, keeps deploys small, and still yields real 404s for junk URLs because the edge holds the key list. |
| Fail-closed scientific rendering | A wrong-but-plausible figure is worse than a visible error; every data path validates identity and refuses to substitute. |
| True 404s everywhere | Missing content returns the 404 page with a 404 status — content never falls through to an application shell. |
| Grammar-gated source route | The route can only ever address the fixed source collection; the URL space is closed by construction. |
