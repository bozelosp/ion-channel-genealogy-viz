import Link from "next/link";

// Exported as 404.html; the CDN serves it with a true 404 status for any
// object that does not exist, so missing content never falls through to an
// application shell.

export default function NotFound() {
  return (
    <div className="flex min-h-screen items-center justify-center bg-slate-50 px-4 text-slate-900 dark:bg-slate-950 dark:text-slate-100">
      <div className="max-w-md space-y-4 text-center">
        <p className="text-sm font-semibold uppercase tracking-wide text-slate-500">404</p>
        <h1 className="text-3xl font-bold">Page not found</h1>
        <p className="text-sm text-slate-600 dark:text-slate-400">
          The requested page does not exist on the Ion Channel Genealogy explorer.
        </p>
        <div className="flex justify-center gap-4 text-sm">
          <Link href="/" className="text-blue-600 hover:text-blue-700 dark:text-blue-300">
            Home
          </Link>
          <Link href="/visualizer" className="text-blue-600 hover:text-blue-700 dark:text-blue-300">
            Visualizer
          </Link>
          <Link href="/omnimodels" className="text-blue-600 hover:text-blue-700 dark:text-blue-300">
            Omnimodels
          </Link>
        </div>
      </div>
    </div>
  );
}
