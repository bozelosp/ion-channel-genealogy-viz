import Link from "next/link";

export default function Home() {
  return (
    <div className="flex min-h-screen flex-col bg-white dark:bg-slate-950">
      <main className="flex flex-1 flex-col items-center justify-center px-4">
        <h1 className="text-4xl font-semibold tracking-tight text-slate-900 dark:text-white sm:text-5xl">
          Ion Channel Genealogy
        </h1>
        <p className="mt-4 max-w-lg text-center text-base text-slate-500 dark:text-slate-400">
          Explore the source code similarity relationships between ion channel models across ModelDB.
        </p>
        <div className="mt-10 flex gap-4">
          <Link
            href="/visualizer"
            className="rounded-md bg-slate-900 px-5 py-2.5 text-sm font-medium text-white transition-colors hover:bg-slate-800 dark:bg-white dark:text-slate-900 dark:hover:bg-slate-100"
          >
            Visualizer
          </Link>
          <Link
            href="/omnimodels"
            className="rounded-md border border-slate-300 px-5 py-2.5 text-sm font-medium text-slate-700 transition-colors hover:bg-slate-50 dark:border-slate-700 dark:text-slate-300 dark:hover:bg-slate-900"
          >
            Omnimodels
          </Link>
        </div>
      </main>

      <footer className="py-6 text-center text-xs text-slate-400 dark:text-slate-600">
        <div className="flex items-center justify-center gap-4">
          <a
            href="https://modeldb.science"
            target="_blank"
            rel="noopener noreferrer"
            className="hover:text-slate-600 dark:hover:text-slate-400"
          >
            ModelDB
          </a>
          <span>·</span>
          <a
            href="https://icg.neurotheory.ox.ac.uk/"
            target="_blank"
            rel="noopener noreferrer"
            className="hover:text-slate-600 dark:hover:text-slate-400"
          >
            ICG
          </a>
        </div>
      </footer>
    </div>
  );
}
