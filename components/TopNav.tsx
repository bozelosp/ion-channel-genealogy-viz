"use client";

import Link from "next/link";
import type { ReactNode } from "react";

export default function TopNav({ right }: { right?: ReactNode }) {
  const linkClass =
    "rounded-lg border border-slate-300 px-3 py-1.5 text-slate-700 transition-colors hover:bg-slate-100 dark:border-slate-700 dark:text-slate-300 dark:hover:bg-slate-800";

  return (
    <nav className="sticky top-0 z-50 border-b border-slate-200 bg-white/70 backdrop-blur dark:border-slate-800 dark:bg-slate-900/70">
      <div className="mx-auto flex h-16 max-w-7xl items-center justify-between px-4">
        <Link href="/" className="text-lg font-semibold text-slate-900 dark:text-slate-100">
          Ion Channel Genealogy
        </Link>
        <div className="flex items-center gap-2 text-sm md:gap-3">
          <Link href="/visualizer" className={linkClass}>
            Visualizer
          </Link>
          <Link href="/omnimodels" className={linkClass}>
            Omnimodels
          </Link>
          {right}
        </div>
      </div>
    </nav>
  );
}
