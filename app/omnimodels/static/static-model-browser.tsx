"use client";

import Link from "next/link";
import { useDeferredValue, useEffect, useMemo, useState } from "react";

const PAGE_SIZE = 60;

interface StaticModelListItem {
  channel: string;
  routeSlug: string;
  title: string;
  suffix: string | null;
  states: string[];
}

export function StaticModelBrowser({ models }: { models: StaticModelListItem[] }) {
  const [query, setQuery] = useState("");
  const [channel, setChannel] = useState("all");
  const [visibleCount, setVisibleCount] = useState(PAGE_SIZE);
  const deferredQuery = useDeferredValue(query.trim().toLowerCase());

  const channels = useMemo(
    () => Array.from(new Set(models.map((model) => model.channel))).sort((a, b) => a.localeCompare(b)),
    [models],
  );

  const filtered = useMemo(() => models.filter((model) => {
    if (channel !== "all" && model.channel !== channel) return false;
    if (!deferredQuery) return true;
    const normalized = deferredQuery.endsWith(".mod") ? deferredQuery.slice(0, -4) : deferredQuery;
    return [model.title, model.routeSlug, model.suffix ?? "", ...model.states]
      .some((value) => value.toLowerCase().includes(normalized));
  }), [channel, deferredQuery, models]);

  useEffect(() => {
    setVisibleCount(PAGE_SIZE);
  }, [channel, deferredQuery]);

  const visible = filtered.slice(0, visibleCount);

  return (
    <div className="space-y-6">
      <div className="grid gap-4 rounded-xl border border-slate-200 bg-white p-4 sm:grid-cols-[minmax(0,1fr)_12rem] dark:border-slate-800 dark:bg-slate-900">
        <label>
          <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">Search reports</span>
          <input
            type="search"
            value={query}
            onChange={(event) => setQuery(event.target.value)}
            placeholder="Model, suffix, state, or slug"
            className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-3 py-2 text-sm text-slate-900 focus:border-blue-500 focus:outline-none focus:ring-2 focus:ring-blue-500/50 dark:border-slate-700 dark:bg-slate-950 dark:text-slate-100"
          />
        </label>
        <label>
          <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">Channel</span>
          <select
            value={channel}
            onChange={(event) => setChannel(event.target.value)}
            className="mt-1 w-full rounded-lg border border-slate-300 bg-white px-3 py-2 text-sm text-slate-900 focus:border-blue-500 focus:outline-none focus:ring-2 focus:ring-blue-500/50 dark:border-slate-700 dark:bg-slate-950 dark:text-slate-100"
          >
            <option value="all">All channels</option>
            {channels.map((value) => <option key={value} value={value}>{value}</option>)}
          </select>
        </label>
      </div>

      <p className="text-sm text-slate-600 dark:text-slate-400" aria-live="polite">
        Showing <strong>{visible.length.toLocaleString()}</strong> of <strong>{filtered.length.toLocaleString()}</strong> matching reports
        {filtered.length !== models.length ? ` (${models.length.toLocaleString()} total)` : ""}.
      </p>

      <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-3">
        {visible.map((item) => (
          <Link
            key={`${item.channel}/${item.routeSlug}`}
            href={`/omnimodels/static/${item.channel.toLowerCase()}/${item.routeSlug}`}
            prefetch={false}
            className="rounded-lg border border-slate-200 bg-white p-4 text-sm shadow-sm transition hover:border-blue-400 hover:shadow-md dark:border-slate-800 dark:bg-slate-900"
          >
            <div className="font-semibold text-slate-800 dark:text-slate-100">{item.title}</div>
            <div className="mt-1 text-xs text-slate-500 dark:text-slate-400">
              {item.channel} · Suffix: {item.suffix ?? "n/a"}
            </div>
            <div className="text-xs text-slate-500 dark:text-slate-400">
              States: {item.states.length ? item.states.join(", ") : "n/a"}
            </div>
          </Link>
        ))}
      </div>

      {visible.length < filtered.length && (
        <div className="flex justify-center">
          <button
            type="button"
            onClick={() => setVisibleCount((count) => Math.min(filtered.length, count + PAGE_SIZE))}
            className="rounded-lg border border-slate-300 bg-white px-5 py-2.5 text-sm font-medium text-slate-700 hover:bg-slate-50 dark:border-slate-700 dark:bg-slate-900 dark:text-slate-200 dark:hover:bg-slate-800"
          >
            Show {Math.min(PAGE_SIZE, filtered.length - visible.length).toLocaleString()} more
          </button>
        </div>
      )}

      {filtered.length === 0 && (
        <p className="rounded-lg border border-dashed border-slate-300 p-8 text-center text-sm text-slate-500 dark:border-slate-700">
          No generated reports match those filters.
        </p>
      )}
    </div>
  );
}
