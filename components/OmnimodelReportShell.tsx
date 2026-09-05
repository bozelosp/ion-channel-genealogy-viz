"use client";

import Link from "next/link";
import { useEffect, useState } from "react";
import { OmnimodelGeneratedReport } from "@/components/OmnimodelGeneratedReport";
import type { GeneratedOmnimodelReport } from "@/lib/omnimodels";

// The static export serves every model URL from this one client shell: an
// edge function validates the requested model against the published key
// list and rewrites valid URLs here, so the shell only has to resolve the
// original pathname against the compact lookup index and fetch that model's
// report object. Resolution is fail-closed: any missing or mismatched data
// renders an explicit error, never substitute content.

const SHELL_INDEX_URL = "/data/omnimodel-shell-index.json";

type ShellIndexEntry = [
  channel: string,
  originalSlug: string,
  routeSlug: string,
  title: string,
  imageNames: string[],
];

interface ShellIndex {
  count: number;
  models: Record<string, ShellIndexEntry>;
  byRoute: Record<string, string>;
}

type ShellState =
  | { phase: "loading" }
  | { phase: "picker" }
  | { phase: "not-found"; requested: string }
  | { phase: "error"; message: string }
  | { phase: "ready"; report: GeneratedOmnimodelReport; imageUrls: string[]; channel: string };

let shellIndexPromise: Promise<ShellIndex> | null = null;

function loadShellIndex(): Promise<ShellIndex> {
  if (!shellIndexPromise) {
    shellIndexPromise = (async () => {
      const response = await fetch(SHELL_INDEX_URL);
      if (!response.ok) {
        throw new Error(`Shell index request failed (${response.status})`);
      }
      const index = (await response.json()) as ShellIndex;
      if (!index || typeof index.count !== "number" || !index.models || !index.byRoute) {
        throw new Error("Shell index has an unexpected shape");
      }
      return index;
    })();
    shellIndexPromise.catch(() => {
      shellIndexPromise = null;
    });
  }
  return shellIndexPromise;
}

function decodeSegment(segment: string): string | null {
  try {
    return decodeURIComponent(segment);
  } catch {
    return null;
  }
}

function resolveShellKey(pathname: string, index: ShellIndex): string | null {
  const segments = pathname.split("/").filter(Boolean).map(decodeSegment);
  if (segments.some((segment) => segment === null)) return null;
  const parts = segments as string[];

  if (parts[0]?.toLowerCase() !== "omnimodels") return null;

  // /omnimodels/static/<channel>/<routeSlug>
  if (parts.length === 4 && parts[1].toLowerCase() === "static") {
    const routeKey = `${parts[2].toLowerCase()}/${parts[3].toLowerCase()}`;
    return index.byRoute[routeKey] ?? null;
  }

  // /omnimodels/<channel>/<originalSlug>
  if (parts.length === 3 && !["static", "report"].includes(parts[1].toLowerCase())) {
    const key = `${parts[1].toLowerCase()}/${parts[2].toLowerCase()}`;
    return key in index.models ? key : null;
  }

  return null;
}

function buildImageUrls(entry: ShellIndexEntry): string[] {
  const [channel, originalSlug, , , imageNames] = entry;
  return imageNames.map(
    (name) =>
      `/data/omnimodels/${encodeURIComponent(channel)}/${encodeURIComponent(originalSlug)}/images/${encodeURIComponent(name)}`,
  );
}

export function OmnimodelReportShell() {
  const [state, setState] = useState<ShellState>({ phase: "loading" });

  useEffect(() => {
    const controller = new AbortController();

    async function resolve() {
      const pathname = window.location.pathname.replace(/\/+$/, "") || "/";
      if (pathname === "/omnimodels/report") {
        setState({ phase: "picker" });
        return;
      }

      try {
        const index = await loadShellIndex();
        const key = resolveShellKey(pathname, index);
        const entry = key ? index.models[key] : undefined;
        if (!entry) {
          document.title = "Omnimodel not found";
          setState({ phase: "not-found", requested: pathname });
          return;
        }

        const [channel, originalSlug, routeSlug, title] = entry;
        const reportUrl = `/data/omnimodel-reports/${encodeURIComponent(channel.toLowerCase())}/${encodeURIComponent(routeSlug)}.json`;
        const response = await fetch(reportUrl, { signal: controller.signal });
        if (!response.ok) {
          throw new Error(`Report request failed (${response.status})`);
        }
        const report = (await response.json()) as GeneratedOmnimodelReport;
        if (
          report.channel?.toLowerCase() !== channel.toLowerCase() ||
          report.routeSlug?.toLowerCase() !== routeSlug.toLowerCase() ||
          report.originalSlug?.toLowerCase() !== originalSlug.toLowerCase() ||
          !Array.isArray(report.sections)
        ) {
          throw new Error("Report identity does not match the requested model");
        }

        document.title = `${title} – ${report.channel} omnimodel`;
        setState({ phase: "ready", report, imageUrls: buildImageUrls(entry), channel: report.channel });
      } catch (error) {
        if (controller.signal.aborted) return;
        console.error("Omnimodel report resolution failed:", error);
        setState({
          phase: "error",
          message: "The requested scientific report could not be loaded. No substitute data is being shown.",
        });
      }
    }

    resolve();
    return () => controller.abort();
  }, []);

  if (state.phase === "loading") {
    return (
      <div className="space-y-4" aria-busy="true" aria-live="polite">
        <div className="h-4 w-40 animate-pulse rounded bg-slate-200 dark:bg-slate-800" />
        <div className="h-9 w-2/3 animate-pulse rounded bg-slate-200 dark:bg-slate-800" />
        <div className="h-64 animate-pulse rounded-md border border-slate-200 bg-slate-100 dark:border-slate-800 dark:bg-slate-900" />
      </div>
    );
  }

  if (state.phase === "picker") {
    return (
      <div className="space-y-4">
        <h1 className="text-3xl font-semibold">Omnimodel report viewer</h1>
        <p className="text-sm text-slate-600 dark:text-slate-400">
          This page renders a single model report. Pick a model from the library to view its report.
        </p>
        <Link href="/omnimodels" className="text-sm text-blue-600 hover:text-blue-700 dark:text-blue-300">
          Browse omnimodels →
        </Link>
      </div>
    );
  }

  if (state.phase === "not-found") {
    return (
      <div className="space-y-4">
        <h1 className="text-3xl font-semibold">Omnimodel not found</h1>
        <p className="text-sm text-slate-600 dark:text-slate-400">
          No generated report exists for <code className="font-mono">{state.requested}</code>.
        </p>
        <Link href="/omnimodels" className="text-sm text-blue-600 hover:text-blue-700 dark:text-blue-300">
          ← All channels
        </Link>
      </div>
    );
  }

  if (state.phase === "error") {
    return (
      <div className="space-y-4" role="alert">
        <h1 className="text-3xl font-semibold">Report unavailable</h1>
        <p className="text-sm text-red-700 dark:text-red-400">{state.message}</p>
        <Link href="/omnimodels" className="text-sm text-blue-600 hover:text-blue-700 dark:text-blue-300">
          ← All channels
        </Link>
      </div>
    );
  }

  return (
    <div className="space-y-8">
      <Link
        href={`/omnimodels/${encodeURIComponent(state.channel)}`}
        className="inline-flex items-center text-sm text-blue-600 hover:text-blue-700 dark:text-blue-300 dark:hover:text-blue-200"
      >
        ← {state.channel} omnimodels
      </Link>
      <OmnimodelGeneratedReport report={state.report} imageUrls={state.imageUrls} />
    </div>
  );
}
