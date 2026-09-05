import { listGeneratedOmnimodels } from "@/lib/omnimodels";
import { StaticModelBrowser } from "./static-model-browser";

export const dynamic = "force-static";

export default function OmnimodelStaticIndexPage() {
  const reports = listGeneratedOmnimodels().map((report) => ({
    channel: report.channel,
    routeSlug: report.routeSlug,
    title: report.title,
    suffix: report.summary.suffix,
    states: report.summary.states,
  }));

  return (
    <div className="space-y-8">
      <header className="space-y-2">
        <h1 className="text-4xl font-bold tracking-tight">Omnimodels (Static)</h1>
        <p className="text-sm text-slate-600 dark:text-slate-400">
          Browse all {reports.length.toLocaleString()} generated scientific reports without loading every card at once.
        </p>
      </header>
      <StaticModelBrowser models={reports} />
    </div>
  );
}
