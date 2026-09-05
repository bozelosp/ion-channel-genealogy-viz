import Link from "next/link";
import { getOmnimodelManifest } from "@/lib/omnimodels";

export const dynamic = "force-static";

export default async function OmnimodelsIndexPage() {
  const manifest = await getOmnimodelManifest();

  return (
    <div className="space-y-10">
      <header>
        <h1 className="text-4xl font-bold tracking-tight">Omnimodels</h1>
      </header>

      <section className="grid gap-6 sm:grid-cols-2">
        {manifest.map((channel) => (
          <Link
            key={channel.channel}
            href={`/omnimodels/${encodeURIComponent(channel.channel)}`}
            className="group rounded-xl border border-slate-200 bg-white p-6 shadow-sm transition hover:border-blue-400 hover:shadow-lg dark:border-slate-800 dark:bg-slate-900"
          >
            <div className="flex items-center justify-between">
              <h2 className="text-2xl font-semibold capitalize">
                {channel.channel}
              </h2>
              <span className="rounded-full bg-blue-50 px-3 py-1 text-sm font-medium text-blue-600 group-hover:bg-blue-500 group-hover:text-white dark:bg-blue-900/40 dark:text-blue-300">
                {channel.models.length} models
              </span>
            </div>
          </Link>
        ))}
      </section>
    </div>
  );
}
