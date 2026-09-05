import Link from "next/link";
import { notFound } from "next/navigation";
import { getChannelSummary, getOmnimodelManifest } from "@/lib/omnimodels";
import { ModelGrid } from "./model-grid";

type RouteParams = { channel: string };
type ChannelPageContext = { params: Promise<RouteParams> };

export const dynamic = "force-static";
export const dynamicParams = false;

export async function generateStaticParams(): Promise<RouteParams[]> {
  const channels = await getOmnimodelManifest();
  return channels.map(({ channel }) => ({ channel }));
}

export async function generateMetadata(ctx: ChannelPageContext) {
  const { channel } = await ctx.params;
  return {
    title: `${channel} Omnimodels`,
    description: `Markdown library for the ${channel} omnimodel fits.`,
  };
}

export default async function ChannelPage(ctx: ChannelPageContext) {
  const { channel } = await ctx.params;
  const channelSummary = await getChannelSummary(channel);

  if (!channelSummary) {
    notFound();
  }

  return (
    <div className="space-y-8">
      <div className="space-y-3">
        <Link
          href="/omnimodels"
          className="inline-flex items-center text-sm text-blue-600 hover:text-blue-700 dark:text-blue-300 dark:hover:text-blue-200"
        >
          ← All channels
        </Link>
        <div className="flex flex-wrap items-end justify-between gap-4">
          <div>
            <h1 className="text-3xl font-bold capitalize">{channelSummary.channel} omnimodels</h1>
            <p className="mt-2 max-w-2xl text-sm text-slate-600 dark:text-slate-400">
              {channelSummary.models.length} markdown exports generated from the omnimodel fitting pipeline.
            </p>
          </div>
          <div className="rounded-md border border-slate-200 bg-white px-4 py-2 text-sm dark:border-slate-800 dark:bg-slate-900">
            <span className="font-semibold text-slate-700 dark:text-slate-200">Source path</span>
            <div className="mt-1 font-mono text-xs text-slate-500 dark:text-slate-400">
              {channelSummary.channel}
            </div>
          </div>
        </div>
      </div>

      <ModelGrid
        channel={channelSummary.channel}
        models={channelSummary.models.map((model) => ({
          slug: model.slug,
          title: model.title,
        }))}
      />
    </div>
  );
}
