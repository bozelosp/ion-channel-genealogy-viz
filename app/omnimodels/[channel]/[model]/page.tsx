import { notFound } from "next/navigation";
import Link from "next/link";
import { OmnimodelGeneratedReport } from "@/components/OmnimodelGeneratedReport";
import {
  findOmnimodelByOriginalSlug,
  listOmnimodelImageUrls,
  readOmnimodelReportByModel,
} from "@/lib/omnimodels";

type RouteParams = { channel: string; model: string };
type ModelPageContext = { params: Promise<RouteParams> };

export const dynamic = "force-static";
export const dynamicParams = true;
export const revalidate = 31_536_000;

export function generateStaticParams(): RouteParams[] {
  return [];
}

export async function generateMetadata(ctx: ModelPageContext) {
  const params = await ctx.params;
  const summary = findOmnimodelByOriginalSlug(params.channel, params.model);
  if (!summary) return { title: "Omnimodel not found" };
  return {
    title: `${summary.title} – ${summary.channel} omnimodel`,
    description: `Generated scientific report for the ${summary.title} model in the ${summary.channel} omnimodel library.`,
  };
}

export default async function ModelPage(ctx: ModelPageContext) {
  const params = await ctx.params;
  const report = await readOmnimodelReportByModel(params.channel, params.model);
  if (!report) notFound();

  const imageUrls = await listOmnimodelImageUrls(report.channel, report.originalSlug);
  return (
    <div className="space-y-8">
      <Link
        href={`/omnimodels/${encodeURIComponent(report.channel)}`}
        className="inline-flex items-center text-sm text-blue-600 hover:text-blue-700 dark:text-blue-300 dark:hover:text-blue-200"
      >
        ← {report.channel} omnimodels
      </Link>
      <OmnimodelGeneratedReport report={report} imageUrls={imageUrls} />
    </div>
  );
}
