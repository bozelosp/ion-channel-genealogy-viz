import { notFound } from "next/navigation";
import type { Metadata } from "next";
import { OmnimodelGeneratedReport } from "@/components/OmnimodelGeneratedReport";
import {
  findGeneratedOmnimodelSummary,
  listOmnimodelImageUrls,
  readGeneratedOmnimodelReport,
} from "@/lib/omnimodels";

type RouteParams = { channel: string; slug: string };
type PageContext = { params: Promise<RouteParams> };

export const dynamic = "force-static";
export const dynamicParams = true;
export const revalidate = 31_536_000;

// Reports are generated and cached on first request instead of expanding every report during builds.
export function generateStaticParams(): RouteParams[] {
  return [];
}

export async function generateMetadata(ctx: PageContext): Promise<Metadata> {
  const params = await ctx.params;
  const report = findGeneratedOmnimodelSummary(params.channel, params.slug);
  if (!report) return { title: "Omnimodel not found" };
  return {
    title: `${report.title} – ${report.channel} omnimodel`,
    description: `Generated omnimodel report for ${report.title} in the ${report.channel} channel class.`,
  };
}

export default async function OmnimodelStaticPage(ctx: PageContext) {
  const params = await ctx.params;
  const report = await readGeneratedOmnimodelReport(params.channel, params.slug);
  if (!report) notFound();

  const imageUrls = await listOmnimodelImageUrls(report.channel, report.originalSlug);
  return <OmnimodelGeneratedReport report={report} imageUrls={imageUrls} />;
}
