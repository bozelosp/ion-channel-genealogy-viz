import fs from "node:fs/promises";
import path from "node:path";
import { cache } from "react";
import omnimodelIndexJson from "@/public/data/omnimodel-index.json";

const PUBLIC_OMNIMODEL_ROOT = path.resolve(process.cwd(), "public", "data", "omnimodels");
const GENERATED_REPORT_ROOT = path.resolve(process.cwd(), "public", "data", "omnimodel-reports");

export interface GeneratedOmnimodelSection {
  key: string;
  title: string;
  value: unknown;
}

export interface GeneratedOmnimodelSummary {
  channel: string;
  originalSlug: string;
  routeSlug: string;
  title: string;
  summary: {
    suffix: string | null;
    gmaxName: string | null;
    states: string[];
    gates: Record<string, number | null> | null;
  };
}

export interface GeneratedOmnimodelReport extends GeneratedOmnimodelSummary {
  sections: GeneratedOmnimodelSection[];
}

export interface OmnimodelModelSummary {
  channel: string;
  slug: string;
  routeSlug: string;
  title: string;
  absoluteDir: string;
}

export interface OmnimodelChannelSummary {
  channel: string;
  models: OmnimodelModelSummary[];
}

const OMNIMODEL_INDEX = omnimodelIndexJson as GeneratedOmnimodelSummary[];
const byGeneratedRoute = new Map<string, GeneratedOmnimodelSummary>();
const byOriginalSlug = new Map<string, GeneratedOmnimodelSummary>();
const channelSummaries = new Map<string, OmnimodelChannelSummary>();

function routeKey(channel: string, slug: string): string {
  return `${channel.toLowerCase()}/${slug.toLowerCase()}`;
}

for (const item of OMNIMODEL_INDEX) {
  byGeneratedRoute.set(routeKey(item.channel, item.routeSlug), item);
  byOriginalSlug.set(routeKey(item.channel, item.originalSlug), item);

  const existing = channelSummaries.get(item.channel) ?? { channel: item.channel, models: [] };
  existing.models.push({
    channel: item.channel,
    slug: item.originalSlug,
    routeSlug: item.routeSlug,
    title: item.title,
    absoluteDir: path.join(PUBLIC_OMNIMODEL_ROOT, item.channel, item.originalSlug),
  });
  channelSummaries.set(item.channel, existing);
}

for (const channel of channelSummaries.values()) {
  channel.models.sort((a, b) => a.slug.localeCompare(b.slug));
}

const manifest = Array.from(channelSummaries.values()).sort((a, b) =>
  a.channel.localeCompare(b.channel),
);

export function listGeneratedOmnimodels(): GeneratedOmnimodelSummary[] {
  return OMNIMODEL_INDEX;
}

export function findGeneratedOmnimodelSummary(
  channel: string,
  routeSlug: string,
): GeneratedOmnimodelSummary | undefined {
  return byGeneratedRoute.get(routeKey(channel, routeSlug));
}

export function findOmnimodelByOriginalSlug(
  channel: string,
  originalSlug: string,
): GeneratedOmnimodelSummary | undefined {
  return byOriginalSlug.get(routeKey(channel, originalSlug));
}

export const readGeneratedOmnimodelReport = cache(async (
  channel: string,
  routeSlug: string,
): Promise<GeneratedOmnimodelReport | undefined> => {
  const summary = findGeneratedOmnimodelSummary(channel, routeSlug);
  if (!summary) return undefined;

  const reportPath = path.join(
    GENERATED_REPORT_ROOT,
    summary.channel.toLowerCase(),
    `${summary.routeSlug}.json`,
  );

  try {
    const report = JSON.parse(await fs.readFile(reportPath, "utf8")) as GeneratedOmnimodelReport;
    if (routeKey(report.channel, report.routeSlug) !== routeKey(summary.channel, summary.routeSlug)) {
      throw new Error(`Omnimodel report identity mismatch for ${summary.channel}/${summary.routeSlug}`);
    }
    return report;
  } catch (error) {
    if ((error as NodeJS.ErrnoException).code === "ENOENT") return undefined;
    throw error;
  }
});

export async function getOmnimodelManifest(): Promise<OmnimodelChannelSummary[]> {
  return manifest;
}

export async function getChannelSummary(
  channel: string,
): Promise<OmnimodelChannelSummary | undefined> {
  const normalized = channel.toLowerCase();
  return manifest.find((item) => item.channel.toLowerCase() === normalized);
}

export async function getModelSummary(
  channel: string,
  modelSlug: string,
): Promise<OmnimodelModelSummary | undefined> {
  const item = findOmnimodelByOriginalSlug(channel, modelSlug);
  if (!item) return undefined;
  return {
    channel: item.channel,
    slug: item.originalSlug,
    routeSlug: item.routeSlug,
    title: item.title,
    absoluteDir: path.join(PUBLIC_OMNIMODEL_ROOT, item.channel, item.originalSlug),
  };
}

export async function readOmnimodelReportByModel(
  channel: string,
  modelSlug: string,
): Promise<GeneratedOmnimodelReport | undefined> {
  const item = findOmnimodelByOriginalSlug(channel, modelSlug);
  return item ? readGeneratedOmnimodelReport(item.channel, item.routeSlug) : undefined;
}

const IMAGE_EXTENSION = /\.(?:png|jpe?g|gif|svg|webp)$/i;

export const listOmnimodelImageUrls = cache(async (
  channel: string,
  modelSlug: string,
): Promise<string[]> => {
  const item = findOmnimodelByOriginalSlug(channel, modelSlug);
  if (!item) return [];

  const imagesDir = path.join(PUBLIC_OMNIMODEL_ROOT, item.channel, item.originalSlug, "images");
  try {
    const entries = await fs.readdir(imagesDir, { withFileTypes: true });
    return entries
      .filter((entry) => entry.isFile() && IMAGE_EXTENSION.test(entry.name))
      .map((entry) => entry.name)
      .sort((a, b) => a.localeCompare(b))
      .map((filename) => buildOmnimodelAssetUrl(item.channel, item.originalSlug, `images/${filename}`));
  } catch (error) {
    if ((error as NodeJS.ErrnoException).code === "ENOENT") return [];
    throw error;
  }
});

function encodeSegment(segment: string): string {
  return encodeURIComponent(segment);
}

function normalizeAssetPath(assetPath: string): string[] {
  const segments = assetPath.replace(/\\/g, "/").split("/").filter(Boolean);
  if (segments.length === 0) return [];
  const imagesIndex = segments.map((segment) => segment.toLowerCase()).lastIndexOf("images");
  return imagesIndex >= 0 ? segments.slice(imagesIndex) : ["images", ...segments];
}

export function buildOmnimodelAssetUrl(
  channel: string,
  modelSlug: string,
  assetPath: string,
): string {
  const parts = normalizeAssetPath(assetPath);
  if (parts.length === 0) return assetPath;
  return `/data/omnimodels/${encodeSegment(channel)}/${encodeSegment(modelSlug)}/${parts.map(encodeSegment).join("/")}`;
}
