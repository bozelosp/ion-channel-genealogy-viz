import { NextResponse } from "next/server";
import { getOmnimodelManifest } from "@/lib/omnimodels";

export const dynamic = "force-static";

export async function GET() {
  const manifest = await getOmnimodelManifest();
  const simplified = manifest.map((channel) => ({
    channel: channel.channel,
    models: channel.models.map((model) => ({
      slug: model.slug,
    })),
  }));

  return NextResponse.json(
    { channels: simplified },
    { headers: { "Cache-Control": "public, max-age=3600, s-maxage=31536000, stale-while-revalidate=86400" } },
  );
}
