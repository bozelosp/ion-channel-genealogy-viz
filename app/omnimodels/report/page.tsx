import type { Metadata } from "next";
import { OmnimodelReportShell } from "@/components/OmnimodelReportShell";

// Static shell that serves every /omnimodels/<channel>/<model> and
// /omnimodels/static/<channel>/<slug> URL in the exported deployment. The CDN
// router only rewrites URLs whose model key exists in the published key list,
// so unknown models receive a true 404 before ever reaching this shell.

export const dynamic = "force-static";

export const metadata: Metadata = {
  title: "Omnimodel report",
  robots: { index: false, follow: true },
};

export default function OmnimodelReportShellPage() {
  return <OmnimodelReportShell />;
}
