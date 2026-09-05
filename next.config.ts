import type { NextConfig } from "next";

// STATIC_EXPORT=1 builds the fully static production deployment. Browser
// security headers are then attached at the CDN, and the dev-only server
// features are excluded by scripts/build-static.mjs.
const staticExport = process.env.STATIC_EXPORT === "1";

const developmentEvalSource = process.env.NODE_ENV === "development" ? " 'unsafe-eval'" : "";
const contentSecurityPolicy = [
  "default-src 'self'",
  "base-uri 'none'",
  "object-src 'none'",
  "frame-ancestors 'none'",
  "form-action 'self'",
  `script-src 'self' 'unsafe-inline'${developmentEvalSource}`,
  "style-src 'self' 'unsafe-inline'",
  "img-src 'self' data: blob:",
  "font-src 'self' data:",
  "connect-src 'self'",
  "worker-src 'self' blob:",
  "frame-src 'self'",
  "media-src 'none'",
  "manifest-src 'self'",
  "upgrade-insecure-requests",
].join("; ");

const browserSecurityHeaders = [
  { key: "Content-Security-Policy", value: contentSecurityPolicy },
  { key: "Referrer-Policy", value: "strict-origin-when-cross-origin" },
  { key: "Permissions-Policy", value: "accelerometer=(), camera=(), geolocation=(), gyroscope=(), magnetometer=(), microphone=(), payment=(), usb=()" },
  { key: "Strict-Transport-Security", value: "max-age=63072000" },
  { key: "X-Content-Type-Options", value: "nosniff" },
  { key: "X-Frame-Options", value: "DENY" },
  { key: "X-Permitted-Cross-Domain-Policies", value: "none" },
];

const nextConfig: NextConfig = {
  poweredByHeader: false,
  ...(staticExport ? { output: "export" as const } : {}),
  // All plots are pre-generated. The local server needs no image optimizer.
  images: { unoptimized: true },
  ...(staticExport ? {} : {
    async headers() {
      return [
        {
          source: "/:path*",
          headers: browserSecurityHeaders,
        },
        {
          source: "/data/:path*",
          headers: [
            { key: "Cache-Control", value: "public, max-age=3600, s-maxage=31536000, stale-while-revalidate=86400" },
          ],
        },
      ];
    },
  }),
};

export default nextConfig;
