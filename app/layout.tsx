import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "Ion Channel Genealogy Visualizer",
  description: "Interactive visualization and analysis of ion channel models from ModelDB and ICG databases",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body className="antialiased">{children}</body>
    </html>
  );
}
