import type { ReactNode } from "react";
import TopNav from "../../components/TopNav";

export default function OmnimodelsLayout({
  children,
}: {
  children: ReactNode;
}) {
  return (
    <div className="min-h-screen bg-slate-50 text-slate-900 dark:bg-slate-950 dark:text-slate-100">
      <TopNav />
      <main className="mx-auto w-full max-w-6xl px-4 py-10 sm:py-12">{children}</main>
    </div>
  );
}
