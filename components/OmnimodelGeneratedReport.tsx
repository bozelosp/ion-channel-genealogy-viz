import { Fragment, type ReactElement } from "react";
import type { GeneratedOmnimodelReport } from "@/lib/omnimodels";

function isPrimitive(value: unknown): value is string | number | boolean | null {
  return (
    value === null ||
    typeof value === "string" ||
    typeof value === "number" ||
    typeof value === "boolean"
  );
}

function renderArray(value: unknown[]): ReactElement {
  if (value.length === 0) {
    return <p>[]</p>;
  }

  const simple = value.every((item) => isPrimitive(item));

  if (simple && value.length <= 20) {
    return (
      <ul>
        {value.map((item, index) => (
          <li key={index}>{String(item)}</li>
        ))}
      </ul>
    );
  }

  return <pre>{JSON.stringify(value, null, 2)}</pre>;
}

function renderObject(value: Record<string, unknown>): ReactElement {
  const entries = Object.entries(value);

  if (entries.length === 0) {
    return <p>{"{}"}</p>;
  }

  const simple = entries.every(([, item]) => isPrimitive(item));

  if (simple && entries.length <= 20) {
    return (
      <dl>
        {entries.map(([key, item]) => (
          <Fragment key={key}>
            <dt>{key}</dt>
            <dd>{item === null ? "null" : String(item)}</dd>
          </Fragment>
        ))}
      </dl>
    );
  }

  return <pre>{JSON.stringify(value, null, 2)}</pre>;
}

function renderValue(value: unknown): ReactElement {
  if (value === null || value === undefined) {
    return <p><em>None</em></p>;
  }

  if (typeof value === "string" || typeof value === "number" || typeof value === "boolean") {
    return <p>{String(value)}</p>;
  }

  if (Array.isArray(value)) {
    return renderArray(value);
  }

  if (typeof value === "object") {
    return renderObject(value as Record<string, unknown>);
  }

  return <pre>{JSON.stringify(value, null, 2)}</pre>;
}

interface OmnimodelGeneratedReportProps {
  report: GeneratedOmnimodelReport;
  imageUrls?: string[];
}

export function OmnimodelGeneratedReport({ report, imageUrls = [] }: OmnimodelGeneratedReportProps) {
  const { summary } = report;

  return (
    <div className="space-y-8">
      <header className="space-y-2">
        <p className="text-sm text-slate-500">Generated omnimodel report</p>
        <h1 className="text-3xl font-semibold">{report.title}</h1>
        <p className="text-slate-600">
          Channel class: <strong>{report.channel}</strong>
        </p>
        <div className="text-sm text-slate-600">
          <div>Suffix: {summary.suffix ?? <em>unknown</em>}</div>
          <div>GMAX name: {summary.gmaxName ?? <em>unknown</em>}</div>
          <div>
            States: {summary.states.length > 0 ? summary.states.join(", ") : <em>not specified</em>}
          </div>
          <div>
            Gates: {summary.gates ? (
              <span>
                {Object.entries(summary.gates)
                  .map(([key, value]) => `${key}: ${value ?? "null"}`)
                  .join(", ")}
              </span>
            ) : (
              <em>not specified</em>
            )}
          </div>
        </div>
        <div className="text-xs text-slate-400">
          Route slug: <code>{report.routeSlug}</code> · Original slug: <code>{report.originalSlug}</code>
        </div>
      </header>

      {report.sections.map((section) => (
        <section key={section.key} className="space-y-2">
          <h2 className="text-2xl font-semibold">{section.title}</h2>
          <div className="rounded-md border border-slate-200 bg-white p-4 text-sm text-slate-700 dark:border-slate-700 dark:bg-slate-900 dark:text-slate-200">
            {renderValue(section.value)}
          </div>
        </section>
      ))}

      {imageUrls.length > 0 && (
        <section className="space-y-4">
          <h2 className="text-2xl font-semibold">Figures</h2>
          <div className="grid gap-6 lg:grid-cols-2">
            {imageUrls.map((src) => (
              <figure key={src} className="rounded-md border border-slate-200 bg-white p-4 dark:border-slate-700 dark:bg-slate-900">
                {/* Scientific plots retain their native aspect ratio and load only near the viewport. */}
                {/* eslint-disable-next-line @next/next/no-img-element */}
                <img
                  src={src}
                  alt={`Generated plot for ${report.title}`}
                  loading="lazy"
                  decoding="async"
                  className="h-auto w-full"
                />
              </figure>
            ))}
          </div>
        </section>
      )}
    </div>
  );
}
