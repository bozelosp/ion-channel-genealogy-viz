export const OTHER_ION_CLASS_LABEL = 'Other';

const BASE_ION_CLASS_OPTIONS = ['K', 'Na', 'Ca', 'Ih', 'KCa', OTHER_ION_CLASS_LABEL] as const;

/**
 * Controls whether the "Other" ion class should be exposed in the UI.
 * Defaults to `false` when `NEXT_PUBLIC_SHOW_OTHER_CLASS` is not set.
 */
export const SHOW_OTHER_ION_CLASS = process.env.NEXT_PUBLIC_SHOW_OTHER_CLASS === 'true';

export type IonClassName = (typeof BASE_ION_CLASS_OPTIONS)[number];

function normaliseIonClass(raw: unknown): string | null {
  if (typeof raw !== 'string') return null;
  const trimmed = raw.trim();
  return trimmed.length > 0 ? trimmed : null;
}

export function isOtherIonClass(value: unknown): boolean {
  const normalised = normaliseIonClass(value);
  if (!normalised) return false;
  return normalised.toLowerCase() === OTHER_ION_CLASS_LABEL.toLowerCase();
}

export function shouldSuppressIonClass(value?: string | null): boolean {
  return !SHOW_OTHER_ION_CLASS && isOtherIonClass(value);
}

export function formatIonClassForDisplay(value?: string | null): string {
  if (!value || shouldSuppressIonClass(value)) {
    return '—';
  }
  return value;
}

/**
 * Returns the list of ion classes that should be shown to the user, preserving the
 * canonical ordering used throughout the UI.
 */
export function getVisibleIonClassOptions(): string[] {
  if (SHOW_OTHER_ION_CLASS) {
    return Array.from(BASE_ION_CLASS_OPTIONS);
  }
  return Array.from(BASE_ION_CLASS_OPTIONS).filter((cls) => cls !== OTHER_ION_CLASS_LABEL);
}

