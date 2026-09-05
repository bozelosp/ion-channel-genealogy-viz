import path from 'node:path';
import { promises as fs } from 'node:fs';

const UNSAFE_SEGMENT_CHARACTERS = /[\u0000-\u001f\u007f/\\]/;

export function safePathSegment(value, label) {
  const segment = String(value);
  if (
    !segment
    || segment === '.'
    || segment === '..'
    || segment.length > 200
    || UNSAFE_SEGMENT_CHARACTERS.test(segment)
  ) {
    throw new Error(`${label} is not a safe path segment`);
  }
  return segment;
}

export function resolveWithin(root, ...segments) {
  const resolvedRoot = path.resolve(root);
  const resolved = path.resolve(resolvedRoot, ...segments);
  const relative = path.relative(resolvedRoot, resolved);
  if (relative === '..' || relative.startsWith(`..${path.sep}`) || path.isAbsolute(relative)) {
    throw new Error(`Resolved path escapes ${resolvedRoot}`);
  }
  return resolved;
}

async function assertNoSymlinkComponents(root, target, label) {
  const resolvedRoot = path.resolve(root);
  const resolvedTarget = resolveWithin(resolvedRoot, path.relative(resolvedRoot, path.resolve(target)));

  const rootInfo = await fs.lstat(resolvedRoot);
  if (!rootInfo.isDirectory() || rootInfo.isSymbolicLink()) {
    throw new Error(`Protected root must be a real directory, not a symlink`);
  }

  const relative = path.relative(resolvedRoot, resolvedTarget);
  let current = resolvedRoot;
  for (const segment of relative.split(path.sep)) {
    current = path.join(current, segment);
    try {
      const info = await fs.lstat(current);
      if (info.isSymbolicLink()) {
        throw new Error(`${label} must not traverse a symbolic link`);
      }
      if (current !== resolvedTarget && !info.isDirectory()) {
        throw new Error(`${label} has a non-directory path component`);
      }
    } catch (error) {
      if (error && typeof error === 'object' && error.code === 'ENOENT') {
        break;
      }
      throw error;
    }
  }
  return resolvedTarget;
}

export async function ensureRealDirectory(directoryPath, label) {
  const resolved = path.resolve(directoryPath);
  await fs.mkdir(resolved, { recursive: true });
  const info = await fs.lstat(resolved);
  if (!info.isDirectory() || info.isSymbolicLink()) {
    throw new Error(`${label} must be a real directory, not a symlink`);
  }
  return resolved;
}

export async function assertWritableWithin(root, target, label) {
  return assertNoSymlinkComponents(root, target, label);
}

export async function assertRemovableWithin(root, target, label) {
  const resolvedRoot = path.resolve(root);
  const resolvedTarget = await assertNoSymlinkComponents(root, target, label);
  if (resolvedTarget === resolvedRoot) {
    throw new Error(`${label} must not be the protected root itself`);
  }
  return resolvedTarget;
}
