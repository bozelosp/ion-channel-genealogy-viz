import {
  closeSync,
  constants,
  fstatSync,
  lstatSync,
  openSync,
  readSync,
} from 'node:fs';
import { TextDecoder } from 'node:util';

export const MAX_SUPERMODEL_JSON_BYTES = 32 * 1024 * 1024;
export const MAX_SUPERMODEL_JSON_CONTAINER_ENTRIES = 100_000;
export const MAX_SUPERMODEL_JSON_DEPTH = 100;
export const MAX_SUPERMODEL_JSON_VALUES = 4_000_000;

const STRICT_UTF8_DECODER = new TextDecoder('utf-8', { fatal: true });

function sourceError(message) {
  const error = new Error(`Invalid supermodel JSON source: ${message}`);
  error.code = 'ERR_INVALID_SUPERMODEL_SOURCE';
  return error;
}

function assertJsonResourceBounds(sourceText) {
  let values = 0;
  let inString = false;
  let escaped = false;
  const containers = [];

  const chargeValue = () => {
    values += 1;
    if (values > MAX_SUPERMODEL_JSON_VALUES) {
      throw sourceError(`source exceeds the ${MAX_SUPERMODEL_JSON_VALUES.toLocaleString('en-US')}-value limit`);
    }
  };

  for (let index = 0; index < sourceText.length; index += 1) {
    const character = sourceText[index];
    if (inString) {
      if (escaped) {
        escaped = false;
      } else if (character === '\\') {
        escaped = true;
      } else if (character === '"') {
        inString = false;
      }
      continue;
    }

    if (character === '"') {
      inString = true;
      chargeValue();
      continue;
    }
    if (character === '{' || character === '[') {
      containers.push({ opening: character, commas: 0 });
      if (containers.length > MAX_SUPERMODEL_JSON_DEPTH) {
        throw sourceError(`source exceeds the nesting-depth limit of ${MAX_SUPERMODEL_JSON_DEPTH}`);
      }
      chargeValue();
      continue;
    }
    if (character === '}' || character === ']') {
      const container = containers.pop();
      const expectedOpening = character === '}' ? '{' : '[';
      if (!container || container.opening !== expectedOpening) {
        throw sourceError('source has unbalanced JSON containers');
      }
      continue;
    }
    if (character === ',') {
      const container = containers[containers.length - 1];
      if (container) {
        container.commas += 1;
        if (container.commas + 1 > MAX_SUPERMODEL_JSON_CONTAINER_ENTRIES) {
          throw sourceError(
            `source exceeds the ${MAX_SUPERMODEL_JSON_CONTAINER_ENTRIES.toLocaleString('en-US')}-entry container limit`,
          );
        }
      }
      continue;
    }
    if (
      character === '-'
      || (character >= '0' && character <= '9')
      || character === 't'
      || character === 'f'
      || character === 'n'
    ) {
      chargeValue();
      while (index + 1 < sourceText.length && !/[\s,\]}:]/.test(sourceText[index + 1])) {
        index += 1;
      }
    }
  }

  if (containers.length !== 0) {
    throw sourceError('source has unbalanced JSON containers');
  }
}

export function loadSupermodelJson(sourcePath) {
  const sourceInfo = lstatSync(sourcePath);
  if (sourceInfo.isSymbolicLink() || !sourceInfo.isFile()) {
    throw sourceError('source must be a regular, non-symlink file');
  }
  if (sourceInfo.size > MAX_SUPERMODEL_JSON_BYTES) {
    throw sourceError('source exceeds the 32 MiB limit');
  }
  if (typeof constants.O_NOFOLLOW !== 'number') {
    throw sourceError('this platform does not support no-follow file opens');
  }

  const descriptor = openSync(sourcePath, constants.O_RDONLY | constants.O_NOFOLLOW);
  try {
    const openedInfo = fstatSync(descriptor);
    if (!openedInfo.isFile()) {
      throw sourceError('opened source is not a regular file');
    }
    if (openedInfo.size > MAX_SUPERMODEL_JSON_BYTES) {
      throw sourceError('opened source exceeds the 32 MiB limit');
    }

    const contents = Buffer.allocUnsafe(openedInfo.size);
    let offset = 0;
    while (offset < contents.length) {
      const bytesRead = readSync(descriptor, contents, offset, contents.length - offset, null);
      if (bytesRead === 0) {
        throw sourceError('source changed while it was being read');
      }
      offset += bytesRead;
    }
    if (readSync(descriptor, Buffer.allocUnsafe(1), 0, 1, null) !== 0) {
      throw sourceError('source changed while it was being read');
    }

    let sourceText;
    try {
      sourceText = STRICT_UTF8_DECODER.decode(contents);
    } catch {
      throw sourceError('source is not valid UTF-8');
    }
    assertJsonResourceBounds(sourceText);

    try {
      return JSON.parse(sourceText, (_key, value) => {
        if (typeof value === 'number' && !Number.isFinite(value)) {
          throw sourceError('all numbers must be finite');
        }
        return value;
      });
    } catch (error) {
      if (error?.code === 'ERR_INVALID_SUPERMODEL_SOURCE') {
        throw error;
      }
      throw sourceError(`source is not valid JSON: ${error.message}`);
    }
  } finally {
    closeSync(descriptor);
  }
}
