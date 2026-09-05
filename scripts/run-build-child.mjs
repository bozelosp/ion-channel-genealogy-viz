// Register before loading any generator/build code. If the coordinator dies
// between spawn and registration, either our live record prevents takeover,
// or the missing unique owner record stops us before we can touch build data.
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';

const [ownerPath, token, script, ...args] = process.argv.slice(2);
const record = path.join(path.dirname(ownerPath), `child-${process.pid}-${token}`);
// The parent removes this record only after 'close'. If the parent dies,
// recovery removes it only after the OS reports this PID is no longer alive.

try {
  if (!fs.existsSync(ownerPath)) throw new Error('Build coordinator no longer owns the lock');
  fs.writeFileSync(record, '', { flag: 'wx', mode: 0o600 });
  if (!fs.existsSync(ownerPath)) throw new Error('Build coordinator no longer owns the lock');
  process.argv = [process.execPath, script, ...args];
  await import(pathToFileURL(script).href);
} catch (error) {
  console.error(`build child: ${error.message}`);
  process.exitCode = 1;
}
