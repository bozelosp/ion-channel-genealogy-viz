import assert from 'node:assert/strict';
import test from 'node:test';

test('changes-only and full views preserve source lines for identical and edited files', async (t) => {
  let reply;
  globalThis.self = { postMessage: (message) => { reply = message; } };
  t.after(() => { delete globalThis.self; });
  await import('../app/visualizer/diff.worker.ts');

  function compare(source, target, changesOnly) {
    self.onmessage({ data: { id: 1, source, target, changesOnly, theme: 'light' } });
    assert.equal(reply.error, undefined);
    return reply.html;
  }

  const source = 'common\nbefore\n';
  const identical = compare(source, source, true);
  assert.match(identical, /source files are identical/);
  assert.match(identical, /<tbody><\/tbody>/);
  assert.match(compare(source, source, false), /&nbsp;common</);

  const changed = compare(source, 'common\nafter\n', true);
  assert.doesNotMatch(changed, /&nbsp;common</);
  assert.match(changed, />2<\/td><td class="diff_sub">-&nbsp;before</);
  assert.match(changed, />2<\/td><td class="diff_add">\+&nbsp;after</);
});
