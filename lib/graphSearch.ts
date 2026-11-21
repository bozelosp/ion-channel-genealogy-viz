interface SearchableNode {
  id: string | number;
  label?: string;
  name?: string;
  original_model?: Record<string, unknown>;
  identical_models?: Array<Record<string, unknown>>;
}

export function searchGraph(
  nodes: readonly SearchableNode[], visibleIds: ReadonlySet<string>, query: string,
): Set<string> {
  const matches = new Set<string>();
  const term = query.trim().toLowerCase();
  if (term.length < 3) return matches;

  const filenameOnly = term.endsWith('.mod') && !term.includes('/');
  const contains = (value: unknown) => typeof value === 'string' && value.toLowerCase().includes(term);
  function matchesModel(model?: Record<string, unknown>): boolean {
    if (!model) return false;
    if (filenameOnly) {
      return typeof model.mod_filename === 'string' && model.mod_filename.toLowerCase() === term;
    }
    return contains(model.unique_modelDB_mod_id) || contains(model.mod_filename) || contains(model.mod_filepath);
  }

  for (const node of nodes) {
    const id = String(node.id);
    if (!visibleIds.has(id)) continue;
    const matchesNode = !filenameOnly && (contains(id) || contains(node.label) || contains(node.name));
    if (matchesNode || matchesModel(node.original_model) ||
        (Array.isArray(node.identical_models) && node.identical_models.some(matchesModel))) {
      matches.add(id);
    }
  }
  return matches;
}
