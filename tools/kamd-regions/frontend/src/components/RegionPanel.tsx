import { useState } from "react";
import { PALETTE } from "../colors";
import type { Doc } from "../state";
import { positionsOf, summary } from "../state";
import type { Residues, SystemInfo, Validation } from "../types";
import type { Action } from "../state";

interface Props {
  doc: Doc;
  res: Residues;
  system: SystemInfo;
  validation: Validation | null;
  selection: number[];
  dispatch: (a: Action) => void;
  onSelectRegion: (id: number) => void;
  onSuggest: (method: string, params: Record<string, unknown>) => void;
  busy: boolean;
}

export function RegionPanel({ doc, res, system, validation, selection, dispatch, onSelectRegion, onSuggest, busy }: Props) {
  const [checked, setChecked] = useState<Set<number>>(new Set());
  const [method, setMethod] = useState<string>(system.membrane ? "tm" : "gnm");
  const [params, setParams] = useState<Record<string, Record<string, unknown>>>({});
  const [showParams, setShowParams] = useState(false);
  const [confirming, setConfirming] = useState(false); // in-page confirm: browser dialogs may be blocked
  const [palette, setPalette] = useState<number | null>(null);
  const toggle = (id: number) =>
    setChecked((s) => {
      const n = new Set(s);
      if (n.has(id)) n.delete(id);
      else n.add(id);
      return n;
    });
  const methodInfo = system.methods[method] as { label: string; params: Record<string, number | boolean> } | undefined;
  const post = system.methods["postprocess"] as Record<string, number | boolean>;
  const allParams = { ...(methodInfo?.params ?? {}), ...post };
  const cur = params[method] ?? {};
  const selRegion = (() => {
    const ids = new Set(selection.map((p) => doc.assign[p]).filter((v) => v));
    return ids.size === 1 ? [...ids][0] : null;
  })();

  return (
    <aside className="panel regions" data-testid="region-panel">
      <div className="suggest">
        <div className="row">
          <select value={method} onChange={(e) => setMethod(e.target.value)} data-testid="method-select" aria-label="method">
            <option value="gnm">GNM slow modes</option>
            <option value="contacts">Contact map + Leiden</option>
            <option value="tm" disabled={!system.membrane}>Transmembrane (tm)</option>
          </select>
          <button
            className="primary"
            disabled={busy}
            data-testid="resuggest"
            title={confirming ? "your manual edits will be replaced (undo is available)" : "run the method"}
            onClick={() => {
              if (doc.edited && !confirming) {
                setConfirming(true);
                window.setTimeout(() => setConfirming(false), 4000);
                return;
              }
              setConfirming(false);
              onSuggest(method, cur);
            }}
          >
            {busy ? "Working…" : confirming ? "Replace edits?" : "Re-suggest"}
          </button>
          <button className="ghost" onClick={() => setShowParams((x) => !x)} title="method parameters">⚙</button>
        </div>
        {showParams && (
          <div className="params">
            {Object.entries(allParams).map(([k, v]) => (
              <label key={k}>
                <span>{k}</span>
                {typeof v === "boolean" ? (
                  <input
                    type="checkbox"
                    checked={(cur[k] as boolean | undefined) ?? v}
                    onChange={(e) => setParams({ ...params, [method]: { ...cur, [k]: e.target.checked } })}
                  />
                ) : k === "resolution" ? (
                  <span className="slider">
                    <input
                      type="range"
                      min={0.05}
                      max={2}
                      step={0.05}
                      value={(cur[k] as number | undefined) ?? v}
                      onChange={(e) => setParams({ ...params, [method]: { ...cur, [k]: Number(e.target.value) } })}
                    />
                    <span>{String((cur[k] as number | undefined) ?? v)}</span>
                  </span>
                ) : (
                  <input
                    type="number"
                    value={String((cur[k] as number | undefined) ?? v)}
                    step="any"
                    onChange={(e) => setParams({ ...params, [method]: { ...cur, [k]: Number(e.target.value) } })}
                  />
                )}
              </label>
            ))}
            <div className="muted small">n_regions 0 = automatic (eigengap); per_group 0 = by size.</div>
          </div>
        )}
      </div>

      <div className="region-tools">
        <button onClick={() => dispatch({ type: "newRegion", positions: selection })} title="new region from the selection (n)" data-testid="new-region">
          + New
        </button>
        <button disabled={checked.size < 2} onClick={() => { dispatch({ type: "mergeRegions", ids: [...checked] }); setChecked(new Set()); }}>
          Merge
        </button>
        <button
          disabled={selRegion === null}
          title="move the selected residues of one region into a new region"
          onClick={() => selRegion !== null && dispatch({ type: "splitRegion", id: selRegion, positions: selection })}
        >
          Split
        </button>
        <button disabled={!checked.size} onClick={() => { dispatch({ type: "deleteRegions", ids: [...checked] }); setChecked(new Set()); }}>
          Delete
        </button>
      </div>

      <div className="region-list">
        {doc.regions.map((r, k) => {
          const pos = positionsOf(doc, r.id);
          const nAtoms = pos.reduce((s, p) => s + res.n_atoms[p], 0);
          const st = validation?.regions[String(r.id)];
          const issues = st?.issues ?? [];
          const level = issues.some((i) => i.level === "error") ? "error" : issues.length ? "warning" : "";
          return (
            <div key={r.id} className={`region ${checked.has(r.id) ? "checked" : ""}`} data-testid={`region-${k + 1}`}>
              <div className="row">
                <input type="checkbox" checked={checked.has(r.id)} onChange={() => toggle(r.id)} aria-label="select region" />
                <button className="swatch" style={{ background: r.color }} onClick={() => setPalette(palette === r.id ? null : r.id)} title="change colour">
                  {k + 1}
                </button>
                <input
                  className="name"
                  value={r.name}
                  onChange={(e) => dispatch({ type: "updateRegion", id: r.id, patch: { name: e.target.value } })}
                  aria-label="region name"
                />
                {level && (
                  <span className={`badge ${level}`} title={issues.map((i) => i.message).join("\n")}>
                    {level === "error" ? "✖" : "⚠"}
                  </span>
                )}
                <button className="ghost tiny" title="move up (group order)" onClick={() => dispatch({ type: "moveRegion", id: r.id, delta: -1 })}>▲</button>
                <button className="ghost tiny" title="move down" onClick={() => dispatch({ type: "moveRegion", id: r.id, delta: 1 })}>▼</button>
              </div>
              {palette === r.id && (
                <div className="palette">
                  {PALETTE.map((c) => (
                    <button key={c} style={{ background: c }} onClick={() => { dispatch({ type: "updateRegion", id: r.id, patch: { color: c } }); setPalette(null); }} />
                  ))}
                  <input type="color" value={r.color} onChange={(e) => dispatch({ type: "updateRegion", id: r.id, patch: { color: e.target.value } })} />
                </div>
              )}
              <div className="row small">
                <a className="summary" onClick={() => onSelectRegion(r.id)} title="select these residues">
                  {summary(pos, res)}
                </a>
              </div>
              <div className="row small">
                <span className="atoms">{nAtoms.toLocaleString()} atoms</span>
                <label className="seg" title="region_type">
                  <button className={r.region_type === 1 ? "on" : ""} onClick={() => dispatch({ type: "updateRegion", id: r.id, patch: { region_type: 1 } })}>
                    solution
                  </button>
                  <button className={r.region_type === 2 ? "on" : ""} onClick={() => dispatch({ type: "updateRegion", id: r.id, patch: { region_type: 2 } })}>
                    TM
                  </button>
                </label>
                <label title="refer_link: share the random direction of another group">
                  link
                  <select value={r.refer_link} onChange={(e) => dispatch({ type: "updateRegion", id: r.id, patch: { refer_link: Number(e.target.value) } })}>
                    <option value={0}>own</option>
                    {doc.regions.map((q, j) => (q.id === r.id ? null : <option key={q.id} value={q.id}>→ {j + 1}</option>))}
                  </select>
                </label>
                {st?.frac_in_slab !== null && st?.frac_in_slab !== undefined && (
                  <span className="muted" title="fraction of residues inside the membrane slab">{Math.round(st.frac_in_slab * 100)}% in slab</span>
                )}
              </div>
            </div>
          );
        })}
        {!doc.regions.length && <div className="muted">No regions. Select residues and press “+ New”.</div>}
      </div>
    </aside>
  );
}
