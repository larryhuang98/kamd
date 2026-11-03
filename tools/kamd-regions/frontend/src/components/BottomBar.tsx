import { useState } from "react";
import type { ExportResult, Settings, SystemInfo, Validation } from "../types";

interface Props {
  validation: Validation | null;
  system: SystemInfo;
  settings: Settings;
  onSettings: (p: Partial<Settings>) => void;
  onExport: (engines: string[]) => Promise<ExportResult | null>;
  onSave: () => void;
  onLoadFile: (f: File) => void;
  onImport: (files: FileList) => void;
  outDir: string;
}

function download(name: string, text: string) {
  const url = URL.createObjectURL(new Blob([text], { type: "text/plain" }));
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function BottomBar(p: Props) {
  const [engines, setEngines] = useState<string[]>(["amber", "gromacs", "openmm"]);
  const [result, setResult] = useState<ExportResult | null>(null);
  const [showSettings, setShowSettings] = useState(false);
  const v = p.validation;
  const blocking = v ? v.errors.filter((e) => e.engines.some((x) => engines.includes(x))) : [];
  const st = p.settings;
  const num = (k: keyof Settings, label: string, step = "any") => (
    <label>
      <span>{label}</span>
      <input
        type="number"
        step={step}
        value={st[k] === null || st[k] === undefined ? "" : String(st[k])}
        onChange={(e) => p.onSettings({ [k]: e.target.value === "" ? null : Number(e.target.value) } as Partial<Settings>)}
      />
    </label>
  );
  return (
    <footer className="panel bottom" data-testid="bottom-bar">
      <div className="validation" data-testid="validation">
        {!v && <span className="muted">validating…</span>}
        {v && !v.errors.length && !v.warnings.length && <div className="ok">✔ No problems found.</div>}
        {v?.errors.map((e, i) => (
          <div key={`e${i}`} className="issue error">
            <b>error</b> <span className="eng">{e.engines.join("/")}</span> {e.message}
          </div>
        ))}
        {v?.warnings.map((e, i) => (
          <div key={`w${i}`} className="issue warning">
            <b>warning</b> <span className="eng">{e.engines.join("/")}</span> {e.message}
          </div>
        ))}
        {v && v.n_unassigned_protein > 0 && (
          <div className="issue info">{v.n_unassigned_protein} protein residues are unassigned (KAMD does not move them).</div>
        )}
      </div>
      <div className="export">
        <div className="row">
          {["amber", "gromacs", "openmm"].map((e) => (
            <label key={e} className="chk">
              <input
                type="checkbox"
                checked={engines.includes(e)}
                onChange={(x) => setEngines(x.target.checked ? [...engines, e] : engines.filter((y) => y !== e))}
              />
              {e === "amber" ? "Amber" : e === "gromacs" ? "GROMACS" : "OpenMM"}
            </label>
          ))}
          <button
            className="primary"
            data-testid="export"
            disabled={!v || blocking.length > 0 || !engines.length}
            title={blocking.length ? "fix the errors first" : "export and verify"}
            onClick={async () => setResult(await p.onExport(engines))}
          >
            Export
          </button>
          <button onClick={p.onSave} data-testid="save">Save project</button>
          <label className="filebtn">
            Load project
            <input type="file" accept=".json" onChange={(e) => e.target.files?.[0] && p.onLoadFile(e.target.files[0])} />
          </label>
          <label className="filebtn" title="&kamd mdin, or kamd.dat + index.ndx, or kamd_openmm.json">
            Import KAMD input
            <input type="file" multiple onChange={(e) => e.target.files && p.onImport(e.target.files)} data-testid="import" />
          </label>
          <button className="ghost" onClick={() => setShowSettings((x) => !x)}>KAMD settings</button>
        </div>
        {showSettings && (
          <div className="params settings">
            {num("pseudobd_cycles", "pseudobd_cycles", "1")}
            {num("nzero_crossings", "nzero_crossings", "1")}
            {num("temperature", "temperature (K)")}
            {num("shell_cut", "shell_cut (Å)")}
            {num("surface_cut", "surface_cut (Å)")}
            {num("membrane_cut", "membrane_cut (Å)")}
            {num("kamd_shell_mode", "kamd_shell_mode (Amber)", "1")}
            {num("seed", "seed", "1")}
            <label>
              <span>membrane slab</span>
              <span className="slab">
                <select
                  value={st.membrane?.axis ?? 0}
                  onChange={(e) => {
                    const ax = Number(e.target.value);
                    const m = p.system.membrane;
                    p.onSettings({ membrane: ax ? { axis: ax, lo: st.membrane?.lo ?? m?.lo ?? -15, hi: st.membrane?.hi ?? m?.hi ?? 15 } : null });
                  }}
                >
                  <option value={0}>auto{p.system.membrane ? ` (${p.system.membrane.axis_name}: ${p.system.membrane.lo.toFixed(1)} … ${p.system.membrane.hi.toFixed(1)})` : " (none)"}</option>
                  <option value={1}>x</option>
                  <option value={2}>y</option>
                  <option value={3}>z</option>
                </select>
                {st.membrane && (
                  <>
                    <input type="number" step="0.5" value={st.membrane.lo} onChange={(e) => p.onSettings({ membrane: { ...st.membrane!, lo: Number(e.target.value) } })} />
                    <input type="number" step="0.5" value={st.membrane.hi} onChange={(e) => p.onSettings({ membrane: { ...st.membrane!, hi: Number(e.target.value) } })} />
                  </>
                )}
              </span>
            </label>
          </div>
        )}
        {result && (
          <div className="result" data-testid="export-result">
            {!result.ok && <div className="issue error">Export blocked for: {result.blocked.join(", ")}</div>}
            {result.verification && (
              <div className={result.verification.ok ? "ok" : "issue error"}>
                {result.verification.ok ? "✔ verified: " : "✖ verification failed: "}
                {result.verification.kamd_config.message}
                {(["amber", "gromacs", "openmm"] as const).map((e) =>
                  result.verification![e].length ? (
                    <span key={e}>
                      {" · "}
                      {e} {result.verification![e].filter((x) => x.ok).length}/{result.verification![e].length} match
                    </span>
                  ) : null,
                )}
              </div>
            )}
            <div className="files">
              {Object.keys(result.files).map((f) => (
                <button key={f} className="link" onClick={() => download(f, result.files[f])} data-testid={`file-${f}`}>
                  ⬇ {f}
                </button>
              ))}
            </div>
            {result.files["kamd_amber.in"] && <pre className="preview" data-testid="amber-preview">{result.files["kamd_amber.in"]}</pre>}
          </div>
        )}
        <div className="muted small">Saving writes to {p.outDir} on the machine running kamd-regions.</div>
      </div>
    </footer>
  );
}

export { download };
