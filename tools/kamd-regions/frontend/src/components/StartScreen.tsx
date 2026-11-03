import { useState } from "react";
import { bridge } from "../desktop";

interface Props {
  version: string;
  recent: string[][];
  error: string | null;
  busy: boolean;
  onOpen: (paths: string[]) => void;
  onCancel?: () => void;
}

const base = (p: string) => p.split(/[\\/]/).pop() ?? p;

export function StartScreen({ version, recent, error, busy, onOpen, onCancel }: Props) {
  const [text, setText] = useState("");
  const native = bridge();
  const pick = async () => {
    const files = await native!.open_dialog(true);
    if (files.length) onOpen(files);
  };
  return (
    <div className="start" data-testid="start-screen">
      <div className="start-card">
        <div className="start-head">
          <h1>KAMD Regions</h1>
          <span className="muted small">version {version}</span>
        </div>
        <p className="muted">
          Choose reference groups for Kinetic Accelerated MD: open a structure, review the suggested
          regions, edit them residue by residue, and export Amber, GROMACS or OpenMM input.
        </p>
        {native ? (
          <button className="primary big" onClick={pick} disabled={busy} data-testid="open-files">
            {busy ? "Loading…" : "Open structure…"}
          </button>
        ) : (
          <div className="path-input">
            <label htmlFor="paths">File paths on this computer, one per line</label>
            <textarea
              id="paths"
              rows={3}
              value={text}
              placeholder={"/path/system.prmtop\n/path/system.rst7"}
              onChange={(e) => setText(e.target.value)}
              data-testid="paths"
            />
            <button
              className="primary"
              disabled={busy || !text.trim()}
              onClick={() => onOpen(text.split("\n").map((s) => s.trim()).filter(Boolean))}
              data-testid="open-paths"
            >
              {busy ? "Loading…" : "Open"}
            </button>
          </div>
        )}
        <div className="muted small">
          Select several files at once: a structure (PDB, mmCIF, gro), or a topology with coordinates (prmtop + rst7/inpcrd,
          tpr + gro). Add an existing KAMD input (&amp;kamd mdin, kamd.dat + index.ndx) or a saved project.json to continue
          from it. A project made on another structure is moved over by sequence alignment.
        </div>
        {error && <div className="issue error" data-testid="open-error">{error}</div>}
        {recent.length > 0 && (
          <div className="recent">
            <div className="label">Recent</div>
            {recent.map((r, i) => (
              <button key={i} className="recent-item" onClick={() => onOpen(r)} disabled={busy} title={r.join("\n")}>
                <span>{r.map(base).join(" + ")}</span>
                <span className="muted small">{r[0].slice(0, Math.max(0, r[0].length - base(r[0]).length - 1))}</span>
              </button>
            ))}
          </div>
        )}
        {onCancel && (
          <button className="ghost" onClick={onCancel}>
            Back to the current system
          </button>
        )}
      </div>
    </div>
  );
}
