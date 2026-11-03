import { useCallback, useEffect, useMemo, useReducer, useRef, useState } from "react";
import { api, isDemo, type AppInfo } from "./api";
import { ErrorBoundary } from "./components/ErrorBoundary";
import { StartScreen } from "./components/StartScreen";
import { bridge, saveText, whenBridge } from "./desktop";
import { unassignedColor } from "./colors";
import { BottomBar } from "./components/BottomBar";
import { RegionPanel } from "./components/RegionPanel";
import { SequenceBar, type Row } from "./components/SequenceBar";
import { Viewer3D } from "./components/Viewer3D";
import { decode, encode, historyReducer, residueLabel, type Base, type Doc, type History } from "./state";
import type { ProjectDoc, StructureData, SystemInfo, Validation } from "./types";
import type { PickEvent } from "./viewer/MolViewer";

const emptyDoc: Doc = {
  regions: [],
  assign: new Int32Array(0),
  settings: {} as Doc["settings"],
  method: {},
  edited: false,
};

function useTheme(): [boolean, (d: boolean) => void] {
  const media = window.matchMedia?.("(prefers-color-scheme: dark)");
  const [dark, setDark] = useState<boolean>(() => {
    try {
      const s = localStorage.getItem("kamd-theme");
      if (s) return s === "dark";
    } catch {
      /* storage unavailable */
    }
    return !!media?.matches;
  });
  useEffect(() => {
    document.documentElement.dataset.theme = dark ? "dark" : "light";
    try {
      localStorage.setItem("kamd-theme", dark ? "dark" : "light");
    } catch {
      /* ignore */
    }
  }, [dark]);
  return [dark, setDark];
}

export default function App() {
  const [system, setSystem] = useState<SystemInfo | null>(null);
  const [structure, setStructure] = useState<StructureData | null>(null);
  const [base, setBase] = useState<Base>({ inputs: {}, n_residues: 0, residue_names: [] });
  const [hist, dispatch] = useReducer(historyReducer, { past: [], present: emptyDoc, future: [] } as History);
  const doc = hist.present;
  const [selected, setSelected] = useState<Uint8Array>(new Uint8Array(0));
  const [anchor, setAnchor] = useState<number | null>(null);
  const [preview, setPreview] = useState<Int32Array | null>(null);
  const [validation, setValidation] = useState<Validation | null>(null);
  const [notes, setNotes] = useState<string[]>([]);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [menu, setMenu] = useState<{ x: number; y: number } | null>(null);
  const [hover, setHover] = useState<{ pos: number; x: number; y: number } | null>(null);
  const [hover3d, setHover3d] = useState<string | null>(null);
  const [focusPos, setFocusPos] = useState<number | null>(null);
  const [outDir, setOutDir] = useState("");
  const [dark, setDark] = useTheme();

  const [appInfo, setAppInfo] = useState<AppInfo | null>(null);
  const [showStart, setShowStart] = useState(false);
  const [opening, setOpening] = useState(false);
  const [openError, setOpenError] = useState<string | null>(null);
  const [native, setNative] = useState(false);

  // ------------------------------------------------------------------ load
  useEffect(() => {
    whenBridge().then((b) => setNative(!!b));
    api
      .app()
      .then((info) => {
        setAppInfo(info);
        if (info.error) setOpenError(info.error);
        if (info.loaded) loadSession();
      })
      .catch((e) => setError(String(e)));
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const openPaths = async (paths: string[]) => {
    setOpening(true);
    setOpenError(null);
    try {
      await api.open(paths);
      window.location.reload();
    } catch (e) {
      setOpenError(String(e).replace(/^Error: /, ""));
      setOpening(false);
    }
  };

  const loadSession = () => {
    api
      .session()
      .then((s) => {
        setSystem(s.system);
        setBase({ inputs: s.project.inputs, n_residues: s.project.n_residues, residue_names: s.project.residue_names, residue_table: s.project.residue_table });
        dispatch({ type: "replace", doc: decode(s.project, s.system.residues) });
        setSelected(new Uint8Array(s.system.residues.index.length));
        setValidation(s.validation);
        setNotes(s.system.notes);
        setOutDir(s.out_dir);
        return api.structure();
      })
      .then(setStructure)
      .catch((e) => setError(String(e)));
  };

  const res = system?.residues;
  const posOf = useMemo(() => {
    const m = new Map<number, number>();
    res?.index.forEach((g, i) => m.set(g, i));
    return m;
  }, [res]);

  const rows: Row[] = useMemo(() => {
    if (!system || !res) return [];
    return system.chains.map((c) => {
      const positions: number[] = [];
      for (let g = c.start; g < c.stop; g++) {
        const p = posOf.get(g);
        if (p !== undefined) positions.push(p);
      }
      return { chain: c.chain, positions };
    });
  }, [system, res, posOf]);

  const assign = preview ?? doc.assign;
  const colors = useMemo(() => new Map(doc.regions.map((r) => [r.id, r.color])), [doc.regions]);

  // ------------------------------------------------------------------ validation (debounced)
  const projectDoc = useCallback((d: Doc = doc): ProjectDoc | null => (res ? encode(d, res, base) : null), [doc, res, base]);
  const vTimer = useRef<number>();
  useEffect(() => {
    if (!res || !doc.regions) return;
    window.clearTimeout(vTimer.current);
    vTimer.current = window.setTimeout(() => {
      const pd = projectDoc();
      if (pd)
        api
          .validate(pd)
          .then(setValidation)
          .catch((e) => setError(String(e)));
    }, 250);
  }, [doc, res, projectDoc]);

  // ------------------------------------------------------------------ selection
  const selList = useMemo(() => {
    const out: number[] = [];
    selected.forEach((v, i) => {
      if (v) out.push(i);
    });
    return out;
  }, [selected]);
  const selGlobal = useMemo(() => new Set(selList.map((p) => res!.index[p])), [selList, res]);

  const select = useCallback(
    (positions: number[], mode: "replace" | "add" | "toggle") => {
      setSelected((s) => {
        const n = mode === "replace" ? new Uint8Array(s.length) : s.slice();
        for (const p of positions) n[p] = mode === "toggle" ? (n[p] ? 0 : 1) : 1;
        return n;
      });
    },
    [],
  );

  const onPick = useCallback(
    (e: PickEvent) => {
      const pos = posOf.get(e.residue);
      if (pos === undefined) return;
      if (e.shift && anchor !== null) {
        const [a, b] = [Math.min(anchor, pos), Math.max(anchor, pos)];
        select(Array.from({ length: b - a + 1 }, (_, k) => a + k), e.ctrl ? "add" : "replace");
      } else if (e.ctrl) {
        select([pos], "toggle");
        setAnchor(pos);
      } else {
        select([pos], "replace");
        setAnchor(pos);
      }
      setFocusPos(pos);
    },
    [posOf, anchor, select],
  );

  const onHover3d = useCallback(
    (g: number | null) => {
      if (g === null || !res) return setHover3d(null);
      const pos = posOf.get(g);
      setHover3d(pos === undefined ? null : residueLabel(res, pos, { ...doc, assign }));
    },
    [posOf, res, doc, assign],
  );

  const colorOf = useCallback(
    (g: number) => {
      const p = posOf.get(g);
      const v = p === undefined ? 0 : assign[p];
      return v ? colors.get(v) ?? unassignedColor(dark) : unassignedColor(dark);
    },
    [posOf, assign, colors, dark],
  );

  // ------------------------------------------------------------------ commands
  const assignSel = useCallback(
    (region: number) => {
      if (!selList.length) return;
      dispatch({ type: "assign", positions: selList, region });
    },
    [selList],
  );

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      const t = e.target as HTMLElement;
      if (t && (t.tagName === "INPUT" || t.tagName === "SELECT" || t.tagName === "TEXTAREA")) return;
      const mod = e.ctrlKey || e.metaKey;
      if (mod && e.key.toLowerCase() === "z") {
        e.preventDefault();
        dispatch({ type: e.shiftKey ? "redo" : "undo" });
      } else if (mod && e.key.toLowerCase() === "y") {
        e.preventDefault();
        dispatch({ type: "redo" });
      } else if (!mod && /^[1-9]$/.test(e.key)) {
        const k = Number(e.key) - 1;
        if (k < doc.regions.length) assignSel(doc.regions[k].id);
        else if (k === doc.regions.length && selList.length) dispatch({ type: "newRegion", positions: selList });
      } else if (!mod && (e.key === "0" || e.key === "Delete" || e.key === "Backspace")) {
        assignSel(0);
      } else if (!mod && e.key === "n" && selList.length) {
        dispatch({ type: "newRegion", positions: selList });
      } else if (e.key === "Escape") {
        setSelected((s) => new Uint8Array(s.length));
        setMenu(null);
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [doc.regions, assignSel, selList]);

  const applyServerProject = (p: ProjectDoc, v: Validation, n: string[], keepHistory = true) => {
    if (!res) return;
    dispatch({ type: "replace", doc: decode(p, res), keepHistory });
    setValidation(v);
    if (n.length) setNotes(n);
  };

  const onSuggest = async (method: string, params: Record<string, unknown>) => {
    setBusy(true);
    try {
      const r = await api.suggest(method, params, doc.settings);
      applyServerProject(r.project, r.validation, r.notes.length ? r.notes : [`${method}: ${r.project.regions.length} regions suggested`]);
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  };

  const onExport = async (engines: string[]) => {
    const pd = projectDoc();
    if (!pd) return null;
    try {
      return await api.export(pd, engines);
    } catch (e) {
      setError(String(e));
      return null;
    }
  };

  const onSave = async () => {
    const pd = projectDoc();
    if (!pd) return;
    const b = bridge();
    if (b) {
      const path = await b.save_dialog("project.json");
      if (!path) return;
      try {
        const r = await api.save(pd, path);
        setNotes([`project saved to ${r.path}`]);
      } catch (e) {
        setError(String(e));
      }
      return;
    }
    try {
      const r = await api.save(pd);
      setNotes([isDemo ? "project.json downloaded" : `project saved to ${r.path}`]);
      saveText("project.json", JSON.stringify(pd, null, 1));
    } catch (e) {
      setError(String(e));
    }
  };

  const onLoadFile = async (f: File) => {
    try {
      const r = await api.loadProject(JSON.parse(await f.text()));
      applyServerProject(r.project, r.validation, r.notes.length ? r.notes : [`loaded ${f.name}`]);
    } catch (e) {
      setError(String(e));
    }
  };

  const onImport = async (files: FileList) => {
    try {
      const list = Array.from(files);
      const texts = await Promise.all(list.map((f) => f.text()));
      let kind = "";
      let text = "";
      let ndx: string | undefined;
      list.forEach((f, i) => {
        const t = texts[i];
        if (/&kamd/i.test(t)) [kind, text] = ["amber", t];
        else if (f.name.endsWith(".ndx")) ndx = t;
        else if (f.name.endsWith(".json")) [kind, text] = ["openmm", t];
        else if (/refer[-_]groups/.test(t)) [kind, text] = ["gromacs", t];
      });
      if (!kind) throw new Error("no &kamd namelist, kamd.dat or kamd_openmm.json among the files");
      if (kind === "gromacs" && !ndx) throw new Error("select kamd.dat together with its index.ndx");
      const r = await api.importText(kind, text, ndx);
      applyServerProject(r.project, r.validation, r.notes.length ? r.notes : ["imported"]);
    } catch (e) {
      setError(String(e));
    }
  };

  const onExportFolder = async (engines: string[]) => {
    const b = bridge();
    const pd = projectDoc();
    if (!b || !pd) return null;
    const dir = await b.folder_dialog();
    if (!dir) return null;
    try {
      return await api.export(pd, engines, true, dir);
    } catch (e) {
      setError(String(e));
      return null;
    }
  };

  // menu bar of the desktop app (and keyboard-free access to the same commands)
  useEffect(() => {
    (window as unknown as { kamdMenu?: (c: string) => void }).kamdMenu = (cmd: string) => {
      if (cmd === "open") setShowStart(true);
      else if (cmd === "save") onSave();
      else if (cmd === "load") document.getElementById("load-project")?.click();
      else if (cmd === "import") document.getElementById("import-input")?.click();
      else if (cmd === "export") document.getElementById("export-folder")?.click();
      else if (cmd === "undo") dispatch({ type: "undo" });
      else if (cmd === "redo") dispatch({ type: "redo" });
      else if (cmd === "theme") setDark(!dark);
    };
  });

  // state for the desktop self-test and the e2e tests
  useEffect(() => {
    (window as unknown as { __kamdState?: unknown }).__kamdState = appInfo
      ? {
          screen: !appInfo.loaded || showStart ? "start" : system ? "main" : "loading",
          desktop: native,
          regions: doc.regions.length,
          atoms: system?.n_atoms ?? 0,
        }
      : null;
  }, [appInfo, showStart, system, doc.regions.length, native]);

  // ------------------------------------------------------------------ render
  if (appInfo && (!appInfo.loaded || showStart))
    return (
      <StartScreen
        version={appInfo.version}
        recent={appInfo.recent}
        error={openError}
        busy={opening}
        onOpen={openPaths}
        onCancel={appInfo.loaded ? () => setShowStart(false) : undefined}
      />
    );
  if (error && !system) return <div className="fatal">Could not load: {error}</div>;
  if (!system || !res) return <div className="loading">Loading system…</div>;
  const membrane = doc.settings.membrane ?? system.membrane;

  return (
    <div className="app" onClick={() => menu && setMenu(null)}>
      <header className="top">
        <b>KAMD Regions</b>
        <span className="muted small">
          {Object.values(system.files).map((f) => f.split("/").pop()).join(" + ")} · {system.n_atoms.toLocaleString()} atoms ·{" "}
          {res.index.length} assignable residues
          {system.membrane ? ` · membrane ${system.membrane.axis_name} ${system.membrane.lo.toFixed(1)}…${system.membrane.hi.toFixed(1)} Å (${system.membrane.source})` : ""}
          {!system.has_topology && <span className="warn-chip"> no topology: Amber numbering may differ</span>}
        </span>
        <span className="spacer" />
        {!isDemo && (
          <button className="ghost" onClick={() => setShowStart(true)} title="open another structure" data-testid="open-other">
            Open…
          </button>
        )}
        <button className="ghost" disabled={!hist.past.length} onClick={() => dispatch({ type: "undo" })} title="undo (ctrl+z)">↶</button>
        <button className="ghost" disabled={!hist.future.length} onClick={() => dispatch({ type: "redo" })} title="redo (ctrl+shift+z)">↷</button>
        <button className="ghost" onClick={() => setDark(!dark)} title="toggle theme">{dark ? "☀" : "☾"}</button>
      </header>
      {isDemo && (
        <div className="demo-banner" data-testid="demo-banner">
          Offline snapshot: suggestions were computed by kamd-regions and are embedded in this page. Selecting,
          editing, undo and validation work here; switching methods shows the precomputed default results.
          Exporting Amber/GROMACS/OpenMM inputs needs the local program (<code>kamd-regions serve</code>).
        </div>
      )}
      {(notes.length > 0 || error) && (
        <div className="notes">
          {error && <div className="issue error">{error}</div>}
          {notes.map((n, i) => (
            <div key={i}>{n}</div>
          ))}
          <button className="ghost tiny" onClick={() => { setNotes([]); setError(null); }}>dismiss</button>
        </div>
      )}
      <div className="main">
        <div className="center">
          <ErrorBoundary label="The 3D view">
          <Viewer3D
            structure={structure}
            res={res}
            colorOf={colorOf}
            selected={selGlobal}
            membrane={membrane}
            dark={dark}
            onPick={onPick}
            onHover={onHover3d}
            hoverLabel={hover3d}
          />
          </ErrorBoundary>
          <SequenceBar
            res={res}
            rows={rows}
            assign={assign}
            colors={colors}
            selected={selected}
            dark={dark}
            focusPos={focusPos}
            onSelect={select}
            anchor={anchor}
            onAnchor={setAnchor}
            onBoundaryPreview={setPreview}
            onBoundaryCommit={(a) => dispatch({ type: "setAssign", assign: a })}
            onContext={(x, y) => setMenu({ x, y })}
            onHover={(pos, x, y) => setHover(pos === null ? null : { pos, x, y })}
          />
        </div>
        <RegionPanel
          doc={doc}
          res={res}
          system={system}
          validation={validation}
          selection={selList}
          dispatch={dispatch}
          onSelectRegion={(id) => {
            const ps: number[] = [];
            doc.assign.forEach((v, i) => v === id && ps.push(i));
            select(ps, "replace");
            if (ps.length) setFocusPos(ps[0]);
          }}
          onSuggest={onSuggest}
          busy={busy}
        />
      </div>
      <BottomBar
        validation={validation}
        system={system}
        settings={doc.settings}
        onSettings={(patch) => dispatch({ type: "settings", patch })}
        onExport={onExport}
        onSave={onSave}
        onExportFolder={native ? onExportFolder : undefined}
        onLoadFile={onLoadFile}
        onImport={onImport}
        outDir={outDir}
      />
      {hover && (
        <div className="tooltip" style={{ left: hover.x + 12, top: hover.y + 14 }}>
          {residueLabel(res, hover.pos, { ...doc, assign })}
        </div>
      )}
      {menu && (
        <div className="menu" style={{ left: menu.x, top: menu.y }} data-testid="context-menu" onClick={(e) => e.stopPropagation()}>
          <div className="muted small">{selList.length} residue(s)</div>
          {doc.regions.map((r, k) => (
            <button key={r.id} onClick={() => { assignSel(r.id); setMenu(null); }}>
              <span className="dot" style={{ background: r.color }} /> {k + 1}. {r.name}
            </button>
          ))}
          <button onClick={() => { dispatch({ type: "newRegion", positions: selList }); setMenu(null); }}>+ new region</button>
          <button onClick={() => { assignSel(0); setMenu(null); }}>unassign</button>
        </div>
      )}
    </div>
  );
}
