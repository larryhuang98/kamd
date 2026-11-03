import { useCallback, useEffect, useRef, useState } from "react";
import { textOn, unassignedColor } from "../colors";
import type { Residues } from "../types";

export interface Row {
  chain: string;
  positions: number[]; // assignable positions in this row, in sequence order
}

interface Props {
  res: Residues;
  rows: Row[];
  assign: Int32Array;
  colors: Map<number, string>;
  selected: Uint8Array;
  dark: boolean;
  focusPos: number | null;
  onSelect: (positions: number[], mode: "replace" | "add" | "toggle") => void;
  onAnchor: (pos: number) => void;
  anchor: number | null;
  onBoundaryPreview: (a: Int32Array | null) => void;
  onBoundaryCommit: (a: Int32Array) => void;
  onContext: (x: number, y: number) => void;
  onHover: (pos: number | null, x: number, y: number) => void;
}

const LABEL_W = 46;
const TM_H = 6;
const RULER_H = 13;
const CELL_H = 18;
const SS_H = 7;
const GAP = 10;
const ROW_H = TM_H + RULER_H + CELL_H + SS_H + GAP;

type Drag =
  | { kind: "select"; row: number; anchorJ: number; mode: "replace" | "add" }
  | { kind: "boundary"; row: number; j: number; left: number; right: number; orig: Int32Array; minJ: number; maxJ: number }
  | null;

export function SequenceBar(p: Props) {
  const wrap = useRef<HTMLDivElement>(null);
  const canvas = useRef<HTMLCanvasElement>(null);
  const [zoom, setZoom] = useState(1);
  const [width, setWidth] = useState(800);
  const [scroll, setScroll] = useState(0);
  const drag = useRef<Drag>(null);
  const [cursor, setCursor] = useState("default");
  const cw = Math.round(11 * zoom * 10) / 10;
  const maxLen = Math.max(1, ...p.rows.map((r) => r.positions.length));
  const totalW = LABEL_W + maxLen * cw + 20;
  const height = Math.max(ROW_H, p.rows.length * ROW_H);

  // row/column lookup for every position
  const where = useRef(new Map<number, [number, number]>());
  useEffect(() => {
    const m = new Map<number, [number, number]>();
    p.rows.forEach((r, i) => r.positions.forEach((pos, j) => m.set(pos, [i, j])));
    where.current = m;
  }, [p.rows]);

  useEffect(() => {
    const el = wrap.current;
    if (!el) return;
    const ro = new ResizeObserver(() => setWidth(el.clientWidth));
    ro.observe(el);
    setWidth(el.clientWidth);
    return () => ro.disconnect();
  }, []);

  // scroll a focused residue into view
  useEffect(() => {
    if (p.focusPos === null || !wrap.current) return;
    const w = where.current.get(p.focusPos);
    if (!w) return;
    const x = LABEL_W + w[1] * cw;
    const el = wrap.current;
    if (x - el.scrollLeft < LABEL_W + 10 || x - el.scrollLeft > el.clientWidth - 30) {
      el.scrollLeft = Math.max(0, x - el.clientWidth / 2);
    }
  }, [p.focusPos, cw]);

  const draw = useCallback(() => {
    const c = canvas.current;
    if (!c) return;
    const dpr = window.devicePixelRatio || 1;
    c.width = Math.floor(width * dpr);
    c.height = Math.floor(height * dpr);
    c.style.width = `${width}px`;
    c.style.height = `${height}px`;
    const g = c.getContext("2d")!;
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    const css = getComputedStyle(document.documentElement);
    const bg = css.getPropertyValue("--panel").trim() || (p.dark ? "#1e1f22" : "#fff");
    const fg = css.getPropertyValue("--fg").trim() || (p.dark ? "#ddd" : "#222");
    const muted = css.getPropertyValue("--muted").trim() || "#888";
    const tmColor = css.getPropertyValue("--tm").trim() || "#b8860b";
    const hl = p.dark ? "#ffffff" : "#000000";
    g.fillStyle = bg;
    g.fillRect(0, 0, width, height);
    const un = unassignedColor(p.dark);
    const j0 = Math.max(0, Math.floor(scroll / cw) - 1);
    const j1 = Math.ceil((scroll + width - LABEL_W) / cw) + 1;
    g.textBaseline = "middle";
    g.textAlign = "center";
    p.rows.forEach((row, i) => {
      const y0 = i * ROW_H;
      const yCell = y0 + TM_H + RULER_H;
      const end = Math.min(j1, row.positions.length);
      for (let j = j0; j < end; j++) {
        const pos = row.positions[j];
        const x = LABEL_W + j * cw - scroll;
        const v = p.assign[pos];
        const col = v ? p.colors.get(v) ?? un : un;
        g.fillStyle = col;
        g.fillRect(x, yCell, Math.ceil(cw), CELL_H);
        if (p.res.tm[pos]) {
          g.fillStyle = tmColor;
          g.fillRect(x, y0 + 1, Math.ceil(cw), TM_H - 2);
        }
        if (cw >= 7) {
          g.fillStyle = textOn(col);
          g.font = `${Math.min(13, cw + 1)}px ui-monospace, Menlo, Consolas, monospace`;
          g.fillText(p.res.one[pos], x + cw / 2, yCell + CELL_H / 2 + 0.5);
        }
        const num = p.res.number[pos];
        if (num % 10 === 0) {
          g.fillStyle = muted;
          g.fillRect(x + cw / 2, y0 + TM_H + RULER_H - 3, 1, 3);
          g.font = "10px system-ui, sans-serif";
          g.fillText(String(num), x + cw / 2, y0 + TM_H + 5);
        }
        // secondary structure
        const ss = p.res.ss[pos];
        const ySS = yCell + CELL_H + 1;
        g.fillStyle = fg;
        if (ss === "H") g.fillRect(x, ySS, Math.ceil(cw), SS_H - 1);
        else if (ss === "E") {
          const nxt = j + 1 < row.positions.length ? p.res.ss[row.positions[j + 1]] : "";
          if (nxt === "E") g.fillRect(x, ySS + 1, Math.ceil(cw), SS_H - 3);
          else {
            g.beginPath();
            g.moveTo(x, ySS);
            g.lineTo(x + cw, ySS + (SS_H - 1) / 2);
            g.lineTo(x, ySS + SS_H - 1);
            g.fill();
          }
        } else if (ss) g.fillRect(x, ySS + SS_H / 2 - 1, Math.ceil(cw), 1);
        if (p.selected[pos]) {
          g.fillStyle = hl;
          g.fillRect(x, yCell - 2, Math.ceil(cw), 2);
          g.fillRect(x, yCell + CELL_H, Math.ceil(cw), 2);
          g.globalAlpha = 0.18;
          g.fillRect(x, yCell, Math.ceil(cw), CELL_H);
          g.globalAlpha = 1;
        }
        // boundary handle between j and j+1
        if (j + 1 < row.positions.length) {
          const w = p.assign[row.positions[j + 1]];
          if (w !== v && (w || v)) {
            const bx = x + cw;
            g.fillStyle = hl;
            g.fillRect(bx - 1, yCell - 3, 2, CELL_H + 6);
            g.fillRect(bx - 3, yCell - 3, 6, 3);
            g.fillRect(bx - 3, yCell + CELL_H, 6, 3);
          }
        }
      }
      // chain label (drawn last so it covers scrolled cells)
      g.fillStyle = bg;
      g.fillRect(0, y0, LABEL_W, ROW_H);
      g.fillStyle = fg;
      g.textAlign = "left";
      g.font = "bold 12px system-ui, sans-serif";
      g.fillText(row.chain || "-", 6, yCell + CELL_H / 2);
      if (row.positions.some((q) => p.res.tm[q])) {
        g.fillStyle = tmColor;
        g.font = "9px system-ui, sans-serif";
        g.fillText("TM", 6, y0 + TM_H / 2 + 1);
      }
      g.textAlign = "center";
      if (i > 0) {
        g.fillStyle = muted;
        g.globalAlpha = 0.4;
        g.fillRect(0, y0 - GAP / 2, width, 1);
        g.globalAlpha = 1;
      }
    });
  }, [p.rows, p.assign, p.colors, p.selected, p.dark, p.res, width, height, scroll, cw]);

  useEffect(() => {
    draw();
  }, [draw]);

  // ------------------------------------------------------------------ hit testing
  const hit = (e: { clientX: number; clientY: number }) => {
    const r = canvas.current!.getBoundingClientRect();
    const x = e.clientX - r.left;
    const y = e.clientY - r.top;
    const row = Math.floor(y / ROW_H);
    const fx = (x - LABEL_W + scroll) / cw;
    const j = Math.floor(fx);
    const yIn = y - row * ROW_H;
    const inCells = yIn >= TM_H + RULER_H - 4 && yIn <= TM_H + RULER_H + CELL_H + 4;
    return { x, y, row, j, fx, inCells };
  };

  const boundaryAt = (row: number, fx: number): number | null => {
    const r = p.rows[row];
    if (!r) return null;
    const b = Math.round(fx); // boundary between b-1 and b
    if (Math.abs(fx - b) * cw > 4 || b <= 0 || b >= r.positions.length) return null;
    const a1 = p.assign[r.positions[b - 1]];
    const a2 = p.assign[r.positions[b]];
    return a1 !== a2 && (a1 || a2) ? b : null;
  };

  const runBounds = (row: number, b: number, assign: Int32Array) => {
    const pos = p.rows[row].positions;
    const L = assign[pos[b - 1]];
    const R = assign[pos[b]];
    let s = b - 1;
    while (s > 0 && assign[pos[s - 1]] === L) s--;
    let e = b;
    while (e + 1 < pos.length && assign[pos[e + 1]] === R) e++;
    return { L, R, minJ: s + 1, maxJ: e }; // new boundary in [s+1, e]
  };

  const onDown = (e: React.MouseEvent) => {
    if (e.button !== 0) return;
    const h = hit(e);
    const row = p.rows[h.row];
    if (!row) return;
    if (h.inCells) {
      const b = boundaryAt(h.row, h.fx);
      if (b !== null) {
        const rb = runBounds(h.row, b, p.assign);
        drag.current = { kind: "boundary", row: h.row, j: b, left: rb.L, right: rb.R, orig: p.assign.slice(), minJ: rb.minJ, maxJ: rb.maxJ };
        e.preventDefault();
        return;
      }
    }
    if (h.j < 0 || h.j >= row.positions.length) return;
    const pos = row.positions[h.j];
    const ctrl = e.ctrlKey || e.metaKey;
    if (e.shiftKey && p.anchor !== null) {
      const w = where.current.get(p.anchor);
      if (w && w[0] === h.row) {
        const [a, b2] = [Math.min(w[1], h.j), Math.max(w[1], h.j)];
        p.onSelect(row.positions.slice(a, b2 + 1), ctrl ? "add" : "replace");
        return;
      }
    }
    p.onAnchor(pos);
    if (ctrl) {
      p.onSelect([pos], "toggle");
      drag.current = { kind: "select", row: h.row, anchorJ: h.j, mode: "add" };
    } else {
      p.onSelect([pos], "replace");
      drag.current = { kind: "select", row: h.row, anchorJ: h.j, mode: "replace" };
    }
    e.preventDefault();
  };

  const onMove = (e: React.MouseEvent) => {
    const h = hit(e);
    const d = drag.current;
    if (!d) {
      const row = p.rows[h.row];
      const onB = h.inCells && boundaryAt(h.row, h.fx) !== null;
      setCursor(onB ? "ew-resize" : row && h.j >= 0 && h.j < row.positions.length ? "pointer" : "default");
      if (row && h.j >= 0 && h.j < row.positions.length && h.inCells) p.onHover(row.positions[h.j], e.clientX, e.clientY);
      else p.onHover(null, 0, 0);
      return;
    }
    const row = p.rows[d.row];
    if (d.kind === "select") {
      const j = Math.max(0, Math.min(row.positions.length - 1, h.j));
      const [a, b] = [Math.min(d.anchorJ, j), Math.max(d.anchorJ, j)];
      p.onSelect(row.positions.slice(a, b + 1), d.mode === "add" ? "add" : "replace");
    } else {
      const nb = Math.max(d.minJ, Math.min(d.maxJ, Math.round(h.fx)));
      const a = d.orig.slice();
      if (nb > d.j) for (let k = d.j; k < nb; k++) a[row.positions[k]] = d.left;
      else for (let k = nb; k < d.j; k++) a[row.positions[k]] = d.right;
      p.onBoundaryPreview(a);
    }
  };

  const onUp = (e: React.MouseEvent) => {
    const d = drag.current;
    drag.current = null;
    if (d?.kind === "boundary") {
      const h = hit(e);
      const row = p.rows[d.row];
      const nb = Math.max(d.minJ, Math.min(d.maxJ, Math.round(h.fx)));
      const a = d.orig.slice();
      if (nb > d.j) for (let k = d.j; k < nb; k++) a[row.positions[k]] = d.left;
      else for (let k = nb; k < d.j; k++) a[row.positions[k]] = d.right;
      p.onBoundaryPreview(null);
      if (nb !== d.j) p.onBoundaryCommit(a);
    }
  };

  useEffect(() => {
    const up = (e: MouseEvent) => {
      if (drag.current) onUp(e as unknown as React.MouseEvent);
    };
    window.addEventListener("mouseup", up);
    return () => window.removeEventListener("mouseup", up);
  });

  const onContext = (e: React.MouseEvent) => {
    e.preventDefault();
    const h = hit(e);
    const row = p.rows[h.row];
    if (row && h.j >= 0 && h.j < row.positions.length) {
      const pos = row.positions[h.j];
      if (!p.selected[pos]) p.onSelect([pos], "replace");
    }
    p.onContext(e.clientX, e.clientY);
  };

  const onWheel = (e: React.WheelEvent) => {
    if (e.ctrlKey || e.metaKey) {
      e.preventDefault();
      setZoom((z) => Math.max(0.4, Math.min(3, z * (e.deltaY < 0 ? 1.15 : 1 / 1.15))));
    }
  };

  return (
    <div className="seqbar">
      <div className="seqbar-tools">
        <span className="muted">Sequence</span>
        <button title="zoom out" onClick={() => setZoom((z) => Math.max(0.4, z / 1.25))}>−</button>
        <button title="zoom in" onClick={() => setZoom((z) => Math.min(3, z * 1.25))}>+</button>
        <span className="muted small">drag to select · drag black bars to move a boundary · 1–9 assign · 0 unassign · right-click for menu · ctrl+wheel zoom</span>
      </div>
      <div
        ref={wrap}
        className="seqbar-scroll"
        data-testid="seqbar"
        onScroll={(e) => setScroll((e.target as HTMLDivElement).scrollLeft)}
        onWheel={onWheel}
      >
        <div style={{ width: totalW, height: 1 }} />
        <canvas
          ref={canvas}
          style={{ position: "sticky", left: 0, top: 0, cursor, display: "block", marginTop: -1 }}
          onMouseDown={onDown}
          onMouseMove={onMove}
          onMouseLeave={() => p.onHover(null, 0, 0)}
          onContextMenu={onContext}
          data-cw={cw}
          data-labelw={LABEL_W}
          data-rowh={ROW_H}
          data-celly={TM_H + RULER_H}
        />
      </div>
    </div>
  );
}
