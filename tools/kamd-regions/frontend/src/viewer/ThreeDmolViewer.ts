import * as $3Dmol from "3dmol";
import type { MolViewer, PickEvent, Representation, SlabSpec, ViewerData } from "./MolViewer";

type Atom = { index: number; x: number; y: number; z: number; ss?: string; ssbegin?: boolean; ssend?: boolean; resi: number; chain: string };

export class ThreeDmolViewer implements MolViewer {
  private v: any;
  private solute: any = null;
  private lipids: any = null;
  private ions: any = null;
  private resOf: number[] = [];
  private color: (r: number) => string = () => "#cccccc";
  private rep: Representation = "cartoon";
  private selected = new Set<number>();
  private highlight = "#000000";
  private shapes: any[] = [];
  private lipidsVisible = false;
  private ionsVisible = false;
  private pickCb: (e: PickEvent) => void = () => {};
  private hoverCb: (r: number | null) => void = () => {};
  private surfaceTimer: number | undefined;

  constructor(el: HTMLElement, background: string) {
    this.v = $3Dmol.createViewer(el, { backgroundColor: background, antialias: true } as any);
    this.v.setHoverDuration(80);
  }

  load(d: ViewerData): void {
    this.v.clear();
    this.resOf = d.soluteRes;
    this.solute = this.v.addModel(d.solutePdb, "pdb", { keepH: false, assignBonds: true });
    // use the backend's DSSP so the cartoon matches the sequence bar
    const atoms: Atom[] = this.solute.selectedAtoms({});
    let prevKey = "";
    let prevAtom: Atom | null = null;
    for (const a of atoms) {
      const r = this.resOf[a.index];
      const s = d.ss.get(r) ?? "";
      a.ss = s === "H" ? "h" : s === "E" ? "s" : "c";
      const key = `${a.chain}|${a.ss}`;
      a.ssbegin = false;
      a.ssend = false;
      if (prevAtom && key !== prevKey) prevAtom.ssend = true;
      if (key !== prevKey) a.ssbegin = true;
      prevKey = key;
      prevAtom = a;
    }
    if (prevAtom) prevAtom.ssend = true;
    this.lipids = d.lipidPdb ? this.v.addModel(d.lipidPdb, "pdb", { keepH: false }) : null;
    this.ions = d.ionPdb ? this.v.addModel(d.ionPdb, "pdb", { keepH: false }) : null;
    this.solute.setClickable({}, true, (atom: Atom, _v: unknown, ev: MouseEvent) => {
      const e = ev as MouseEvent | undefined;
      this.pickCb({ residue: this.resOf[atom.index], shift: !!e?.shiftKey, ctrl: !!(e?.ctrlKey || e?.metaKey) });
    });
    this.solute.setHoverable(
      {},
      true,
      (atom: Atom) => this.hoverCb(this.resOf[atom.index]),
      () => this.hoverCb(null),
    );
    this.restyle();
    this.v.zoomTo({ model: this.solute });
    this.v.render();
    // the container may have been laid out after the viewer was created
    window.setTimeout(() => this.resize(), 50);
  }

  orientNormalUp(axis: number): void {
    // side view: membrane normal vertical on screen
    if (axis === 3) this.v.rotate(-90, "x");
    else if (axis === 1) this.v.rotate(90, "z");
    this.v.zoomTo({ model: this.solute });
    this.v.render();
  }

  private colorfunc = (atom: Atom) => this.color(this.resOf[atom.index]);

  private restyle(): void {
    if (!this.solute) return;
    const m = this.solute;
    const sel = { model: m };
    if (this.rep === "cartoon") this.v.setStyle(sel, { cartoon: { colorfunc: this.colorfunc, arrows: true } });
    else if (this.rep === "sticks") this.v.setStyle(sel, { stick: { colorfunc: this.colorfunc, radius: 0.2 } });
    else this.v.setStyle(sel, { cartoon: { colorfunc: this.colorfunc, opacity: 0.6 } });
    if (this.selected.size) {
      const pred = (a: Atom) => this.selected.has(this.resOf[a.index]);
      this.v.addStyle({ model: m, predicate: pred }, { stick: { color: this.highlight, radius: 0.28 } });
    }
    if (this.lipids) {
      this.v.setStyle({ model: this.lipids }, this.lipidsVisible ? { line: { colorscheme: "grayCarbon", opacity: 0.7 } } : {});
    }
    if (this.ions) {
      this.v.setStyle({ model: this.ions }, this.ionsVisible ? { sphere: { colorscheme: "element", radius: 0.9 } } : {});
    }
    this.refreshSurface();
    this.v.render();
  }

  private refreshSurface(): void {
    window.clearTimeout(this.surfaceTimer);
    this.v.removeAllSurfaces();
    if (this.rep !== "surface" || !this.solute) return;
    // surfaces are slow; debounce so drags stay responsive
    this.surfaceTimer = window.setTimeout(() => {
      this.v.removeAllSurfaces();
      this.v.addSurface("VDW", { colorfunc: this.colorfunc, opacity: 0.92 } as any, { model: this.solute });
      this.v.render();
    }, 150);
  }

  setColors(color: (r: number) => string): void {
    this.color = color;
    this.restyle();
  }

  setRepresentation(rep: Representation): void {
    this.rep = rep;
    this.restyle();
  }

  setSelection(res: Set<number>, highlight: string): void {
    this.selected = res;
    this.highlight = highlight;
    // the surface does not depend on the selection; avoid rebuilding it
    const rep = this.rep;
    if (rep === "surface") {
      this.v.setStyle({ model: this.solute }, { cartoon: { colorfunc: this.colorfunc, opacity: 0.6 } });
      if (res.size) {
        const pred = (a: Atom) => this.selected.has(this.resOf[a.index]);
        this.v.addStyle({ model: this.solute, predicate: pred }, { stick: { color: highlight, radius: 0.28 } });
      }
      this.v.render();
    } else this.restyle();
  }

  setMembrane(slab: SlabSpec | null, visible: boolean, color: string): void {
    for (const s of this.shapes) this.v.removeShape(s);
    this.shapes = [];
    if (slab && visible && this.solute) {
      const atoms: Atom[] = this.solute.selectedAtoms({});
      const lo = [Infinity, Infinity, Infinity];
      const hi = [-Infinity, -Infinity, -Infinity];
      for (const a of atoms) {
        const p = [a.x, a.y, a.z];
        for (let k = 0; k < 3; k++) {
          lo[k] = Math.min(lo[k], p[k]);
          hi[k] = Math.max(hi[k], p[k]);
        }
      }
      const ax = slab.axis - 1;
      const pad = 12;
      for (const level of [slab.lo, slab.hi]) {
        const c = [0, 1, 2].map((k) => (k === ax ? level : (lo[k] + hi[k]) / 2));
        const dim = [0, 1, 2].map((k) => (k === ax ? 0.4 : hi[k] - lo[k] + 2 * pad));
        this.shapes.push(
          this.v.addBox({
            center: { x: c[0], y: c[1], z: c[2] },
            dimensions: { w: dim[0], h: dim[1], d: dim[2] },
            color,
            opacity: 0.28,
          }),
        );
      }
    }
    this.v.render();
  }

  setLipidsVisible(v: boolean): void {
    this.lipidsVisible = v;
    this.restyle();
  }

  setIonsVisible(v: boolean): void {
    this.ionsVisible = v;
    this.restyle();
  }

  setBackground(color: string): void {
    this.v.setBackgroundColor(color);
    this.v.render();
  }

  onPick(cb: (e: PickEvent) => void): void {
    this.pickCb = cb;
  }

  onHover(cb: (r: number | null) => void): void {
    this.hoverCb = cb;
  }

  zoomTo(res?: number[]): void {
    if (!this.solute) return;
    if (res && res.length) {
      const s = new Set(res);
      this.v.zoomTo({ model: this.solute, predicate: (a: Atom) => s.has(this.resOf[a.index]) }, 400);
    } else this.v.zoomTo({ model: this.solute }, 400);
  }

  resize(): void {
    this.v.resize();
    if (this.solute) this.v.zoomTo({ model: this.solute });
    this.v.render();
  }

  dispose(): void {
    this.v.clear();
  }
}
