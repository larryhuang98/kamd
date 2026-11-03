import { useEffect, useRef, useState } from "react";
import type { MolViewer, PickEvent, Representation } from "../viewer/MolViewer";
import { ThreeDmolViewer } from "../viewer/ThreeDmolViewer";
import type { Membrane, Residues, StructureData } from "../types";

interface Props {
  structure: StructureData | null;
  res: Residues;
  colorOf: (residue: number) => string;
  selected: Set<number>;
  membrane: Membrane | { axis: number; lo: number; hi: number } | null;
  dark: boolean;
  onPick: (e: PickEvent) => void;
  onHover: (residue: number | null) => void;
  hoverLabel: string | null;
}

export function Viewer3D(p: Props) {
  const el = useRef<HTMLDivElement>(null);
  const viewer = useRef<MolViewer | null>(null);
  const [rep, setRep] = useState<Representation>("cartoon");
  const [lipids, setLipids] = useState(false);
  const [ions, setIons] = useState(false);
  const [slab, setSlab] = useState(true);
  const bg = p.dark ? "#16171a" : "#ffffff";

  useEffect(() => {
    if (!el.current) return;
    const v = new ThreeDmolViewer(el.current, bg);
    viewer.current = v;
    const ro = new ResizeObserver(() => v.resize());
    ro.observe(el.current);
    return () => {
      ro.disconnect();
      v.dispose();
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  useEffect(() => {
    const v = viewer.current;
    if (!v || !p.structure) return;
    const ss = new Map<number, string>();
    p.res.index.forEach((g, i) => ss.set(g, p.res.ss[i]));
    v.load({
      solutePdb: p.structure.solute_pdb,
      soluteRes: p.structure.solute_res,
      ss,
      lipidPdb: p.structure.lipid_pdb,
      ionPdb: p.structure.ion_pdb,
    });
    v.setColors(p.colorOf);
    if (p.membrane) v.orientNormalUp(p.membrane.axis);
    (window as any).__kamdViewerReady = true;
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [p.structure]);

  useEffect(() => {
    viewer.current?.onPick(p.onPick);
    viewer.current?.onHover(p.onHover);
  }, [p.onPick, p.onHover]);

  useEffect(() => viewer.current?.setColors(p.colorOf), [p.colorOf]);
  useEffect(() => viewer.current?.setSelection(p.selected, p.dark ? "#ffffff" : "#111111"), [p.selected, p.dark]);
  useEffect(() => viewer.current?.setRepresentation(rep), [rep]);
  useEffect(() => viewer.current?.setLipidsVisible(lipids), [lipids, p.structure]);
  useEffect(() => viewer.current?.setIonsVisible(ions), [ions, p.structure]);
  useEffect(() => viewer.current?.setBackground(bg), [bg]);
  useEffect(
    () => viewer.current?.setMembrane(p.membrane ? { axis: p.membrane.axis, lo: p.membrane.lo, hi: p.membrane.hi } : null, slab, p.dark ? "#9fb7d9" : "#7a93b8"),
    [p.membrane, slab, p.structure, p.dark],
  );

  return (
    <section className="viewer-wrap panel">
      <div className="viewer-tools">
        <select value={rep} onChange={(e) => setRep(e.target.value as Representation)} aria-label="representation">
          <option value="cartoon">cartoon</option>
          <option value="surface">surface</option>
          <option value="sticks">sticks</option>
        </select>
        {p.membrane && (
          <label className="chk">
            <input type="checkbox" checked={slab} onChange={(e) => setSlab(e.target.checked)} /> slab
          </label>
        )}
        {p.structure?.lipid_pdb && (
          <label className="chk">
            <input type="checkbox" checked={lipids} onChange={(e) => setLipids(e.target.checked)} /> lipids
          </label>
        )}
        {p.structure?.ion_pdb && (
          <label className="chk">
            <input type="checkbox" checked={ions} onChange={(e) => setIons(e.target.checked)} /> ions
          </label>
        )}
        <button className="ghost" onClick={() => viewer.current?.zoomTo(p.selected.size ? [...p.selected] : undefined)}>zoom</button>
        <span className="hover-label" data-testid="hover-label">{p.hoverLabel ?? ""}</span>
      </div>
      <div ref={el} className="viewer" data-testid="viewer" />
    </section>
  );
}
