// Okabe-Ito (colour-blind safe) without black, then tab20.  Same order as the backend.
export const OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"];
export const TAB20 = [
  "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896",
  "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7",
  "#bcbd22", "#dbdb8d", "#17becf", "#9edae5",
];
export const PALETTE = [...OKABE_ITO, ...TAB20];

export function colorFor(i: number): string {
  return PALETTE[i % PALETTE.length];
}

export function unassignedColor(dark: boolean): string {
  return dark ? "#5c5f66" : "#d4d4d4";
}

export function luminance(hex: string): number {
  const h = hex.replace("#", "");
  const n = parseInt(h.length === 3 ? h.split("").map((c) => c + c).join("") : h, 16);
  const r = ((n >> 16) & 255) / 255;
  const g = ((n >> 8) & 255) / 255;
  const b = (n & 255) / 255;
  const lin = (c: number) => (c <= 0.03928 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);
  return 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b);
}

export function textOn(hex: string): string {
  return luminance(hex) > 0.35 ? "#111" : "#fff";
}
