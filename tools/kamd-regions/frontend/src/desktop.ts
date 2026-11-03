// Desktop mode (kamd-regions app): the page runs in a native window and talks to Python
// through pywebview's bridge (window.pywebview.api). In a normal browser the bridge is absent
// and the page falls back to browser downloads and path input.

export interface Bridge {
  open_dialog(multiple?: boolean): Promise<string[]>;
  folder_dialog(): Promise<string | null>;
  save_dialog(name: string): Promise<string | null>;
  save_text(name: string, text: string): Promise<string | null>;
  reveal(path: string): Promise<boolean>;
  quit(): Promise<void>;
}

export function bridge(): Bridge | null {
  return ((window as unknown as { pywebview?: { api?: Bridge } }).pywebview?.api as Bridge) ?? null;
}

/** Resolves with the bridge once pywebview has injected it, or null after `ms` in a browser. */
export function whenBridge(ms = 1500): Promise<Bridge | null> {
  const b = bridge();
  if (b && typeof b.open_dialog === "function") return Promise.resolve(b);
  return new Promise((resolve) => {
    const done = () => resolve(bridge());
    window.addEventListener("pywebviewready", done, { once: true });
    window.setTimeout(() => resolve(bridge()), ms);
  });
}

/** Save text through a native dialog (desktop) or a browser download. Returns the path or name. */
export async function saveText(name: string, text: string): Promise<string | null> {
  const b = bridge();
  if (b) return b.save_text(name, text);
  const url = URL.createObjectURL(new Blob([text], { type: "text/plain" }));
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
  return name;
}
