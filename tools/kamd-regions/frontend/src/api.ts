import type { ExportResult, ProjectDoc, StructureData, SystemInfo, Validation } from "./types";

async function req<T>(path: string, body?: unknown): Promise<T> {
  const r = await fetch(path, {
    method: body === undefined ? "GET" : "POST",
    headers: body === undefined ? undefined : { "Content-Type": "application/json" },
    body: body === undefined ? undefined : JSON.stringify(body),
  });
  if (!r.ok) {
    let msg = `${r.status} ${r.statusText}`;
    try {
      const j = await r.json();
      if (j.detail) msg = typeof j.detail === "string" ? j.detail : JSON.stringify(j.detail);
    } catch {
      /* not JSON */
    }
    throw new Error(msg);
  }
  return (await r.json()) as T;
}

export const api = {
  session: () =>
    req<{ system: SystemInfo; project: ProjectDoc; validation: Validation; out_dir: string }>("api/session"),
  structure: () => req<StructureData>("api/structure"),
  suggest: (method: string, params: Record<string, unknown>, settings: unknown) =>
    req<{ project: ProjectDoc; notes: string[]; validation: Validation }>("api/suggest", { method, params, settings }),
  validate: (project: ProjectDoc) => req<Validation>("api/validate", { project }),
  export: (project: ProjectDoc, engines: string[], write = false, out_dir?: string) =>
    req<ExportResult>("api/export", { project, engines, write, out_dir }),
  save: (project: ProjectDoc, path?: string) => req<{ path: string }>("api/project/save", { project, path }),
  importText: (kind: string, text: string, ndx_text?: string) =>
    req<{ project: ProjectDoc; notes: string[]; validation: Validation }>("api/import", { kind, text, ndx_text }),
  loadProject: (project: unknown) =>
    req<{ project: ProjectDoc; notes: string[]; validation: Validation }>("api/project/load", { project }),
};
