import { demoData, validateLocal } from "./demo";
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

export interface AppInfo {
  version: string;
  desktop: boolean;
  loaded: boolean;
  recent: string[][];
  out_dir: string;
  error?: string | null; // a file given on the command line could not be opened
}

const live = {
  app: () => req<AppInfo>("api/app"),
  open: (paths: string[]) => req<{ ok: boolean; notes: string[] }>("api/open", { paths }),
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

// Offline snapshot (kamd-regions snapshot): answer from the embedded data, validate in the browser.
const demo = demoData();
const OFFLINE = "This is an offline snapshot page; run kamd-regions locally to";

export const isDemo = demo !== null;

export const api: typeof live = demo
  ? {
      app: async () => ({ version: demo.session.system.version, desktop: false, loaded: true, recent: [], out_dir: "-" }),
      open: async () => {
        throw new Error(`${OFFLINE} open other files.`);
      },
      session: async () => demo.session,
      structure: async () => demo.structure,
      suggest: async (method, params) => {
        const r = demo.suggest[method];
        if (!r) throw new Error(`${OFFLINE} run '${method}'.`);
        if (params && Object.keys(params).length)
          throw new Error(`${OFFLINE} change method parameters; the snapshot has the default-parameter results only.`);
        return r;
      },
      validate: async (project) => validateLocal(demo.session.system, project),
      export: async () => {
        throw new Error(`${OFFLINE} export and verify the Amber/GROMACS/OpenMM inputs.`);
      },
      save: async () => ({ path: "(downloaded; the snapshot cannot write to disk)" }),
      importText: async () => {
        throw new Error(`${OFFLINE} import KAMD inputs.`);
      },
      loadProject: async (project) => ({
        project: project as ProjectDoc,
        notes: [],
        validation: validateLocal(demo.session.system, project as ProjectDoc),
      }),
    }
  : live;
