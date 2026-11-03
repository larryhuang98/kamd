import { defineConfig } from "@playwright/test";

// The e2e test starts the real backend on the Amber test system (and VDAC for screenshots).
// KAMD_REGIONS_CMD can point to a specific executable (default: kamd-regions on PATH).
const cmd = process.env.KAMD_REGIONS_CMD ?? "kamd-regions";
const repo = "../../..";
const out = "test-results/server-out";

export default defineConfig({
  testDir: "e2e",
  timeout: 90_000,
  use: {
    viewport: { width: 1400, height: 900 },
    deviceScaleFactor: process.env.SCREENSHOTS ? 2 : 1,
    launchOptions: {
      // PW_CHROMIUM: use a preinstalled Chromium instead of the one matching this Playwright version
      executablePath: process.env.PW_CHROMIUM || undefined,
      args: ["--use-gl=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"],
    },
    acceptDownloads: true,
  },
  webServer: [
    {
      command: `${cmd} serve ${repo}/amber/test/prmtop ${repo}/amber/test/restrt ${repo}/amber/test/md_test.in --no-browser --port 8799 --out-dir ${out}`,
      url: "http://127.0.0.1:8799/api/session",
      timeout: 120_000,
      reuseExistingServer: false,
    },
    {
      command: `${cmd} serve ../tests/data/3emn.pdb --no-browser --port 8798 --out-dir ${out}-vdac`,
      url: "http://127.0.0.1:8798/api/session",
      timeout: 120_000,
      reuseExistingServer: false,
    },
    {
      command: `${cmd} serve --no-browser --port 8797 --out-dir ${out}-empty`,
      url: "http://127.0.0.1:8797/api/app",
      timeout: 120_000,
      reuseExistingServer: false,
      env: { KAMD_REGIONS_CONFIG: "test-results/cfg-empty" },
    },
  ],
});
