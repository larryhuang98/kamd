import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// The production build goes into the Python package so `pip install` ships the GUI.
export default defineConfig({
  plugins: [react()],
  base: "./",
  build: {
    outDir: "../kamd_regions/static",
    emptyOutDir: true,
    chunkSizeWarningLimit: 4000,
  },
  server: { proxy: { "/api": "http://127.0.0.1:8765" } },
});
