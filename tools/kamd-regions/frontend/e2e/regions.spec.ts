import { expect, test, type Page } from "@playwright/test";
import { readFileSync } from "node:fs";

const AMBER = "http://127.0.0.1:8799/";
const VDAC = "http://127.0.0.1:8798/";
const SHOTS = process.env.SCREENSHOTS ? "../docs/img/" : null;

async function ready(page: Page) {
  await expect(page.getByTestId("seqbar")).toBeVisible();
  await page.waitForFunction(() => (window as any).__kamdViewerReady === true);
  await page.waitForTimeout(800); // first WebGL frame
}

// centre of the cell for position j in sequence row `row`
async function cell(page: Page, j: number, row = 0) {
  const c = page.getByTestId("seqbar").locator("canvas");
  const box = (await c.boundingBox())!;
  const d = await c.evaluate((el) => ({ ...(el as HTMLCanvasElement).dataset }));
  const cw = Number(d.cw), lw = Number(d.labelw), rh = Number(d.rowh), cy = Number(d.celly);
  return { x: box.x + lw + j * cw + cw / 2, y: box.y + row * rh + cy + 9, cw };
}

test("drag-select, new region, re-suggest, undo, export", async ({ page }) => {
  page.on("dialog", (d) => d.accept());
  await page.goto(AMBER);
  await ready(page);
  // the imported &kamd groups
  await expect(page.locator(".region-list .region")).toHaveCount(4);
  await expect(page.getByTestId("region-1")).toContainText("A:2-5");
  await expect(page.getByTestId("region-3")).toContainText("A:23-30");
  await expect(page.getByTestId("validation")).not.toContainText("error");

  // drag over residues 10-15 (positions 9..14) and make a new region with the keyboard shortcut
  const a = await cell(page, 9);
  const b = await cell(page, 14);
  await page.mouse.move(a.x, a.y);
  await page.mouse.down();
  await page.mouse.move((a.x + b.x) / 2, a.y, { steps: 4 });
  await page.mouse.move(b.x, b.y, { steps: 4 });
  await page.mouse.up();
  await page.keyboard.press("5"); // 5 = one more than the 4 regions -> new region
  await expect(page.getByTestId("region-5")).toContainText("A:10-15");
  await expect(page.getByTestId("region-5")).toContainText("atoms");
  if (SHOTS) await page.screenshot({ path: SHOTS + "amber_test.png" });

  // reassign residue 16 to region 5 via the context menu
  const c16 = await cell(page, 15);
  await page.mouse.click(c16.x, c16.y);
  await page.mouse.click(c16.x, c16.y, { button: "right" });
  await page.getByTestId("context-menu").getByText("5. region").click();
  await expect(page.getByTestId("region-5")).toContainText("A:10-16");

  // move the boundary between region 2 (A:18-20) and region 3 (A:23-30)? first between 5 and unassigned
  const edge = await cell(page, 15);
  await page.mouse.move(edge.x + edge.cw / 2, edge.y);
  await page.mouse.down();
  await page.mouse.move(edge.x + edge.cw / 2 + edge.cw, edge.y, { steps: 5 });
  await page.mouse.up();
  await expect(page.getByTestId("region-5")).toContainText("A:10-17");

  // undo / redo
  await page.keyboard.press("Control+z");
  await expect(page.getByTestId("region-5")).toContainText("A:10-16");
  await page.keyboard.press("Control+Shift+z");
  await expect(page.getByTestId("region-5")).toContainText("A:10-17");

  // export and check the files
  await page.getByTestId("export").click();
  await expect(page.getByTestId("export-result")).toContainText("verified");
  const pre = page.getByTestId("amber-preview");
  await expect(pre).toContainText('refer_mask(1) = ":2-5"');
  await expect(pre).toContainText('refer_mask(5) = ":10-17"');
  await expect(pre).toContainText('solvent_mask    = ":40-9692"');
  const dl = page.waitForEvent("download");
  await page.getByTestId("file-index.ndx").click();
  const ndx = readFileSync(await (await dl).path(), "utf8");
  // residues 10-17 = atoms 166..284 (1-based; residue 18 starts at atom 285)
  const g5 = ndx.split("[ region_5 ]")[1].split("[")[0].trim().split(/\s+/).map(Number);
  expect(g5[0]).toBe(166);
  expect(g5[g5.length - 1]).toBe(284);
  expect(g5.length).toBe(284 - 166 + 1);

  // re-suggest with GNM replaces the edited regions (after the confirm dialog)
  await page.getByTestId("method-select").selectOption("gnm");
  await page.getByTestId("resuggest").click();
  await expect(page.locator(".region-list .region")).toHaveCount(2);
  await expect(page.getByTestId("region-1")).toContainText("A:9-19");
  await page.keyboard.press("Control+z");
  await expect(page.locator(".region-list .region")).toHaveCount(5);
});

test("3D click selects in the sequence bar", async ({ page }) => {
  await page.goto(AMBER);
  await ready(page);
  // clicking a sequence cell highlights it in 3D; hover text comes from the 3D view
  const c = await cell(page, 2);
  await page.mouse.click(c.x, c.y);
  await page.keyboard.press("0");
  await expect(page.getByTestId("region-1")).toContainText("A:2, A:4-5");
});

test("membrane protein: tm groups and slab", async ({ page }) => {
  await page.goto(VDAC);
  await ready(page);
  await expect(page.locator("header")).toContainText("membrane z -11.5");
  await expect(page.getByTestId("region-1").locator("input.name")).toHaveValue("TM 1");
  await expect(page.getByTestId("validation")).toContainText("need a membrane bath");
  await expect(page.getByTestId("export")).toBeDisabled();
  if (SHOTS) {
    await page.screenshot({ path: SHOTS + "vdac_tm.png" });
    await page.locator("header button[title='toggle theme']").click();
    await page.waitForTimeout(600);
    await page.screenshot({ path: SHOTS + "vdac_tm_dark.png" });
    await page.locator("header button[title='toggle theme']").click();
  }
});
