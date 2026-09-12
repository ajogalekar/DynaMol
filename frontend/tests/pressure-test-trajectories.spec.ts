import { expect, type Page } from '@playwright/test';
import { readFileSync } from 'node:fs';
import { stat, writeFile } from 'node:fs/promises';
import { test } from './testWorkspace';

type Case = { name: string; pdbId: string; jobId: string; datasetId: string; baseURL: string };
type Atom = {
  index: number;
  name: string;
  residue: string;
  resid: number;
  element: string;
  category: string;
};
type Series = {
  label: string;
  kind: string;
  unit: string;
  output_atoms: number[];
  values: number[];
  times_ps: number[];
};
const manifestPath = process.env.DYNAMOL_PRESSURE_TEST_MANIFEST;
const rawCases = manifestPath ? JSON.parse(readFileSync(manifestPath, 'utf8')).cases : [];
if (!Array.isArray(rawCases)) throw new Error('Pressure-test manifest must contain cases[]');
const cases: Case[] = rawCases.map((raw: Record<string, string>) => {
  const entry = {
    name: raw.name || raw.pdbId || raw.jobId,
    pdbId: raw.pdbId || raw.pdb_id,
    jobId: raw.jobId || raw.job_id,
    datasetId: raw.datasetId || raw.dataset_id,
    baseURL: raw.baseURL || raw.base_url,
  };
  const url = new URL(entry.baseURL);
  if (
    !['localhost', '127.0.0.1', '[::1]'].includes(url.hostname) ||
    !url.port ||
    ['8765', '5173', '4173'].includes(url.port)
  )
    throw new Error('Pressure-test cases require a dedicated isolated loopback server');
  if (![entry.jobId, entry.datasetId].every((value) => /^[a-zA-Z0-9_-]+$/.test(value)))
    throw new Error('Each pressure-test case requires a completed jobId and datasetId');
  return entry;
});
const escape = (value: string) => value.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
async function pick(page: Page, atom: Atom) {
  await page
    .getByRole('textbox', { name: 'Find an atom', exact: true })
    .fill(`${atom.residue}${atom.resid} · ${atom.name}`);
  await page
    .locator('.atom-search-results')
    .getByRole('button', {
      name: new RegExp(`#${atom.index + 1} · ${escape(atom.element)}$`),
    })
    .click();
}
// Explicitly opt in: ordinary suites must not reset any workspace for this audit.
if (!cases.length)
  test.describe.skip('pressure-test manifest required', () => {
    test('set DYNAMOL_PRESSURE_TEST_MANIFEST to completed outputs', async () => {});
  });

for (const entry of cases)
  test.describe(`${entry.pdbId}: ${entry.name}`, () => {
    test.use({ baseURL: entry.baseURL });
    test('completed native trajectory opens, plots, plays, scrubs and downloads', async ({
      page,
      request,
    }, info) => {
      test.setTimeout(240_000);
      const browserErrors: string[] = [],
        consoleErrors: string[] = [],
        forbiddenRequests: string[] = [];
      const checks: string[] = [],
        posthoc: unknown[] = [];
      const evidence: Record<string, unknown> = {
        case: entry,
        startedAt: new Date().toISOString(),
        checks,
        posthoc,
        browserErrors,
        consoleErrors,
        forbiddenRequests,
      };
      page.on('pageerror', (error) => browserErrors.push(error.message));
      page.on('console', (message) => {
        if (message.type() === 'error') consoleErrors.push(message.text());
      });
      page.on('request', (req) => {
        if (
          req.method() === 'POST' &&
          /^\/api\/jobs(?:$|\/[^/]+\/(resume|cancel)$)/.test(new URL(req.url()).pathname)
        )
          forbiddenRequests.push(`${req.method()} ${req.url()}`);
      });
      try {
        const responses = await Promise.all([
          request.get(`/api/jobs/${entry.jobId}`),
          request.get(`/api/datasets/${entry.datasetId}`),
          request.get(`/api/jobs/${entry.jobId}/measurements`),
        ]);
        for (const response of responses) expect(response.ok(), await response.text()).toBeTruthy();
        const [job, dataset, snapshot] = await Promise.all(responses.map((r) => r.json()));
        expect(job.status).toBe('completed');
        expect(job.dataset_id).toBe(entry.datasetId);
        expect(dataset.n_frames).toBe(101);
        expect(dataset.time_unit).toBe('ps');
        expect(dataset.times_ps).toHaveLength(101);
        expect(dataset.times_ps[0]).toBeCloseTo(0, 4);
        expect(dataset.times_ps[100]).toBeCloseTo(100, 3);
        expect(snapshot.output_dataset_id).toBe(entry.datasetId);
        const series = snapshot.measurements as Series[];
        expect(series.map((m) => m.kind).sort()).toEqual(['angle', 'dihedral', 'distance']);
        for (const m of series) {
          expect(m.values).toHaveLength(101);
          expect(m.times_ps).toHaveLength(101);
          expect(m.values.every(Number.isFinite)).toBeTruthy();
          expect(m.output_atoms.length).toBe(
            ({ distance: 2, angle: 3, dihedral: 4 } as Record<string, number>)[m.kind],
          );
        }
        const waterAtoms = (dataset.atoms as Atom[]).filter((a) => a.category === 'water').length;
        expect(waterAtoms).toBeGreaterThan(0);
        evidence.nativeOutput = {
          engine: job.engine,
          jobName: job.name,
          atoms: dataset.n_atoms,
          frames: dataset.n_frames,
          firstTimePs: dataset.times_ps[0],
          lastTimePs: dataset.times_ps[100],
          waterAtoms,
        };
        checks.push(
          'completed output: 101 finite saved samples, 0–100 ps and three mapped measurements',
        );
        await page.goto('/');
        await expect(
          page.getByRole('button', { name: 'Save snapshot', exact: true }),
        ).toBeEnabled();
        await expect(page.locator('.canvas-loading')).toHaveCount(0);
        await page.getByRole('button', { name: 'New simulation', exact: true }).click();
        const studio = page.getByRole('dialog', { name: 'Set molecules in motion.' });
        const card = studio.locator('.job-card').filter({
          has: page.locator('.job-top strong', {
            hasText: new RegExp(`^${escape(job.name)}$`),
          }),
        });
        await expect(card).toHaveCount(1);
        await expect(card.locator('.job-status')).toContainText('completed');
        const downloadEvent = page.waitForEvent('download');
        await card.getByRole('link', { name: 'Files', exact: true }).click();
        const download = await downloadEvent;
        expect(await download.failure()).toBeNull();
        const downloadPath = await download.path();
        expect(downloadPath).toBeTruthy();
        const bytes = (await stat(downloadPath!)).size;
        expect(bytes).toBeGreaterThan(100);
        expect(download.suggestedFilename()).toContain(entry.jobId);
        evidence.download = { filename: download.suggestedFilename(), bytes };
        checks.push('real native job Files archive downloaded');
        await card.getByRole('button', { name: 'Open trajectory', exact: true }).click();
        await expect(studio).toHaveCount(0);
        await expect(page.getByRole('button', { name: 'Explore', exact: true })).toHaveClass(
          'active',
        );
        await expect(page.locator('.structure-card h2')).toHaveText(dataset.name);
        await expect(page.locator('.canvas-loading')).toHaveCount(0, { timeout: 90_000 });
        await expect(page.locator('.measurement-tab')).toHaveCount(3);
        await expect(page.locator('.timeline-meta')).toContainText('101 saved frames');
        await expect(page.locator('.plot-live-status')).toHaveCount(0);
        for (const m of series) {
          const tab = page
            .locator('.measurement-tab')
            .filter({ has: page.getByRole('button', { name: m.label, exact: true }) });
          await expect(tab).toBeVisible();
          await expect(
            tab.getByRole('button', { name: `Hide ${m.label} in view`, exact: true }),
          ).toHaveAttribute('aria-pressed', 'true');
          await tab.getByRole('button', { name: m.label, exact: true }).click();
          await expect(
            page.getByRole('img', {
              name: `${m.label}, plotted across 101 saved frames`,
              exact: true,
            }),
          ).toBeVisible();
        }
        checks.push('Open trajectory enters Explore and restores all three visible saved plots');
        const slider = page.getByRole('slider', { name: 'Trajectory frame', exact: true });
        await expect(slider).toHaveAttribute('max', '100');
        await slider.fill('0');
        await expect(page.locator('.time-display b')).toHaveText('0.00');
        await slider.fill('100');
        await expect(page.locator('.time-display b')).toHaveText('100.00');
        await page.getByRole('button', { name: 'First frame', exact: true }).click();
        await expect(slider).toHaveValue('0');
        await page.getByRole('button', { name: 'Last frame', exact: true }).click();
        await expect(slider).toHaveValue('100');
        await page.getByRole('button', { name: 'First frame', exact: true }).click();
        await page.getByRole('button', { name: 'Play trajectory', exact: true }).click();
        await expect.poll(async () => Number(await slider.inputValue())).toBeGreaterThan(2);
        await page.getByRole('button', { name: 'Pause trajectory', exact: true }).click();
        const paused = await slider.inputValue();
        await page.waitForTimeout(350); // Several animation ticks verify an actual pause.
        await expect(slider).toHaveValue(paused);
        evidence.playback = { pausedFrame: Number(paused), first: 0, last: 100 };
        checks.push('movie advances and pauses; slider and endpoint buttons reach both ends');
        await slider.fill('0');
        await expect(page.locator('.canvas-loading')).toHaveCount(0);
        await expect(page.getByRole('button', { name: 'Hide water', exact: true })).toHaveAttribute(
          'aria-pressed',
          'true',
        );
        await page.getByRole('button', { name: 'Hide water', exact: true }).click();
        await expect(page.getByRole('button', { name: 'Show water', exact: true })).toHaveAttribute(
          'aria-pressed',
          'false',
        );
        await page.getByRole('button', { name: 'Show water', exact: true }).click();
        await expect(page.getByRole('button', { name: 'Hide water', exact: true })).toHaveAttribute(
          'aria-pressed',
          'true',
        );
        await page.getByRole('button', { name: 'Hide water', exact: true }).click();
        await page.waitForTimeout(600);
        const canvas = page.locator('.molecular-viewer canvas').first();
        const beforeZoom = await canvas.screenshot();
        await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
        await page.waitForTimeout(500);
        const zoomedIn = await canvas.screenshot();
        expect(zoomedIn.equals(beforeZoom)).toBeFalsy();
        await page.getByRole('button', { name: 'Zoom out', exact: true }).click();
        await page.waitForTimeout(500);
        expect((await canvas.screenshot()).equals(zoomedIn)).toBeFalsy();
        checks.push('water toggles both ways and both zoom buttons change the molecular canvas');
        for (const original of series) {
          await page.getByRole('button', { name: 'Measure', exact: true }).click();
          await page
            .locator('.measure-types')
            .getByRole('button', { name: new RegExp(`${original.kind}$`, 'i') })
            .click();
          const resultEvent = page.waitForResponse(
            (r) =>
              r.url().endsWith(`/api/datasets/${entry.datasetId}/measurements`) &&
              r.request().method() === 'POST' &&
              r.request().postDataJSON()?.kind === original.kind,
          );
          for (const index of original.output_atoms) {
            const atom = (dataset.atoms as Atom[]).find((a) => a.index === index);
            expect(atom, `mapped output atom ${index}`).toBeTruthy();
            await pick(page, atom!);
          }
          const result = await resultEvent;
          expect(result.ok(), await result.text()).toBeTruthy();
          const measured = await result.json();
          expect(measured.values).toHaveLength(101);
          expect(measured.unit).toBe(original.unit);
          const maxDifference = Math.max(
            ...measured.values.map((value: number, index: number) => {
              expect(Number.isFinite(value)).toBeTruthy();
              const delta = Math.abs(value - original.values[index]);
              return original.kind === 'dihedral' ? Math.abs(((delta + 180) % 360) - 180) : delta;
            }),
          );
          expect(maxDifference).toBeLessThan(0.001);
          await expect(page.locator('.measurement-draft')).toHaveCount(0);
          await expect(page.locator('.measurement-tab')).toHaveCount(4);
          await expect(page.locator('.plot-readout strong')).toHaveText(
            `${measured.values[0].toFixed(2)}${measured.unit}`,
          );
          await expect(page.locator('.chart-wrap svg')).toHaveAttribute(
            'aria-label',
            /plotted across 101 saved frames$/,
          );
          posthoc.push({
            kind: original.kind,
            atoms: original.output_atoms,
            samples: 101,
            maxDifference,
            unit: original.unit,
          });
          await page
            .locator('.measurement-tab.active')
            .getByRole('button', { name: /^Remove / })
            .click();
          await expect(page.locator('.measurement-tab')).toHaveCount(3);
        }
        checks.push('three real UI posthoc plots agree with the saved native tracking samples');
        await expect
          .poll(async () => (await (await request.get('/api/workspace')).json()).state.dataset_id)
          .toBe(entry.datasetId);
        expect(browserErrors).toEqual([]);
        expect(consoleErrors).toEqual([]);
        expect(forbiddenRequests).toEqual([]);
        evidence.passed = true;
      } finally {
        evidence.finishedAt = new Date().toISOString();
        await page
          .screenshot({ path: info.outputPath('trajectory-audit.png'), fullPage: true })
          .catch(() => {});
        await writeFile(
          info.outputPath('trajectory-audit.json'),
          JSON.stringify(evidence, null, 2),
        );
        await info.attach('trajectory-audit', {
          body: JSON.stringify(evidence, null, 2),
          contentType: 'application/json',
        });
      }
    });
  });
