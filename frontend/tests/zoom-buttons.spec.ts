import { test } from './testWorkspace';
import { expect, type Page, type TestInfo } from '@playwright/test';
import { writeFile } from 'node:fs/promises';

type FrameSample = { ms: number; colored: number; width: number; height: number; image: string };
type ZoomInput = { at: number } & ({ label: string } | { pinchDeltaY: number });

/** Sample the real WebGL image on every animation frame. Final-state-only
 * screenshots missed the camera crossing the molecule during an animated zoom. */
async function sampleZoom(page: Page, clicks: ZoomInput[]) {
  return page.evaluate(async (clicks): Promise<FrameSample[]> => {
    const source = document.querySelector<HTMLCanvasElement>('.molecular-viewer__canvas canvas')!;
    const canvas = document.createElement('canvas');
    canvas.width = 480;
    canvas.height = Math.round((source.height / source.width) * canvas.width);
    const ctx = canvas.getContext('2d', { willReadFrequently: true })!;
    const start = performance.now();
    let next = 0;
    const samples: FrameSample[] = [];
    return new Promise((resolve) => {
      const sample = (now: number) => {
        const ms = now - start;
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.drawImage(source, 0, 0, canvas.width, canvas.height);
        const pixels = ctx.getImageData(0, 0, canvas.width, canvas.height).data;
        let colored = 0,
          minX = canvas.width,
          maxX = 0,
          minY = canvas.height,
          maxY = 0;
        for (let i = 0; i < pixels.length; i += 4) {
          const max = Math.max(pixels[i], pixels[i + 1], pixels[i + 2]);
          const min = Math.min(pixels[i], pixels[i + 1], pixels[i + 2]);
          if (max <= 75 || max - min <= 35) continue;
          colored++;
          const x = (i / 4) % canvas.width;
          const y = Math.floor(i / 4 / canvas.width);
          minX = Math.min(minX, x);
          maxX = Math.max(maxX, x);
          minY = Math.min(minY, y);
          maxY = Math.max(maxY, y);
        }
        samples.push({
          ms,
          colored,
          width: maxX - minX,
          height: maxY - minY,
          image: canvas.toDataURL('image/png'),
        });
        while (next < clicks.length && ms >= clicks[next].at) {
          const input = clicks[next];
          if ('label' in input) {
            document
              .querySelector<HTMLButtonElement>(`button[aria-label="${input.label}"]`)!
              .click();
          } else {
            source.dispatchEvent(
              new WheelEvent('wheel', {
                bubbles: true,
                cancelable: true,
                ctrlKey: true,
                deltaY: input.pinchDeltaY,
              }),
            );
          }
          next++;
        }
        if (ms >= clicks.at(-1)!.at + 550) resolve(samples);
        else requestAnimationFrame(sample);
      };
      requestAnimationFrame(sample);
    });
  }, clicks);
}

async function evidence(info: TestInfo, samples: FrameSample[], name: string) {
  const data = JSON.stringify(
    samples.map(({ image: _, ...sample }) => sample),
    null,
    2,
  );
  await writeFile(info.outputPath(`${name}-frames.json`), data);
  await info.attach(`${name}-frames.json`, {
    body: data,
    contentType: 'application/json',
  });
  const maximum = samples.reduce((best, sample) => (sample.colored > best.colored ? sample : best));
  for (const [label, sample] of [
    ['start', samples[0]],
    ['largest', maximum],
    ['end', samples.at(-1)!],
  ] as const) {
    await writeFile(
      info.outputPath(`${name}-${label}.png`),
      Buffer.from(sample.image.split(',')[1], 'base64'),
    );
    await info.attach(`${name}-${label}.png`, {
      body: Buffer.from(sample.image.split(',')[1], 'base64'),
      contentType: 'image/png',
    });
  }
}

test('zoom buttons change molecular scale smoothly throughout each animation', async ({
  page,
}, info) => {
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await page.waitForTimeout(700);
  for (const direction of ['in', 'out'] as const) {
    const samples = await sampleZoom(page, [{ at: 50, label: `Zoom ${direction}` }]);
    await evidence(info, samples, direction);
    expect(samples.length).toBeGreaterThan(12);
    const start = samples[0].colored,
      end = samples.at(-1)!.colored;
    expect(start).toBeGreaterThan(500);
    expect(direction === 'in' ? end / start : start / end).toBeGreaterThan(1.15);
    // Every intermediate frame must stay within the start/end molecular size.
    // A small rasterization tolerance permits antialiasing at changing edges.
    for (const sample of samples) {
      expect(
        sample.colored,
        `${direction}: transient enlargement at ${sample.ms.toFixed(1)} ms`,
      ).toBeLessThanOrEqual(Math.max(start, end) * 1.08);
      expect(
        sample.colored,
        `${direction}: transient disappearance at ${sample.ms.toFixed(1)} ms`,
      ).toBeGreaterThanOrEqual(Math.min(start, end) * 0.92);
    }
  }
  expect(errors).toEqual([]);
});

test('pinching interrupts an unfinished button zoom without a later rebound', async ({
  page,
}, info) => {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.waitForTimeout(700);
  const samples = await sampleZoom(page, [
    { at: 50, label: 'Zoom in' },
    { at: 120, pinchDeltaY: 30 },
  ]);
  await evidence(info, samples, 'pinch-interruption');
  const settled = samples.filter((sample) => sample.ms >= 180);
  const final = settled.at(-1)!;
  expect(final.width / samples[0].width).toBeGreaterThan(0.72);
  expect(final.width / samples[0].width).toBeLessThan(0.87);
  for (const sample of settled) {
    expect(Math.abs(sample.width - final.width)).toBeLessThanOrEqual(1);
    expect(Math.abs(sample.height - final.height)).toBeLessThanOrEqual(1);
    // NGL refines edge antialiasing after motion stops; the extent must remain
    // fixed even when the count of edge pixels changes slightly.
    expect(Math.abs(sample.colored - final.colored)).toBeLessThan(final.colored * 0.08);
  }
});

test('rapid repeated zoom clicks accumulate without overlapping camera jumps', async ({
  page,
}, info) => {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.waitForTimeout(700);
  const samples = await sampleZoom(
    page,
    [50, 95, 140].map((at) => ({ at, label: 'Zoom out' })),
  );
  await evidence(info, samples, 'repeated-out');
  const start = samples[0].colored,
    end = samples.at(-1)!.colored;
  // Three 0.8 scale factors should remain three steps even if clicked before
  // the preceding animation finishes. Pixel area scales by 0.8 ** 6.
  expect(end / start).toBeGreaterThan(0.18);
  expect(end / start).toBeLessThan(0.38);
  for (let index = 1; index < samples.length; index++) {
    expect(samples[index].colored).toBeLessThanOrEqual(start * 1.05);
    expect(samples[index].colored - samples[index - 1].colored).toBeLessThanOrEqual(start * 0.02);
  }
});
