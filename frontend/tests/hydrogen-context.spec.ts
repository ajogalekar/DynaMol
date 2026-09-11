import { test } from './testWorkspace';
import { expect, type Page, type TestInfo } from '@playwright/test';

type ScenePixels = { neutral: number; white: number; ribbon: number };

async function scenePixels(page: Page, info: TestInfo, name: string): Promise<ScenePixels> {
  // NGL builds representations asynchronously after a visibility change.
  await page.waitForTimeout(400);
  const image = await page.locator('.molecular-viewer__canvas canvas').screenshot();
  await info.attach(name, { body: image, contentType: 'image/png' });
  return page.evaluate(async (base64) => {
    const bitmap = await createImageBitmap(
      await (await fetch(`data:image/png;base64,${base64}`)).blob(),
    );
    const canvas = document.createElement('canvas');
    canvas.width = bitmap.width;
    canvas.height = bitmap.height;
    const context = canvas.getContext('2d')!;
    context.drawImage(bitmap, 0, 0);
    const pixels = context.getImageData(0, 0, canvas.width, canvas.height).data;
    const counts = { neutral: 0, white: 0, ribbon: 0 };
    for (let i = 0; i < pixels.length; i += 4) {
      const [r, g, b] = [pixels[i], pixels[i + 1], pixels[i + 2]];
      const min = Math.min(r, g, b),
        max = Math.max(r, g, b);
      if (max - min < 12 && min > 65 && max < 195) counts.neutral++;
      if (max - min < 12 && min >= 195) counts.white++;
      if (b > 95 && g > r + 20 && b > r + 25) counts.ribbon++;
    }
    return counts;
  }, image.toString('base64'));
}

test('polar hydrogen display keeps the carbon skeleton and ribbon through visibility changes', async ({
  page,
  request,
}, info) => {
  const fixture = await request.get('/api/datasets/demo');
  expect(fixture.ok()).toBeTruthy();
  const dataset = await fixture.json();
  expect(
    dataset.atoms.filter((atom: { element: string }) => atom.element === 'H').length,
  ).toBeGreaterThan(300);
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.getByRole('button', { name: 'Auto rotate', exact: true })).toHaveAttribute(
    'aria-pressed',
    'false',
  );
  await page.getByRole('button', { name: 'First frame', exact: true }).click();
  await page.getByRole('button', { name: 'Fit molecule', exact: true }).click();
  await page.getByRole('button', { name: 'Hidden', exact: true }).click();
  const hidden = await scenePixels(page, info, 'hydrogens-hidden');
  await page.getByRole('button', { name: 'Polar only', exact: true }).click();
  const polar = await scenePixels(page, info, 'polar-hydrogens-with-context');
  await page.getByRole('button', { name: 'All', exact: true }).click();
  const all = await scenePixels(page, info, 'all-hydrogens-with-context');
  await page.getByRole('button', { name: 'Polar only', exact: true }).click();
  const polarAgain = await scenePixels(page, info, 'polar-hydrogens-restored');
  await page.getByRole('button', { name: 'Hide protein & nucleic acids', exact: true }).click();
  const proteinHidden = await scenePixels(page, info, 'protein-hidden');
  await page.getByRole('button', { name: 'Show protein & nucleic acids', exact: true }).click();
  const proteinShown = await scenePixels(page, info, 'protein-restored');
  await page.getByRole('button', { name: 'Hidden', exact: true }).click();
  const hiddenAgain = await scenePixels(page, info, 'hydrogens-hidden-restored');
  console.log(
    JSON.stringify({ hidden, polar, all, polarAgain, proteinHidden, proteinShown, hiddenAgain }),
  );
  // The carbon skeleton is neutral grey, whereas ribbons are teal/purple and
  // hydrogen spheres are predominantly white. The previous disconnected H/parent
  // overlay retained only 10% of All's neutral pixels in Polar mode. With the
  // complete heavy-atom graph it retains 51% (All adds shaded hydrogen pixels
  // too). A generous 35% threshold separates those outcomes without requiring
  // an exact image match across graphics drivers.
  expect(
    polar.neutral - hidden.neutral,
    'Polar mode must render a substantial connected carbon skeleton.',
  ).toBeGreaterThan(1000);
  expect((polar.neutral - hidden.neutral) / (all.neutral - hidden.neutral)).toBeGreaterThan(0.35);
  expect(
    polar.ribbon,
    'The protein ribbon must remain visible behind atomic context.',
  ).toBeGreaterThan(hidden.ribbon * 0.7);
  expect(polar.white).toBeGreaterThan(hidden.white + 20);
  expect(all.white, 'All hydrogens must add visible non-polar hydrogens.').toBeGreaterThan(
    polar.white * 1.5,
  );
  expect(proteinHidden.neutral).toBeLessThan(polar.neutral * 0.1);
  expect(proteinHidden.ribbon).toBeLessThan(polar.ribbon * 0.1);
  for (const restored of [polarAgain, proteinShown]) {
    for (const key of ['neutral', 'white', 'ribbon'] as const) {
      expect(Math.abs(restored[key] - polar[key])).toBeLessThan(Math.max(5, polar[key] * 0.02));
    }
  }
  for (const key of ['neutral', 'white', 'ribbon'] as const) {
    expect(Math.abs(hiddenAgain[key] - hidden[key])).toBeLessThan(Math.max(5, hidden[key] * 0.02));
  }
  await expect(page.getByRole('alert')).toHaveCount(0);
});
