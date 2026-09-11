import { expect, test, type Locator, type Page } from '@playwright/test';

async function expectOnscreen(page: Page, locator: Locator) {
  await expect(locator).toBeVisible();
  const bounds = await locator.boundingBox();
  const viewport = await page.evaluate(() => ({ width: innerWidth, height: innerHeight }));
  expect(bounds).not.toBeNull();
  expect(bounds!.x).toBeGreaterThanOrEqual(0);
  expect(bounds!.y).toBeGreaterThanOrEqual(0);
  expect(bounds!.x + bounds!.width).toBeLessThanOrEqual(viewport.width + 1);
  expect(bounds!.y + bounds!.height).toBeLessThanOrEqual(viewport.height + 1);
}

for (const viewport of [
  { width: 1440, height: 1000 },
  { width: 1280, height: 600 },
  { width: 547, height: 638 },
  { width: 390, height: 844 },
]) {
  test(`fullscreen keeps the toolbar and its actions usable at ${viewport.width}×${viewport.height}`, async ({
    page,
  }, info) => {
    await page.setViewportSize(viewport);
    const errors: string[] = [];
    page.on('pageerror', (error) => errors.push(error.message));
    await page.goto('/');
    await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
    await page.getByRole('button', { name: 'Full screen', exact: true }).click();
    await expect
      .poll(() => page.evaluate(() => document.fullscreenElement?.classList.contains('app-shell')))
      .toBe(true);
    await expect(
      page.getByRole('button', { name: 'Exit full screen', exact: true }),
    ).toHaveAttribute('aria-pressed', 'true');
    const open = page.getByRole('button', { name: 'Open files', exact: true });
    const simulate = page.getByRole('button', { name: 'New simulation', exact: true });
    await expectOnscreen(page, open);
    await expectOnscreen(page, simulate);
    await expectOnscreen(page, page.getByRole('link', { name: 'DynaMol home' }));
    const canvas = await page.locator('.molecular-viewer__canvas canvas').boundingBox();
    expect(canvas!.width).toBeGreaterThan(100);
    expect(canvas!.height).toBeGreaterThan(100);
    await info.attach('fullscreen-toolbar', {
      body: await page.screenshot(),
      contentType: 'image/png',
    });

    await open.click();
    const importer = page.getByRole('dialog', { name: 'Bring your molecules.' });
    await expect(importer).toBeVisible();
    expect(await importer.evaluate((node) => document.fullscreenElement?.contains(node))).toBe(
      true,
    );
    await expectOnscreen(page, page.getByRole('button', { name: 'Close import', exact: true }));
    await page.getByRole('button', { name: 'Close import', exact: true }).click();
    await simulate.click();
    const studio = page.getByRole('dialog', { name: 'Set molecules in motion.' });
    await expect(studio).toBeVisible();
    expect(await studio.evaluate((node) => document.fullscreenElement?.contains(node))).toBe(true);
    await expectOnscreen(page, open);
    await page.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
    await page.getByRole('button', { name: 'Exit full screen', exact: true }).click();
    await expect.poll(() => page.evaluate(() => document.fullscreenElement === null)).toBe(true);
    await expect(page.getByRole('button', { name: 'Full screen', exact: true })).toHaveAttribute(
      'aria-pressed',
      'false',
    );
    await expectOnscreen(page, open);
    await expectOnscreen(page, simulate);
    expect(errors).toEqual([]);
  });
}

test('fullscreen state follows browser exits and handles an absent API', async ({ page }) => {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.getByRole('button', { name: 'Full screen', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Exit full screen', exact: true })).toBeVisible();
  await page.evaluate(() => document.exitFullscreen());
  await expect(page.getByRole('button', { name: 'Full screen', exact: true })).toHaveAttribute(
    'aria-pressed',
    'false',
  );
  await page.evaluate(() =>
    Object.defineProperty(HTMLElement.prototype, 'requestFullscreen', {
      value: undefined,
      configurable: true,
    }),
  );
  await page.getByRole('button', { name: 'Full screen', exact: true }).click();
  await expect(page.getByRole('status')).toContainText(
    'Fullscreen is unavailable in this browser.',
  );
});
