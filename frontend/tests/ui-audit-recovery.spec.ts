import { test } from './testWorkspace';
import { expect } from '@playwright/test';

test('audit recovery: unavailable initial example exposes import and a working fallback example', async ({
  page,
}) => {
  await page.route('**/api/datasets/demo', (route) => route.abort('failed'));
  await page.goto('/');
  await expect(page.getByRole('alert')).toBeVisible();
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeDisabled();
  await expect(page.getByRole('button', { name: /Try the ubiquitin example/ })).toBeVisible();
  await page.getByRole('button', { name: 'Open molecular files', exact: true }).click();
  await expect(page.getByRole('dialog', { name: 'Bring your molecules.' })).toBeVisible();
  await page.getByRole('button', { name: 'Close import', exact: true }).click();
  await page.unroute('**/api/datasets/demo');
  await page.getByRole('button', { name: /Try the ubiquitin example/ }).click();
  await expect(page.locator('.structure-card h2')).toContainText('Ubiquitin');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.getByRole('alert')).toHaveCount(0);
});

test('audit recovery: unavailable fullscreen reports a dismissible notification', async ({
  page,
}) => {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  // Browser capability fault injection; the app must handle rejected permission.
  await page.evaluate(() => {
    HTMLElement.prototype.requestFullscreen = () =>
      Promise.reject(new Error('Unavailable in audit browser'));
  });
  await page.getByRole('button', { name: 'Full screen', exact: true }).click();
  await expect(page.locator('.toast')).toContainText(
    'Fullscreen is unavailable in this browser.',
  );
  await page.getByRole('button', { name: 'Dismiss notification', exact: true }).click();
  await expect(page.locator('.toast')).toHaveCount(0);
});

test('audit recovery: graphics context loss blocks snapshots and reload restores the scene', async ({
  page,
}) => {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  // Dispatch the browser's context-loss event without altering molecular/API data.
  await page.locator('.molecular-viewer__canvas canvas').dispatchEvent('webglcontextlost');
  await expect(page.locator('.molecular-viewer__status--error')).toContainText(
    'graphics context was lost',
  );
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeDisabled();
  await page.reload();
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.getByRole('alert')).toHaveCount(0);
});

test('audit responsive: phone scene, inspector and Studio controls remain reachable', async ({
  page,
}, info) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  for (const label of ['Surface', 'Sticks', 'Ball & stick', 'Ribbons']) {
    await page.getByRole('button', { name: label, exact: true }).click();
    await expect(page.getByRole('button', { name: label, exact: true })).toHaveClass(/selected/);
  }
  await page.getByRole('button', { name: 'Polar only', exact: true }).click();
  await page.getByRole('textbox', { name: 'Find an atom', exact: true }).fill('LYS6');
  await expect(page.locator('.atom-search-results').getByRole('button').first()).toBeVisible();
  await page.getByRole('button', { name: 'Clear atom search', exact: true }).click();
  await page.getByRole('button', { name: 'Close measurement inspector', exact: true }).click();
  await page.getByRole('button', { name: 'Toggle measurement inspector', exact: true }).click();
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  const dialog = page.getByRole('dialog', { name: 'Set molecules in motion.' });
  await dialog.getByRole('tab', { name: 'Fetch', exact: true }).click();
  await dialog.getByRole('textbox', { name: 'PDB accession', exact: true }).fill('1UBQ');
  await dialog.getByRole('button', { name: /Advanced controls/ }).click();
  await dialog.getByRole('spinbutton', { name: /^Random seed/ }).fill('2026');
  await expect(dialog.getByRole('spinbutton', { name: /^Random seed/ })).toHaveValue('2026');
  expect(
    await page.evaluate(() => document.documentElement.scrollWidth - window.innerWidth),
  ).toBeLessThanOrEqual(2);
  await info.attach('phone-studio.png', {
    body: await page.screenshot({ fullPage: true }),
    contentType: 'image/png',
  });
  await dialog.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
  await page.getByRole('button', { name: 'Play trajectory', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Pause trajectory', exact: true })).toBeVisible();
  await page.getByRole('button', { name: 'Pause trajectory', exact: true }).click();
});
