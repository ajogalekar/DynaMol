import { test, expect, type Page, type TestInfo } from '@playwright/test';

type Scene = { width: number; height: number; cx: number; cy: number; count: number };

async function ready(page: Page) {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
  await page.mouse.move(3, 3);
  await page.waitForTimeout(900);
}

// Measure actual rendered molecular geometry without reading NGL/React internals.
async function scene(page: Page): Promise<Scene> {
  await page.waitForTimeout(160);
  const png = await page.locator('.molecular-viewer__canvas canvas').screenshot();
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
    let minX = canvas.width,
      maxX = 0,
      minY = canvas.height,
      maxY = 0;
    let totalX = 0,
      totalY = 0,
      count = 0;
    // Restrict to the central scene: screenshots composite separate DOM corner
    // labels and toolbars over the WebGL canvas, so those are not molecule pixels.
    for (let y = Math.ceil(canvas.height * 0.03); y < canvas.height * 0.9; y++)
      for (let x = Math.ceil(canvas.width * 0.2); x < canvas.width * 0.8; x++) {
        const i = (y * canvas.width + x) * 4;
        if (Math.max(pixels[i], pixels[i + 1], pixels[i + 2]) < 65) continue;
        minX = Math.min(minX, x);
        maxX = Math.max(maxX, x);
        minY = Math.min(minY, y);
        maxY = Math.max(maxY, y);
        totalX += x;
        totalY += y;
        count++;
      }
    return {
      width: maxX - minX + 1,
      height: maxY - minY + 1,
      cx: totalX / count,
      cy: totalY / count,
      count,
    };
  }, png.toString('base64'));
}
async function viewport(page: Page) {
  return page.evaluate(() => ({
    scale: window.visualViewport?.scale,
    width: innerWidth,
    height: innerHeight,
    devicePixelRatio,
    scrollX,
    scrollY,
  }));
}
async function center(page: Page) {
  const box = (await page.locator('.molecular-viewer__canvas canvas').boundingBox())!;
  return { x: box.x + box.width / 2, y: box.y + box.height / 2 };
}
function scaled(after: Scene, before: Scene, factor: number, tolerance = 0.055) {
  expect(before.count).toBeGreaterThan(1000);
  expect(after.width / before.width).toBeGreaterThan(factor - tolerance);
  expect(after.width / before.width).toBeLessThan(factor + tolerance);
  expect(after.height / before.height).toBeGreaterThan(factor - tolerance);
  expect(after.height / before.height).toBeLessThan(factor + tolerance);
}
async function report(info: TestInfo, states: Record<string, unknown>) {
  await info.attach('rendered-geometry.json', {
    body: JSON.stringify(states, null, 2),
    contentType: 'application/json',
  });
}

test('trackpad Ctrl+wheel zooms molecular geometry both ways without scaling the page', async ({
  page,
}, info) => {
  await ready(page);
  const originalViewport = await viewport(page);
  const initial = await scene(page);
  const point = await center(page);
  await page.mouse.move(point.x, point.y);
  await page.keyboard.down('Control');
  await page.mouse.wheel(0, -25);
  await page.keyboard.up('Control');
  const expanded = await scene(page);
  scaled(expanded, initial, Math.exp(0.25));
  await page.keyboard.down('Control');
  await page.mouse.wheel(0, 25);
  await page.keyboard.up('Control');
  const restored = await scene(page);
  scaled(restored, initial, 1);
  expect(await viewport(page)).toEqual(originalViewport);
  await report(info, { initial, expanded, restored, originalViewport });
});

test('WebKit cumulative gesture scales zoom once when duplicate Ctrl+wheel events arrive', async ({
  page,
}, info) => {
  await ready(page);
  const originalViewport = await viewport(page);
  const initial = await scene(page);
  // Chrome cannot originate Safari GestureEvent: test its DOM protocol against
  // actual rendering. Physical macOS WebKit trackpad input remains manual.
  await page.locator('.molecular-viewer__canvas canvas').evaluate((target) => {
    const gesture = (type: string, scale: number) => {
      const event = new Event(type, { bubbles: true, cancelable: true });
      Object.defineProperty(event, 'scale', { value: scale });
      target.dispatchEvent(event);
      if (!event.defaultPrevented) throw new Error(`${type} did not consume page pinch`);
    };
    gesture('gesturestart', 1);
    for (const scale of [1.1, 1.2, 1.3]) {
      gesture('gesturechange', scale);
      target.dispatchEvent(
        new WheelEvent('wheel', { bubbles: true, cancelable: true, ctrlKey: true, deltaY: -20 }),
      );
    }
    gesture('gestureend', 1.3);
  });
  const expanded = await scene(page);
  scaled(expanded, initial, 1.3);
  await page.locator('.molecular-viewer__canvas canvas').evaluate((target) => {
    for (const [type, scale] of [
      ['gesturestart', 1],
      ['gesturechange', 1 / 1.3],
      ['gestureend', 1 / 1.3],
    ] as const) {
      const event = new Event(type, { bubbles: true, cancelable: true });
      Object.defineProperty(event, 'scale', { value: scale });
      target.dispatchEvent(event);
    }
  });
  const restored = await scene(page);
  scaled(restored, initial, 1);
  expect(await viewport(page)).toEqual(originalViewport);
  await report(info, { initial, expanded, restored, originalViewport });
});

test('slow touchscreen pinch continues zooming while its center drifts and pans', async ({
  page,
  context,
}, info) => {
  await ready(page);
  const client = await context.newCDPSession(page);
  await client.send('Emulation.setTouchEmulationEnabled', { enabled: true, maxTouchPoints: 2 });
  const originalViewport = await viewport(page);
  const initial = await scene(page);
  const origin = await center(page);
  const points = (distance: number, shift: number) => [
    { x: origin.x + shift - distance / 2, y: origin.y + shift / 2, id: 1 },
    { x: origin.x + shift + distance / 2, y: origin.y + shift / 2, id: 2 },
  ];
  await client.send('Input.dispatchTouchEvent', { type: 'touchStart', touchPoints: points(80, 0) });
  // Separation changes only 1.5 px/event, below NGL's old 2 px threshold.
  for (let step = 1; step <= 24; step++) {
    await client.send('Input.dispatchTouchEvent', {
      type: 'touchMove',
      touchPoints: points(80 + step * 1.5, step),
    });
  }
  await client.send('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] });
  const expanded = await scene(page);
  scaled(expanded, initial, 116 / 80);
  expect(expanded.cx - initial.cx).toBeGreaterThan(15);
  expect(expanded.cy - initial.cy).toBeGreaterThan(5);
  await client.send('Input.dispatchTouchEvent', {
    type: 'touchStart',
    touchPoints: points(116, 24),
  });
  for (let step = 1; step <= 24; step++) {
    await client.send('Input.dispatchTouchEvent', {
      type: 'touchMove',
      touchPoints: points(116 - step * 1.5, 24 - step),
    });
  }
  await client.send('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] });
  const restored = await scene(page);
  scaled(restored, initial, 1);
  expect(Math.abs(restored.cx - initial.cx)).toBeLessThan(4);
  expect(Math.abs(restored.cy - initial.cy)).toBeLessThan(4);
  expect(await viewport(page)).toEqual(originalViewport);
  await client.detach();
  await report(info, { initial, expanded, restored, originalViewport });
});

test('ordinary unmodified wheel still zooms both ways', async ({ page }, info) => {
  await ready(page);
  const initial = await scene(page);
  const point = await center(page);
  await page.mouse.move(point.x, point.y);
  await page.mouse.wheel(0, -120);
  const expanded = await scene(page);
  expect(expanded.width).toBeGreaterThan(initial.width * 1.04);
  expect(expanded.height).toBeGreaterThan(initial.height * 1.04);
  await page.mouse.wheel(0, 120);
  const contracted = await scene(page);
  expect(contracted.width).toBeLessThan(expanded.width * 0.98);
  expect(contracted.height).toBeLessThan(expanded.height * 0.98);
  await report(info, { initial, expanded, contracted });
});
