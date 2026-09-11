import { defineConfig } from '@playwright/test';

/** Run the real local API and Vite first (./start.sh); no API or viewer mocks. */
export default defineConfig({
  testDir: './tests',
  testMatch: '**/*.spec.ts',
  fullyParallel: false,
  workers: 1,
  timeout: 90_000,
  expect: { timeout: 15_000 },
  use: {
    baseURL: process.env.DYNAMOL_BASE_URL || 'http://127.0.0.1:5173',
    channel: 'chrome',
    headless: true,
    viewport: { width: 1440, height: 1000 },
    screenshot: 'only-on-failure',
    trace: 'retain-on-failure',
    acceptDownloads: true,
  },
  reporter: [['list'], ['json', { outputFile: 'test-results/e2e-results.json' }]],
});
