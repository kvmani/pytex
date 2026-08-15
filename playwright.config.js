import { defineConfig, devices } from '@playwright/test';

/*
 * `PYTEX_BASE_URL` points the suite at a workbench that is already running, and
 * suppresses the managed one.
 *
 * `reuseExistingServer` alone is not enough: it reuses whatever happens to hold
 * the port, including a server started before the working tree changed, and the
 * suite then tests an old build while reporting on the new one. Naming the URL
 * makes that choice explicit rather than incidental.
 */
const baseURL = process.env.PYTEX_BASE_URL ?? 'http://127.0.0.1:8765';

export default defineConfig({
  testDir: './tests/browser',
  fullyParallel: false,
  forbidOnly: Boolean(process.env.CI),
  retries: process.env.CI ? 1 : 0,
  workers: process.env.CI ? 1 : undefined,
  reporter: process.env.CI ? [['line'], ['html', { open: 'never' }]] : 'line',
  use: {
    baseURL,
    trace: 'retain-on-failure',
    screenshot: 'only-on-failure',
  },
  projects: [
    {
      name: 'chromium',
      use: { ...devices['Desktop Chrome'] },
    },
  ],
  webServer: process.env.PYTEX_BASE_URL
    ? undefined
    : {
        command: 'python -m pytex.app --log-level WARNING serve --host 127.0.0.1 --port 8765',
        url: 'http://127.0.0.1:8765/api/manifest',
        reuseExistingServer: !process.env.CI,
        timeout: 60_000,
      },
});
