import { expect, test } from '@playwright/test';

const WORKSPACES = [
  'Crystal Viewer',
  'TEM Solver',
  'Diffraction',
  'XRD',
  'Variants',
  'Texture',
  'Calculator',
];

async function openWorkbench(page) {
  const browserErrors = [];
  page.on('pageerror', (error) => browserErrors.push(`page: ${error.message}`));
  page.on('console', (message) => {
    if (message.type() === 'error') browserErrors.push(`console: ${message.text()}`);
  });

  await page.goto('/');
  await expect(page).toHaveTitle('PyTex Workbench');
  await expect(page.getByRole('tab')).toHaveCount(WORKSPACES.length);
  await expect(page.getByRole('tab', { selected: true })).toHaveText('Crystal Viewer');
  return browserErrors;
}

async function expectNewCompletedCalculation(page, action) {
  const completed = page.locator('.activity__entry--completed');
  const before = await completed.count();
  await action();
  await expect.poll(() => completed.count(), { timeout: 20_000 }).toBeGreaterThan(before);
  await expect(page.locator('.activity__entry--failed')).toHaveCount(0);
}

test('loads every scientific workspace without browser errors', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await expect(page.getByRole('tab')).toHaveText(WORKSPACES);

  for (const workspace of WORKSPACES) {
    await page.getByRole('tab', { name: workspace, exact: true }).click();
    await expect(page.getByRole('tab', { name: workspace, exact: true })).toHaveAttribute(
      'aria-selected',
      'true',
    );
    await expect(page.locator('#stage')).not.toBeEmpty();
    await expect(page.locator('#rail-body')).not.toBeEmpty();
  }

  expect(browserErrors).toEqual([]);
});

test('completes the critical default calculations across panels', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  const journeys = [
    ['Crystal Viewer', 'Build structure'],
    ['Diffraction', 'Simulate pattern'],
    ['XRD', 'Simulate XRD pattern'],
    ['Variants', 'Show variants'],
    ['Texture', 'Build texture'],
  ];

  for (const [workspace, action] of journeys) {
    await page.getByRole('tab', { name: workspace, exact: true }).click();
    await expectNewCompletedCalculation(page, () =>
      page.getByRole('button', { name: action, exact: true }).click(),
    );
  }

  await page.getByRole('tab', { name: 'TEM Solver', exact: true }).click();
  await page.getByRole('button', { name: 'Auto-pick', exact: true }).click();
  await expectNewCompletedCalculation(page, () =>
    page.getByRole('button', { name: 'Index the pattern', exact: true }).click(),
  );

  await page.getByRole('tab', { name: 'Calculator', exact: true }).click();
  await expectNewCompletedCalculation(page, () =>
    page.getByRole('button', { name: 'Calculate', exact: true }).click(),
  );

  expect(browserErrors).toEqual([]);
});

test('surfaces a service failure to both the activity log and the user', async ({ page }) => {
  await openWorkbench(page);
  await page.route('**/api/call', async (route) => {
    const request = route.request().postDataJSON();
    if (request.operation?.startsWith('calc.')) {
      await route.fulfill({
        status: 400,
        contentType: 'application/json',
        body: JSON.stringify({
          ok: false,
          error: {
            code: 'browser_test.injected',
            message: 'Synthetic browser-test failure.',
            hint: 'This response is injected by the Playwright error-path test.',
          },
        }),
      });
      return;
    }
    await route.continue();
  });

  await page.getByRole('tab', { name: 'Calculator', exact: true }).click();
  await expect(page.locator('.toast')).toContainText('Synthetic browser-test failure.');
  await expect(page.locator('.activity__entry--failed')).toHaveCount(1);
  await expect(page.locator('.activity__entry--failed')).toContainText(
    'Synthetic browser-test failure.',
  );
});

test('keeps all workspaces reachable in the narrow responsive layout', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  const browserErrors = await openWorkbench(page);

  for (const workspace of WORKSPACES) {
    await expect(page.getByRole('tab', { name: workspace, exact: true })).toBeVisible();
  }
  await expect(page.locator('.masthead__action-label').first()).toBeHidden();
  const overflow = await page.evaluate(() => ({
    clientWidth: document.documentElement.clientWidth,
    scrollWidth: document.documentElement.scrollWidth,
  }));
  expect(overflow.scrollWidth).toBeLessThanOrEqual(overflow.clientWidth + 1);

  await page.getByRole('tab', { name: 'Calculator', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Calculate', exact: true })).toBeVisible();
  expect(browserErrors).toEqual([]);
});
