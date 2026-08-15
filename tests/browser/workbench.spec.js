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

/*
 * The figure has to be whole on arrival.
 *
 * The failure this guards against is silent: an SVG with a square viewBox on a
 * wide stage takes its width as its height, so the card grows past the bottom
 * of the window and the first thing a user does on every panel is scroll. It
 * looks fine in a screenshot of the top of the page, which is why it survived.
 */
test('shows every figure and its controls without scrolling the stage', async ({ page }) => {
  await page.setViewportSize({ width: 1280, height: 720 });
  const browserErrors = await openWorkbench(page);

  for (const workspace of WORKSPACES.filter((name) => name !== 'Calculator')) {
    await page.getByRole('tab', { name: workspace, exact: true }).click();
    const plot = page.locator('#stage .plot').first();
    await expect(plot).toBeVisible();
    // Every panel runs its first example on mount; the card is only worth
    // measuring once that drawing has arrived in it.
    await expect(plot.locator('svg').first()).toBeVisible({ timeout: 20_000 });
    const fits = await page.evaluate(() => {
      const stage = document.getElementById('stage');
      const card = stage.querySelector('.plot');
      const controls = card.querySelector('.plot__controls');
      const visible = controls && getComputedStyle(controls).display !== 'none';
      const bottom = (node) => Math.round(node.getBoundingClientRect().bottom);
      return {
        card: bottom(card),
        controls: visible ? bottom(controls) : null,
        stage: Math.round(stage.getBoundingClientRect().bottom),
        drawing: Math.round(card.querySelector('svg')?.getBoundingClientRect().height ?? 0),
      };
    });
    expect(fits.card).toBeLessThanOrEqual(fits.stage);
    if (fits.controls !== null) expect(fits.controls).toBeLessThanOrEqual(fits.stage);
    // A card that fits because the drawing collapsed is not a pass.
    expect(fits.drawing).toBeGreaterThan(150);
  }

  expect(browserErrors).toEqual([]);
});

test('zooms below 100% and pans with the pan tool', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await page.getByRole('tab', { name: 'Variants', exact: true }).click();
  await expect(page.locator('#stage .plot svg').first()).toBeVisible({ timeout: 20_000 });
  const zoom = page.locator('.plot__zoom');
  await expect(zoom).toHaveText('100%');

  const box = () => page.locator('#stage .plot svg').getAttribute('viewBox');
  const fitted = await box();

  const zoomOut = page.getByRole('button', { name: 'Zoom out', exact: true });
  for (let click = 0; click < 3; click += 1) await zoomOut.click();
  // Zoom out used to stop at Fit; it now runs below it.
  const percent = Number((await zoom.textContent()).replace('%', ''));
  expect(percent).toBeLessThan(100);

  const before = await box();
  const pan = page.getByRole('button', { name: 'Pan tool', exact: true });
  await pan.click();
  await expect(pan).toHaveAttribute('aria-pressed', 'true');
  const surface = page.locator('#stage .plot svg');
  const rect = await surface.boundingBox();
  await page.mouse.move(rect.x + rect.width / 2, rect.y + rect.height / 2);
  await page.mouse.down();
  await page.mouse.move(rect.x + rect.width / 2 + 90, rect.y + rect.height / 2 + 50, { steps: 6 });
  await page.mouse.up();
  expect(await box()).not.toBe(before);

  await page.getByRole('button', { name: 'Fit', exact: true }).click();
  await expect(zoom).toHaveText('100%');
  // Fit is the whole drawing again, whatever the zoom and pan did in between.
  expect(await box()).toBe(fitted);

  expect(browserErrors).toEqual([]);
});

test('reads the picked spots off the TEM pattern itself', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await page.getByRole('tab', { name: 'TEM Solver', exact: true }).click();
  const overlay = page.locator('.plot__overlay');
  // Nothing picked yet, so there is nothing to report.
  await expect(overlay).toBeHidden();

  await page.getByRole('button', { name: 'Auto-pick', exact: true }).click();
  await expect(overlay).toBeVisible();
  await expect(overlay.locator('tbody tr').first()).toBeVisible();
  const rows = await overlay.locator('tbody tr').count();
  expect(rows).toBeGreaterThan(1);

  const first = overlay.locator('tbody tr').first();
  // The reference spot: its ratio to itself is 1 and it has no angle to itself.
  await expect(first.locator('td').nth(3)).toHaveText('1.000');
  await expect(first.locator('td').nth(4)).toHaveText('—');

  const second = overlay.locator('tbody tr').nth(1);
  const d = Number(await second.locator('td').nth(2).innerText());
  const ratio = Number(await second.locator('td').nth(3).innerText());
  const angle = Number(await second.locator('td').nth(4).innerText());
  expect(d).toBeGreaterThan(0);
  expect(ratio).toBeGreaterThan(0);
  expect(angle).toBeGreaterThanOrEqual(0);
  expect(angle).toBeLessThanOrEqual(180);

  // It reports the picks, so clearing them takes it away again.
  await page.getByRole('button', { name: 'Clear picks', exact: true }).click();
  await expect(overlay).toBeHidden();

  expect(browserErrors).toEqual([]);
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
