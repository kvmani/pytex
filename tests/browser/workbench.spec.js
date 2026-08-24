import { readFileSync } from 'node:fs';

import { expect, test } from '@playwright/test';

const WORKSPACES = [
  'Crystal Viewer',
  'TEM Analysis',
  'XRD',
  'EBSD',
  'Variants',
  'Texture',
  'Calculator',
];

/*
 * Every panel, and the tab path that reaches it.
 *
 * The four transmission-electron panels live under one workspace tab, so naming
 * a panel is no longer the same as naming a tab. One table states the
 * difference; `openPanel` walks it.
 */
const PANEL_PATH = {
  'IPF map': ['EBSD', 'IPF map'],
  GROD: ['EBSD', 'GROD'],
  KAM: ['EBSD', 'KAM'],
  'Scan summary': ['EBSD', 'Scan summary'],
  Distributions: ['EBSD', 'Distributions'],
  'Pole figures': ['EBSD', 'Pole figures'],
  'Kikuchi simulator': ['EBSD', 'Kikuchi simulator'],
  'ECCI workflow': ['EBSD', 'ECCI workflow'],
  'SAED Simulator': ['TEM Analysis', 'SAED Simulator'],
  'TEM Solver': ['TEM Analysis', 'TEM Solver'],
  CBED: ['TEM Analysis', 'CBED'],
  'Composite SAED': ['TEM Analysis', 'Composite SAED'],
};

/** Panels that draw a figure, in tab order. Calculator draws none. */
const FIGURE_PANELS = [
  'Crystal Viewer',
  'Kikuchi simulator',
  'SAED Simulator',
  'TEM Solver',
  'CBED',
  'Composite SAED',
  'XRD',
  'IPF map',
  'Variants',
  'Texture',
];

/** The top-level tab bar, which sub-tabs also live in the `tab` role of. */
function workspaceTabs(page) {
  return page.locator('#tabs').getByRole('tab');
}

function workspaceTab(page, name) {
  return page.locator('#tabs').getByRole('tab', { name, exact: true });
}

/** Open a panel by name, through its workspace tab and any sub-tab. */
async function openPanel(page, name) {
  const [workspace, sub] = PANEL_PATH[name] ?? [name, null];
  await workspaceTab(page, workspace).click();
  if (sub) await page.locator('#subtabs').getByRole('tab', { name: sub, exact: true }).click();
}

/*
 * The TEM stage carries two figures — the pattern and the stereogram beside it —
 * so a bare `.plot` selector now matches both. Every TEM assertion names the
 * pattern by the label its own SVG carries.
 */
const PATTERN_SVG = '#stage svg[aria-label="Diffraction pattern"]';
const PATTERN_PLOT = '#stage .plot:has(svg[aria-label="Diffraction pattern"])';

/** The pattern figure's own copy of a frame control (zoom, cursor, overlay). */
function patternControl(page, selector) {
  return page.locator(PATTERN_PLOT).locator(selector);
}

/** A toolbar button of the pattern figure. Zoom, pan and Fit exist on both. */
function patternButton(page, name) {
  return page.locator(PATTERN_PLOT).getByRole('button', { name, exact: true });
}

/**
 * Fill a Miller-index control, which is a grid of boxes rather than a text field.
 *
 * `rows` is a list of index rows: `[[1, 1, 1], [2, 0, 0]]`. Rows are added with
 * the control's own button when the parameter starts with fewer than are
 * wanted, so the helper drives the control the way a user does rather than
 * reaching past it.
 */
async function fillIndices(control, rows) {
  const existing = await control.locator('.indices__row').count();
  for (let extra = existing; extra < rows.length; extra += 1) {
    await control.getByRole('button', { name: '+ Add row' }).click();
  }
  for (const [rowIndex, indices] of rows.entries()) {
    const boxes = control.locator('.indices__row').nth(rowIndex).locator('.indices__box');
    for (const [boxIndex, value] of indices.entries()) {
      await boxes.nth(boxIndex).fill(String(value));
    }
  }
}

async function openWorkbench(page) {
  const browserErrors = [];
  page.on('pageerror', (error) => browserErrors.push(`page: ${error.message}`));
  page.on('console', (message) => {
    if (message.type() === 'error') browserErrors.push(`console: ${message.text()}`);
  });

  await page.goto('/');
  await expect(page).toHaveTitle('PyTex Workbench');
  await expect(workspaceTabs(page)).toHaveCount(WORKSPACES.length);
  await expect(page.locator('#tabs').getByRole('tab', { selected: true })).toHaveText(
    'Crystal Viewer',
  );
  await dismissTour(page);
  return browserErrors;
}

/*
 * Get past the welcome, which every test meets because Playwright gives each one
 * a fresh profile and the tour is remembered per browser.
 *
 * Skipping it here rather than disabling it in the fixture is deliberate: it
 * means every test in this file re-proves the property that actually matters
 * about a greeting — that one click removes it and the application underneath is
 * immediately usable. A tour that ever failed to get out of the way would fail
 * the whole suite rather than one test nobody runs.
 */
async function dismissTour(page) {
  const tour = page.locator('.tour');
  if (!(await tour.isVisible().catch(() => false))) return;
  await page.locator('.tour__skip').click();
  await expect(tour).toBeHidden();
}

/** Open the message console, which is collapsed on load. */
async function openConsole(page) {
  const toggle = page.locator('#console-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  await expect(page.locator('#console-panel')).toBeVisible();
}

async function expectNewCompletedCalculation(page, action) {
  await openConsole(page);
  const completed = page.locator('.console__entry--success');
  const before = await completed.count();
  await action();
  await expect.poll(() => completed.count(), { timeout: 20_000 }).toBeGreaterThan(before);
  await expect(page.locator('.console__entry--error')).toHaveCount(0);
}

test('greets a first-time visitor and gets out of the way in one click', async ({ page }) => {
  const browserErrors = [];
  page.on('pageerror', (error) => browserErrors.push(`page: ${error.message}`));
  page.on('console', (message) => {
    if (message.type() === 'error') browserErrors.push(`console: ${message.text()}`);
  });

  await page.goto('/');
  const tour = page.locator('.tour');
  await expect(tour).toBeVisible();
  await expect(page.locator('.tour__title')).toHaveText('Welcome to the PyTex Workbench');
  await expect(page.locator('.tour__progress')).toHaveText('Step 1 of 7');

  // Each step points at something that is really on the page. A step whose
  // target had gone missing would be a tour explaining furniture that is not
  // there, which is worse than no tour.
  const targets = [];
  for (let step = 2; step <= 7; step += 1) {
    await page.locator('.tour__actions .button--primary').click();
    await expect(page.locator('.tour__progress')).toHaveText(`Step ${step} of 7`);
    await expect(page.locator('.tour-target')).toHaveCount(1);
    targets.push(await page.locator('.tour-target').getAttribute('id'));
  }
  expect(targets).toEqual(['tabs', 'rail', 'stage', 'open-palette', 'console', 'open-feedback']);

  await page.locator('.tour__actions .button--primary').click();
  await expect(tour).toBeHidden();

  // Seen once, not on every load: the greeting is a first-visit event.
  await page.reload();
  await expect(workspaceTabs(page)).toHaveCount(WORKSPACES.length);
  await expect(tour).toBeHidden();

  // And still reachable on purpose, from the help panel, which is what makes
  // remembering it safe rather than final.
  await page.getByRole('button', { name: 'Help' }).click();
  await page.getByRole('button', { name: 'Show the welcome tour again' }).click();
  await expect(tour).toBeVisible();
  await page.locator('.tour__skip').click();
  await expect(tour).toBeHidden();

  expect(browserErrors).toEqual([]);
});

test('the tour can be skipped from the first step and never returns', async ({ page }) => {
  await page.goto('/');
  await expect(page.locator('.tour')).toBeVisible();
  await page.locator('.tour__skip').click();
  await expect(page.locator('.tour')).toBeHidden();

  // The application underneath must be immediately usable, not merely visible.
  await workspaceTab(page, 'Calculator').click();
  await expect(workspaceTab(page, 'Calculator')).toHaveAttribute('aria-selected', 'true');

  await page.reload();
  await expect(workspaceTabs(page)).toHaveCount(WORKSPACES.length);
  await expect(page.locator('.tour')).toBeHidden();
});

test('sends a feature request and says what became of it', async ({ page }) => {
  const browserErrors = await openWorkbench(page);

  await page.getByRole('button', { name: 'Feedback' }).click();
  const drawer = page.locator('#feedback-drawer');
  await expect(drawer).toBeVisible();
  await expect(drawer.getByRole('heading', { level: 2 })).toHaveText(
    'Tell us what would make PyTex better',
  );

  // An empty note is refused here rather than sent for the server to refuse,
  // and the refusal is an encouragement rather than a validation error.
  await drawer.getByRole('button', { name: 'Send it' }).click();
  await expect(drawer.locator('.field__error')).toContainText('even one line is useful');

  await drawer.locator('#feedback-category').selectOption('feature');
  await drawer
    .locator('#feedback-message')
    .fill('An inverse pole figure of the sample normal would save me an afternoon.');
  await drawer.locator('#feedback-name').fill('A Researcher');
  await drawer.locator('#feedback-email').fill('someone@example.invalid');
  await drawer.locator('#feedback-rating').selectOption('4');
  await drawer.getByRole('button', { name: 'Send it' }).click();

  await expect(drawer.getByRole('heading', { level: 2 })).toHaveText(
    'Thank you — it is on its way',
  );
  // The receipt distinguishes "filed here" from "filed and e-mailed" rather
  // than showing one cheerful message for both.
  await expect(drawer.locator('.field__help')).toContainText('filed locally');
  await expect(drawer.getByRole('button', { name: 'Send another' })).toBeVisible();

  expect(browserErrors).toEqual([]);
});

test('a half-written note survives closing the form', async ({ page }) => {
  await openWorkbench(page);
  const drawer = page.locator('#feedback-drawer');

  await page.getByRole('button', { name: 'Feedback' }).click();
  await drawer.locator('#feedback-message').fill('The legend on the pole figure is');
  await drawer.getByRole('button', { name: 'Not now' }).click();
  await expect(drawer).toBeHidden();

  // People start a note, go back to look at the thing they are describing, and
  // come back. Losing the text at that point loses the report: nobody types it
  // twice.
  await page.getByRole('button', { name: 'Feedback' }).click();
  await expect(drawer.locator('#feedback-message')).toHaveValue(
    'The legend on the pole figure is',
  );
});

test('Miller indices are typed one index to a box', async ({ page }) => {
  const browserErrors = await openWorkbench(page);

  const planes = page.locator('#rail-body .indices--multi').first();
  await expect(planes).toHaveAttribute('aria-label', 'Planes to superimpose (hkl)');

  // Three boxes, each named for the index it holds. Three unlabelled fields
  // would be the same ambiguity in another form.
  const firstRow = planes.locator('.indices__row').first();
  await expect(firstRow.locator('.indices__box')).toHaveCount(3);
  await expect(firstRow.locator('.indices__box').nth(0)).toHaveAttribute('placeholder', 'h');
  await expect(firstRow.locator('.indices__box').nth(1)).toHaveAttribute('placeholder', 'k');
  await expect(firstRow.locator('.indices__box').nth(2)).toHaveAttribute('placeholder', 'l');

  // There is nowhere on the page to put a run of digits: a box takes one index.
  const rowsBefore = await planes.locator('.indices__row').count();
  await planes.getByRole('button', { name: '+ Add row' }).click();
  await expect(planes.locator('.indices__row')).toHaveCount(rowsBefore + 1);

  const added = planes.locator('.indices__row').nth(rowsBefore);
  await added.locator('.indices__box').nth(0).fill('2');
  // A row started and not finished is named for the indices still missing,
  // which the server cannot phrase as well: by then it is one list with a hole.
  await expect(planes.locator('xpath=ancestor::div[@class="field"]').locator('.field__error'))
    .toContainText('k, l are still empty');

  await added.locator('.indices__box').nth(1).fill('0');
  await added.locator('.indices__box').nth(2).fill('0');
  await expect(planes.locator('xpath=ancestor::div[@class="field"]').locator('.field__error'))
    .toBeHidden();

  await expectNewCompletedCalculation(page, () =>
    page.locator('#rail-body').getByRole('button', { name: 'Build structure' }).click(),
  );

  expect(browserErrors).toEqual([]);
});

test('a direction parameter names its boxes u, v and w', async ({ page }) => {
  await openWorkbench(page);
  const directions = page.locator('#rail-body .indices[aria-label*="[uvw]"]').first();
  const row = directions.locator('.indices__row').first();
  await expect(row.locator('.indices__box').nth(0)).toHaveAttribute('placeholder', 'u');
  await expect(row.locator('.indices__box').nth(1)).toHaveAttribute('placeholder', 'v');
  await expect(row.locator('.indices__box').nth(2)).toHaveAttribute('placeholder', 'w');
});

test('all three index boxes fit inside the rail', async ({ page }) => {
  await openWorkbench(page);
  const directions = page.locator('#rail-body .indices[aria-label*="[uvw]"]').first();
  const row = directions.locator('.indices__row').first();
  const boxes = row.locator('.indices__box');
  await expect(boxes).toHaveCount(3);

  // The base rule for a text input is `width: 100%`, and an attribute selector
  // outweighs a bare class: when that rule wins, the first box fills the rail
  // and the other two are pushed off the edge of the screen, unreachable. Each
  // box is therefore narrow, and the third one ends inside the rail.
  const rail = await page.locator('#rail-body').boundingBox();
  for (let index = 0; index < 3; index += 1) {
    const box = await boxes.nth(index).boundingBox();
    expect(box.width).toBeLessThan(72);
    expect(box.x + box.width).toBeLessThanOrEqual(rail.x + rail.width);
  }

  // Four characters - three digits and a sign - fit without scrolling.
  await boxes.nth(0).fill('-100');
  const overflow = await boxes.nth(0).evaluate((node) => node.scrollWidth - node.clientWidth);
  expect(overflow).toBeLessThanOrEqual(1);
});

test('loads every scientific workspace without browser errors', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await expect(workspaceTabs(page)).toHaveText(WORKSPACES);

  for (const workspace of WORKSPACES) {
    await workspaceTab(page, workspace).click();
    await expect(workspaceTab(page, workspace)).toHaveAttribute('aria-selected', 'true');
    await expect(page.locator('#stage')).not.toBeEmpty();
    await expect(page.locator('#rail-body')).not.toBeEmpty();
  }

  // And every panel of the grouped workspace, which the tab bar no longer names.
  for (const panel of Object.keys(PANEL_PATH)) {
    await openPanel(page, panel);
    await expect(page.locator('#subtabs').getByRole('tab', { name: panel, exact: true })).toHaveAttribute(
      'aria-selected',
      'true',
    );
    await expect(page.locator('#stage')).not.toBeEmpty();
  }

  expect(browserErrors).toEqual([]);
});

test('offers the shared CIF phase loader in every structure-aware workspace', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  const structurePanels = [
    'Crystal Viewer',
    'SAED Simulator',
    'XRD',
    'Kikuchi simulator',
    'ECCI workflow',
    'Variants',
    'Texture',
    'Calculator',
  ];

  for (const panel of structurePanels) {
    await openPanel(page, panel);
    await expect(page.getByLabel('Load a CIF crystal structure').first()).toBeVisible();
  }
  expect(browserErrors).toEqual([]);

  // The browser carries the file as phase input; Python owns extension checks
  // and CIF parsing. This deliberately bypasses the file dialog's accept hint
  // to prove that server-side validation still guards every shared picker.
  await openPanel(page, 'XRD');
  await page.getByLabel('Load a CIF crystal structure').setInputFiles({
    name: 'not-a-structure.txt',
    mimeType: 'text/plain',
    buffer: Buffer.from('this is not CIF'),
  });
  await expect(page.locator('.phase-cif')).toContainText('not-a-structure.txt loaded');
  await page.getByRole('button', { name: 'Simulate XRD pattern', exact: true }).click();
  await expect(page.locator('#rail-body .field__error:visible')).toContainText('extension is .txt');

});

test('completes the critical default calculations across panels', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  const journeys = [
    ['Crystal Viewer', 'Build structure'],
    ['Composite SAED', 'Simulate pattern'],
    ['XRD', 'Simulate XRD pattern'],
    ['Variants', 'Show variants'],
    ['Texture', 'Build texture'],
  ];

  for (const [workspace, action] of journeys) {
    await openPanel(page, workspace);
    await expectNewCompletedCalculation(page, () =>
      page.getByRole('button', { name: action, exact: true }).click(),
    );
  }

  await openPanel(page, 'TEM Solver');
  await page.getByRole('button', { name: 'Auto-pick', exact: true }).click();
  await expectNewCompletedCalculation(page, () =>
    page.getByRole('button', { name: 'Index the pattern', exact: true }).click(),
  );

  await workspaceTab(page, 'Calculator').click();
  await expectNewCompletedCalculation(page, () =>
    page.getByRole('button', { name: 'Calculate', exact: true }).click(),
  );

  expect(browserErrors).toEqual([]);
});

test('surfaces a service failure to both the message log and the user', async ({ page }) => {
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
          // The envelope's own narration is what the console renders, so the
          // injected failure must carry one exactly as the server would.
          log: [
            {
              sequence: 90001,
              time: Date.now() / 1000,
              level: 'error',
              message: 'Synthetic browser-test failure.',
              source: request.operation,
            },
          ],
        }),
      });
      return;
    }
    await route.continue();
  });

  await workspaceTab(page, 'Calculator').click();
  await expect(page.locator('.toast')).toContainText('Synthetic browser-test failure.');
  await openConsole(page);
  await expect(page.locator('.console__entry--error')).toHaveCount(1);
  await expect(page.locator('.console__entry--error')).toContainText(
    'Synthetic browser-test failure.',
  );
});

test('the console narrates a session and filters it by severity', async ({ page }) => {
  await openWorkbench(page);

  // Collapsed, the bar still reports: a user who never opens the console must
  // still see that something happened. It shows the *newest* message, so the
  // assertion is that it is not idle rather than that it holds any one line —
  // by the time the first panel has drawn, the start-up notice is history.
  await expect(page.locator('#console-toggle')).not.toContainText('Ready');
  await expect(page.locator('#console-toggle')).toContainText('messages');

  await openConsole(page);
  const stream = page.locator('#console-stream');
  await expect(stream.locator('.console__entry')).not.toHaveCount(0);
  // The start-up notice lives in the stream, where nothing displaces it.
  await expect(stream).toContainText('ready in the web shell');

  // Every entry carries a time, a severity mark and the surface that reported
  // it, so a message can be traced back without guessing.
  const first = stream.locator('.console__entry').first();
  await expect(first.locator('.console__time')).not.toBeEmpty();
  await expect(first.locator('.console__source')).not.toBeEmpty();

  await workspaceTab(page, 'Calculator').click();
  await expect(stream).toContainText('Opened the Calculator workspace.');

  // Errors only: the panel-switch note is info, so it must disappear.
  await page.locator('#console-threshold').selectOption('40');
  await expect(stream).not.toContainText('Opened the Calculator workspace.');
  await page.locator('#console-threshold').selectOption('0');
  await expect(stream).toContainText('Opened the Calculator workspace.');

  // The text filter is what makes a long session searchable.
  await page.locator('#console-search').fill('workspace');
  await expect(stream.locator('.console__entry')).not.toHaveCount(0);
  await page.locator('#console-search').fill('no message says this');
  await expect(stream.locator('.console__entry')).toHaveCount(0);
  await expect(page.locator('.console__empty')).toContainText('No message matches this filter.');
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

  for (const workspace of FIGURE_PANELS) {
    await openPanel(page, workspace);
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

/*
 * The EBSD workspace's forward problem: the pattern behind every indexed point.
 *
 * Three things are worth pinning in a browser. The panel must draw a pattern
 * without a scan, because it is the one view here that needs none. Every band's
 * name must sit *on the screen* — the traces arrive clipped to a generous
 * margin so that a Kossel conic keeps its shape, and a name placed by a
 * fraction along the whole trace would be eaten by the clip without a trace of
 * its own. And the names must run along their bands, the convention every
 * Kikuchi figure in the application shares.
 */
test('the EBSD workspace simulates a Kikuchi pattern and names its bands along them', async ({
  page,
}) => {
  const browserErrors = await openWorkbench(page);
  await openPanel(page, 'Kikuchi simulator');

  const pattern = page.locator('svg[aria-label="Simulated EBSD Kikuchi pattern"]');
  await expect(pattern).toBeVisible({ timeout: 30_000 });
  await expect.poll(() => pattern.locator('polyline').count()).toBeGreaterThan(10);

  const status = page.locator('#stage .plot__status');
  await expect(status).toContainText('band(s)');
  await expect(status).toContainText('intensities kinematic');

  const labels = await page.evaluate(() => {
    const svg = document.querySelector('svg[aria-label="Simulated EBSD Kikuchi pattern"]');
    const box = svg.viewBox.baseVal;
    return [...svg.querySelectorAll('text')]
      .filter((node) => /^\(/.test(node.textContent ?? ''))
      .map((node) => ({
        x: Number(node.getAttribute('x')),
        y: Number(node.getAttribute('y')),
        angle: Number(/rotate\(([-\d.]+)/.exec(node.getAttribute('transform') ?? '')?.[1]),
        width: box.width,
        height: box.height,
      }));
  });

  expect(labels.length).toBeGreaterThan(3);
  for (const label of labels) {
    expect(label.x).toBeGreaterThan(0);
    expect(label.x).toBeLessThan(label.width);
    expect(label.y).toBeGreaterThan(0);
    expect(label.y).toBeLessThan(label.height);
    expect(Number.isFinite(label.angle)).toBe(true);
  }
  // A cube-oriented crystal on a tilted stage gives bands at several slopes, so
  // labels that were all upright would mean the rotation was not applied.
  expect(labels.some((label) => Math.abs(label.angle) > 5)).toBe(true);

  // The geometry is stated beside the picture, including the one number that
  // checks the whole convention.
  const readout = page.locator('#stage .plot__readout-panel');
  await expect(readout).toContainText('Specimen normal at');
  await expect(readout).toContainText('0.364');

  // Turning the names off takes them away rather than leaving them stale.
  await page.locator('#stage').getByRole('button', { name: 'Indices', exact: true }).click();
  await expect
    .poll(() =>
      page.evaluate(() => {
        const svg = document.querySelector('svg[aria-label="Simulated EBSD Kikuchi pattern"]');
        return [...svg.querySelectorAll('text')].filter((node) =>
          /^\(/.test(node.textContent ?? ''),
        ).length;
      }),
    )
    .toBe(0);

  expect(browserErrors).toEqual([]);
});

test('zooms below 100% and pans with the pan tool', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await workspaceTab(page, 'Variants').click();
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
  await openPanel(page, 'TEM Solver');

  /*
   * The numbers are read *while* the plate is worked on, so they are under the
   * drawing rather than over it.
   *
   * They used to be a card pinned to the top-left corner of the figure, painted
   * beneath it so the pattern masked it and a hover raised it. That is right for
   * a picture that is looked at and wrong for one that is picked on: the corner
   * it covered is a corner of the data, and the measurement was hidden behind
   * the very spots it measures. The readout bar can never overlap the drawing,
   * which is the property this asserts.
   */
  const readout = patternControl(page, '.plot__readout');
  await expect(readout).toBeVisible();
  await expect(readout).toContainText('Measured picks');
  // The live cursor lives there too, resting on a dash until the pointer moves.
  await expect(patternControl(page, '.plot__cursor')).toHaveText('—');

  await page.getByRole('button', { name: 'Auto-pick', exact: true }).click();
  await expect(readout.locator('tbody tr').first()).toBeVisible();
  const rows = await readout.locator('tbody tr').count();
  expect(rows).toBeGreaterThan(1);

  // #, R/px, |g|, d, ratio, angle — the three numbers a plate is identified by,
  // plus the two radii they come from.
  await expect(readout.locator('thead th')).toHaveText(['#', 'R / px', '|g| / Å⁻¹', 'd / Å', 'd1/d', '∠ / °']);

  const first = readout.locator('tbody tr').first();
  // The reference spot: its ratio to itself is 1 and it has no angle to itself.
  await expect(first.locator('td').nth(4)).toHaveText('1.000');
  await expect(first.locator('td').nth(5)).toHaveText('—');

  const second = readout.locator('tbody tr').nth(1);
  const g = Number(await second.locator('td').nth(2).innerText());
  const d = Number(await second.locator('td').nth(3).innerText());
  const ratio = Number(await second.locator('td').nth(4).innerText());
  const angle = Number(await second.locator('td').nth(5).innerText());
  expect(g).toBeGreaterThan(0);
  // |g| and d are the same measurement stated two ways, and must agree.
  expect(Math.abs(g * d - 1)).toBeLessThan(0.01);
  expect(ratio).toBeGreaterThan(0);
  expect(angle).toBeGreaterThanOrEqual(0);
  expect(angle).toBeLessThanOrEqual(180);

  // Nothing in the bar overlaps the drawing, at any pick count.
  const clear = await page.evaluate(() => {
    const plot = document
      .querySelector('svg[aria-label="Diffraction pattern"]')
      .closest('.plot');
    const drawing = plot.querySelector('.plot__stage').getBoundingClientRect();
    const bar = plot.querySelector('.plot__readout').getBoundingClientRect();
    const stage = document.getElementById('stage').getBoundingClientRect();
    return {
      below: bar.top >= drawing.bottom - 1,
      drawingHeight: drawing.height,
      insideStage: bar.bottom <= stage.bottom + 1,
    };
  });
  expect(clear.below).toBe(true);
  expect(clear.insideStage).toBe(true);
  // The figure still has room to be a figure.
  expect(clear.drawingHeight).toBeGreaterThan(150);

  // The cursor reports the radius from the beam and the |g| it corresponds to,
  // which is the reading taken while hovering rather than after clicking.
  const cursorText = await page.evaluate(() => {
    const svg = document.querySelector('svg[aria-label="Diffraction pattern"]');
    const box = svg.getBoundingClientRect();
    svg.dispatchEvent(
      new PointerEvent('pointermove', {
        bubbles: true,
        clientX: box.left + box.width * 0.7,
        clientY: box.top + box.height * 0.5,
      }),
    );
    return svg.closest('.plot').querySelector('.plot__cursor').textContent;
  });
  expect(cursorText).toContain('px from beam');
  expect(cursorText).toContain('|g|');

  // It reports the picks, so clearing them returns it to its waiting state.
  await page.getByRole('button', { name: 'Clear picks', exact: true }).click();
  await expect(readout).toContainText('Click the transmitted beam');

  expect(browserErrors).toEqual([]);
});

/* ------------------------------------------------ picking the TEM pattern
 *
 * Three failures these pin, all of them things a user saw on screen and could
 * not have diagnosed:
 *
 * 1. A pick made after zooming or panning landed off by the camera offset,
 *    because the panel converted the pointer itself instead of asking the frame
 *    that owns the viewBox.
 * 2. Two picks produced a lattice through neither of them, because a
 *    rank-deficient least squares answered anyway.
 * 3. A typed coordinate reverted, because focusing the field rebuilt the table
 *    the field was in.
 */

/** Open a practice plate and return a screen<->image coordinate converter. */
async function openPlate(page) {
  await openPanel(page, 'TEM Solver');
  const surface = page.locator(PATTERN_SVG);
  await expect(surface).toBeVisible({ timeout: 20_000 });
  await expect(patternControl(page, '.plot__zoom')).toHaveText('100%');
  return async (x, y) => {
    const geometry = await page.evaluate(() => {
      const svg = document.querySelector('svg[aria-label="Diffraction pattern"]');
      const box = svg.viewBox.baseVal;
      const rect = svg.getBoundingClientRect();
      const scale = Math.min(rect.width / box.width, rect.height / box.height);
      return {
        left: rect.left,
        top: rect.top,
        box: { x: box.x, y: box.y, width: box.width, height: box.height },
        scale,
        offsetX: (rect.width - box.width * scale) / 2,
        offsetY: (rect.height - box.height * scale) / 2,
      };
    });
    return {
      x: geometry.left + geometry.offsetX + (x - geometry.box.x) * geometry.scale,
      y: geometry.top + geometry.offsetY + (y - geometry.box.y) * geometry.scale,
    };
  };
}

/** Every coordinate the pick table is showing, beam first. */
async function pickedCoordinates(page) {
  return page.evaluate(() =>
    [...document.querySelectorAll('.picks__row')].map((row) =>
      [...row.querySelectorAll('.picks__input')].map((input) => Number(input.value)),
    ),
  );
}

test('a pick lands where the pointer is, after zooming and panning', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  const toScreen = await openPlate(page);

  // Move the camera first, then pick. At the fitted view every conversion
  // agrees; the offset is what separates a correct one from a plausible one.
  await patternButton(page, 'Zoom in').click();
  await patternButton(page, 'Zoom in').click();
  const pan = patternButton(page, 'Pan tool');
  await pan.click();
  const rect = await page.locator(PATTERN_SVG).boundingBox();
  await page.mouse.move(rect.x + rect.width / 2, rect.y + rect.height / 2);
  await page.mouse.down();
  await page.mouse.move(rect.x + rect.width / 2 - 70, rect.y + rect.height / 2 + 40, { steps: 6 });
  await page.mouse.up();
  await pan.click();
  await expect(pan).toHaveAttribute('aria-pressed', 'false');

  const target = { x: 512, y: 512 };
  const point = await toScreen(target.x, target.y);
  await page.mouse.click(point.x, point.y);
  await expect(page.locator('.picks__row')).toHaveCount(1);

  const [beam] = await pickedCoordinates(page);
  expect(beam[0]).toBeCloseTo(target.x, 0);
  expect(beam[1]).toBeCloseTo(target.y, 0);

  expect(browserErrors).toEqual([]);
});

test('the view stays where it was put while picking', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  const toScreen = await openPlate(page);
  await patternButton(page, 'Zoom in').click();
  const zoomed = await patternControl(page, '.plot__zoom').textContent();
  expect(Number(zoomed.replace('%', ''))).toBeGreaterThan(100);

  const point = await toScreen(512, 512);
  await page.mouse.click(point.x, point.y);
  await expect(page.locator('.picks__row')).toHaveCount(1);
  // The redraw that follows a pick used to snap the camera back to Fit, which
  // made zooming in to place a spot precisely impossible.
  await expect(patternControl(page, '.plot__zoom')).toHaveText(zoomed);

  expect(browserErrors).toEqual([]);
});

test('two picks lay the lattice through the two picks', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  const toScreen = await openPlate(page);

  // A beam and two reflections of this plate, clicked as a user would.
  const clicks = [
    [512, 512],
    [452.48, 709.3],
    [313.59, 452.73],
  ];
  for (const [x, y] of clicks) {
    const point = await toScreen(x, y);
    await page.mouse.click(point.x, point.y);
  }
  await expect(page.locator('.picks__row')).toHaveCount(3);
  const coordinates = await pickedCoordinates(page);
  for (const [index, [x, y]] of clicks.entries()) {
    expect(coordinates[index][0]).toBeCloseTo(x, 0);
    expect(coordinates[index][1]).toBeCloseTo(y, 0);
  }

  // The beam is where it was clicked and says so, rather than being quietly
  // replaced by a centre that two spots cannot determine.
  const rail = page.locator('.centre-tool');
  await expect(rail).toContainText('held where it was picked', { timeout: 10_000 });
  await expect(rail).toContainText('at least 4 spots');
  await expect(rail.getByRole('button', { name: 'Refine beam from the spots' })).toBeDisabled();
  // Both basis arrows are on picked spots, not on empty nodes.
  await expect(rail).toContainText('a (spot 1');
  await expect(rail).toContainText('b (spot 2');
  await expect(rail).not.toContainText('no pick sits on this node');

  expect(browserErrors).toEqual([]);
});

/*
 * The overlay stops at the edge of the picture.
 *
 * The fitted lattice is generated by walking node indices outwards until the
 * grid certainly covers the image diagonal, so it always overshoots. A viewBox
 * does not contain that overshoot — with `meet` the drawing is letterboxed and
 * anything outside the box is painted over the margins — so the grid ran across
 * blank page beside the pattern, asserting lattice where there is no image.
 *
 * Asserted structurally, because that is where a regression would actually
 * appear: a later overlay appended to the outer element instead of to the
 * clipped group. The lines' own coordinates are checked to still run past the
 * image, so the test cannot pass by the overshoot quietly disappearing and
 * leaving the clip untested.
 */
test('the lattice overlay is clipped to the image', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await openPlate(page);
  await page.getByRole('button', { name: 'Auto-pick', exact: true }).click();
  // The lattice is on by default, so the grid follows the fit without a press.
  // The fit is a round trip, so it arrives after the picks. A handful of lines
  // is the beam crosshair alone; the grid is dozens.
  await expect
    .poll(() => page.locator(`${PATTERN_SVG} line`).count(), { timeout: 20_000 })
    .toBeGreaterThan(20);

  const geometry = await page.evaluate(() => {
    const svg = document.querySelector('svg[aria-label="Diffraction pattern"]');
    const box = svg.viewBox.baseVal;
    const clip = svg.querySelector('clipPath rect');
    let widest = 0;
    for (const line of svg.querySelectorAll('line')) {
      for (const attribute of ['x1', 'x2']) {
        const value = Number(line.getAttribute(attribute));
        widest = Math.max(widest, Math.abs(value - box.width / 2));
      }
    }
    return {
      box: { width: box.width, height: box.height },
      clip: clip && {
        x: Number(clip.getAttribute('x')),
        y: Number(clip.getAttribute('y')),
        width: Number(clip.getAttribute('width')),
        height: Number(clip.getAttribute('height')),
      },
      unclipped: [...svg.children]
        .filter((child) => child.tagName !== 'defs')
        .filter((child) => !child.getAttribute('clip-path'))
        .map((child) => child.tagName),
      widest,
    };
  });

  expect(geometry.clip).toEqual({ x: 0, y: 0, ...geometry.box });
  expect(geometry.unclipped).toEqual([]);
  // The grid really does run past the frame; the clip is what stops it.
  expect(geometry.widest).toBeGreaterThan(geometry.box.width / 2);

  expect(browserErrors).toEqual([]);
});

test('a coordinate typed into the table moves the pick', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await openPlate(page);
  await page.getByRole('button', { name: 'Auto-pick', exact: true }).click();
  await expect(page.locator('.picks__row').first()).toBeVisible();

  const beamX = page.locator('.picks__input').first();
  await expect(beamX).toHaveValue('512.00');
  await beamX.click();
  await beamX.press('Control+a');
  await beamX.type('460');
  await beamX.press('Enter');

  // It must survive the fit that lands a moment later, which used to rebuild
  // the table and put the old value back.
  await expect(beamX).toHaveValue('460.00');
  await page.waitForTimeout(800);
  await expect(beamX).toHaveValue('460.00');
  // With the beam moved off the spots, the fit says so and offers its own.
  await expect(page.locator('.centre-tool')).toContainText('fit says');

  // And the nudge pad moves the same pick, by the step it advertises.
  await page.getByRole('button', { name: 'Move the selected pick right' }).click();
  await expect(beamX).toHaveValue('461.00');

  expect(browserErrors).toEqual([]);
});

test('a pick can be renumbered, promoted and removed from the table', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await openPlate(page);
  await page.getByRole('button', { name: 'Auto-pick', exact: true }).click();
  await expect(page.locator('.picks__row')).toHaveCount(7);

  const before = await pickedCoordinates(page);
  await page.getByRole('button', { name: 'Make spot 1 the transmitted beam' }).click();
  const after = await pickedCoordinates(page);
  // A swap, not a move: the old beam is still a pick, in the row vacated.
  expect(after[0]).toEqual(before[1]);
  expect(after[1]).toEqual(before[0]);

  await page.getByRole('button', { name: 'Remove spot 1' }).click();
  await expect(page.locator('.picks__row')).toHaveCount(6);

  expect(browserErrors).toEqual([]);
});

test('picks can be set from typed coordinates alone', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await openPlate(page);
  await page.getByRole('button', { name: 'Auto-pick', exact: true }).click();
  await expect(page.locator('.picks__row').first()).toBeVisible();

  await page.locator('.picks__io summary').click();
  const area = page.locator('.picks__text');
  await area.fill('512, 512\n452.48, 709.30\n313.59, 452.73');
  await page.getByRole('button', { name: 'Apply these coordinates' }).click();

  await expect(page.locator('.picks__row')).toHaveCount(3);
  expect(await pickedCoordinates(page)).toEqual([
    [512, 512],
    [452.48, 709.3],
    [313.59, 452.73],
  ]);

  expect(browserErrors).toEqual([]);
});

test('a click on an uploaded micrograph snaps to the spot centroid', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await openPlate(page);

  // The micrograph is drawn in the page and handed to the file input as a real
  // File, so this exercises the upload path — including the pixel read that the
  // centroid needs — without a binary fixture in the repository.
  const truth = await page.evaluate(async () => {
    const size = 600;
    const centre = [300, 300];
    const basis = [
      [120, 20],
      [-20, 120],
    ];
    const canvas = document.createElement('canvas');
    canvas.width = size;
    canvas.height = size;
    const context = canvas.getContext('2d');
    context.fillStyle = '#080808';
    context.fillRect(0, 0, size, size);
    for (let m = -3; m <= 3; m += 1) {
      for (let n = -3; n <= 3; n += 1) {
        const x = centre[0] + m * basis[0][0] + n * basis[1][0];
        const y = centre[1] + m * basis[0][1] + n * basis[1][1];
        if (x < 0 || x >= size || y < 0 || y >= size) continue;
        const radius = m === 0 && n === 0 ? 16 : 10;
        const gradient = context.createRadialGradient(x, y, 0, x, y, radius);
        gradient.addColorStop(0, 'rgba(255,255,255,1)');
        gradient.addColorStop(1, 'rgba(255,255,255,0)');
        context.fillStyle = gradient;
        context.beginPath();
        context.arc(x, y, radius, 0, 2 * Math.PI);
        context.fill();
      }
    }
    const blob = await new Promise((resolve) => canvas.toBlob(resolve, 'image/png'));
    const transfer = new DataTransfer();
    transfer.items.add(new File([blob], 'plate.png', { type: 'image/png' }));
    const input = document.querySelector('input[type=file]');
    input.files = transfer.files;
    input.dispatchEvent(new Event('change', { bubbles: true }));
    return { centre, spot: [centre[0] + basis[0][0], centre[1] + basis[0][1]] };
  });
  await expect(page.locator(`${PATTERN_SVG} image`)).toBeVisible();

  const toScreen = async (x, y) =>
    page.evaluate(
      ([px, py]) => {
        const svg = document.querySelector('svg[aria-label="Diffraction pattern"]');
        const box = svg.viewBox.baseVal;
        const rect = svg.getBoundingClientRect();
        const scale = Math.min(rect.width / box.width, rect.height / box.height);
        return {
          x: rect.left + (rect.width - box.width * scale) / 2 + (px - box.x) * scale,
          y: rect.top + (rect.height - box.height * scale) / 2 + (py - box.y) * scale,
        };
      },
      [x, y],
    );

  // Deliberately four pixels off the spot, the way a click actually lands.
  const aimed = await toScreen(truth.spot[0] + 4, truth.spot[1] - 4);
  await page.mouse.click(aimed.x, aimed.y);
  await expect(page.locator('.picks__row')).toHaveCount(1);

  const [beam] = await pickedCoordinates(page);
  // One pass of centre-of-mass lands back near the click, because the window is
  // centred on the click rather than on the spot; the iteration is what makes
  // this assertion pass.
  expect(Math.hypot(beam[0] - truth.spot[0], beam[1] - truth.spot[1])).toBeLessThan(1.5);

  expect(browserErrors).toEqual([]);
});

/*
 * A pattern opened from disk is shown, whatever the frame was doing before.
 *
 * The camera is preserved across a redraw so that picking while zoomed in is
 * possible at all. Preserving it across an *open* is a different thing, and it
 * hid the file: at 448% on a 1024 px plate, a freshly opened 400 px micrograph
 * arrived as an 89 px crop of itself. The assertion is that opening refits.
 */
test('an opened micrograph is fitted rather than inheriting the last view', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await openPlate(page);

  await patternButton(page, 'Zoom in').click();
  await patternButton(page, 'Zoom in').click();
  await expect(patternControl(page, '.plot__zoom')).not.toHaveText('100%');

  const size = await page.evaluate(async () => {
    const side = 400;
    const canvas = document.createElement('canvas');
    canvas.width = side;
    canvas.height = side;
    const context = canvas.getContext('2d');
    context.fillStyle = '#0a0a0a';
    context.fillRect(0, 0, side, side);
    context.fillStyle = '#ffffff';
    context.beginPath();
    context.arc(side / 2, side / 2, 12, 0, 2 * Math.PI);
    context.fill();
    const blob = await new Promise((resolve) => canvas.toBlob(resolve, 'image/png'));
    const transfer = new DataTransfer();
    transfer.items.add(new File([blob], 'opened-plate.png', { type: 'image/png' }));
    const input = document.querySelector('input[type=file]');
    input.files = transfer.files;
    input.dispatchEvent(new Event('change', { bubbles: true }));
    return side;
  });

  await expect(page.locator(`${PATTERN_SVG} image`)).toBeVisible();
  await expect(patternControl(page, '.plot__zoom')).toHaveText('100%');
  // The whole micrograph, not a crop of it: the viewBox is the image itself.
  await expect
    .poll(async () =>
      page.locator(PATTERN_SVG).evaluate((node) => {
        const box = node.viewBox.baseVal;
        return [box.x, box.y, box.width, box.height].join(' ');
      }),
    )
    .toBe(`0 0 ${size} ${size}`);

  // And a redraw of that same pattern still keeps whatever view it is given,
  // because picking a spot precisely means picking it zoomed in.
  await patternButton(page, 'Zoom in').click();
  const zoomed = await patternControl(page, '.plot__zoom').textContent();
  // The element's centre, which is inside the image at any letterboxing: a
  // click in the margin beside a non-square drawing is refused as outside it.
  await page.locator(PATTERN_SVG).click();
  await expect(page.locator('.picks__row')).toHaveCount(1);
  await expect(patternControl(page, '.plot__zoom')).toHaveText(zoomed);

  expect(browserErrors).toEqual([]);
});

/*
 * The simulator: the forward pattern, its stated orientation, and its bands.
 *
 * Three claims, each of which would be invisible in a screenshot. The plate is
 * drawn with the shared drawing, so it looks like the solver's practice plates.
 * The orientation is *stated* beside it rather than implied — and when the
 * orientation asked for is not a zone axis, the deviation is on screen, because
 * a pattern drawn five degrees off [011] is a picture of [011]. And the Kikuchi
 * bands come from the same orientation the spots did, so the band for a plane
 * is as wide as its own spot is far from the beam, which is the relation the
 * overlay exists to teach.
 */
test('the SAED simulator states the orientation it drew and can band it', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await openPanel(page, 'SAED Simulator');

  const drawing = page.locator('#stage svg[aria-label="Simulated diffraction pattern"]');
  await expect(drawing).toBeVisible({ timeout: 20_000 });

  // The first example: zirconium down the basal axis, exactly on the zone.
  const readout = page.locator('#stage .plot__readout');
  await expect(readout).toContainText('Orientation on the beam');
  await expect(readout).toContainText('[0001]');
  await expect(readout).not.toContainText('Off that axis by');
  await expect(page.locator('#stage .plot__status')).toContainText('Zirconium');

  // Bands on request, fetched for the orientation the pattern was drawn from.
  const spotsBefore = await drawing.locator('circle').count();
  await page.locator('#stage').getByRole('button', { name: 'Kikuchi', exact: true }).click();
  await expect(page.locator('#stage .plot__status')).toContainText('Kikuchi band(s)', {
    timeout: 30_000,
  });
  expect(await drawing.locator('polyline').count()).toBeGreaterThan(0);
  // The bands are an addition to the pattern, not a replacement for it.
  expect(await drawing.locator('circle').count()).toBe(spotsBefore);

  // An orientation that is not a zone axis says how far off it is.
  await page.locator('#rail-body summary').filter({ hasText: 'Try an example' }).click();
  await page
    .locator('#rail-body .example')
    .filter({ hasText: 'A measured orientation' })
    .click();
  await expect(readout).toContainText('Off that axis by', { timeout: 30_000 });
  await expect(readout).toContainText('[011]');
  const deviation = await readout.locator('tr', { hasText: 'Off that axis by' }).locator('td').innerText();
  // Thirty, fifty, zero in Bunge is five degrees from [011]: an angle from the
  // Euler convention, not from this program.
  expect(Number(deviation.replace('°', ''))).toBeCloseTo(5.0, 1);

  // And the cursor reports the radius from 000 and the |g| it corresponds to.
  const cursorText = await page.evaluate(() => {
    const svg = document.querySelector('svg[aria-label="Simulated diffraction pattern"]');
    const box = svg.getBoundingClientRect();
    svg.dispatchEvent(
      new PointerEvent('pointermove', {
        bubbles: true,
        clientX: box.left + box.width * 0.72,
        clientY: box.top + box.height * 0.5,
      }),
    );
    return svg.closest('.plot').querySelector('.plot__cursor').textContent;
  });
  expect(cursorText).toContain('px from 000');
  expect(cursorText).toContain('|g|');

  expect(browserErrors).toEqual([]);
});

/*
 * The whole from-disk path, on a pattern whose answer is known.
 *
 * Everything else in this file starts from the gallery, where the plate arrives
 * as coordinates and is drawn as vectors. That never exercises what a user
 * actually does with their own data: choose a file, see it appear, set the
 * calibration for that exposure, click spots on pixels, and index. The tracked
 * fixture makes that testable — it is a simulated zirconium [0001] plate written
 * as a PNG, with a JSON sidecar carrying the answer, so the test can click where
 * the reflections are and require the zone axis back.
 */
test('a pattern file opened from disk is displayed, picked on, and indexed', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await openPlate(page);

  const truth = JSON.parse(readFileSync('fixtures/tem/zr_hcp_basal_saed.json', 'utf-8'));
  await page.setInputFiles('#rail-body input[type="file"]', 'fixtures/tem/zr_hcp_basal_saed.png');

  // It is on screen, whole, and it is the file that was opened.
  const drawn = page.locator(`${PATTERN_SVG} image`);
  await expect(drawn).toBeVisible();
  await expect(drawn).toHaveAttribute('width', String(truth.width_px));
  await expect(patternControl(page, '.plot__zoom')).toHaveText('100%');

  const toScreen = async (x, y) =>
    page.evaluate(
      ([px, py]) => {
        const svg = document.querySelector('svg[aria-label="Diffraction pattern"]');
        const box = svg.viewBox.baseVal;
        const rect = svg.getBoundingClientRect();
        const scale = Math.min(rect.width / box.width, rect.height / box.height);
        return {
          x: rect.left + (rect.width - box.width * scale) / 2 + (px - box.x) * scale,
          y: rect.top + (rect.height - box.height * scale) / 2 + (py - box.y) * scale,
        };
      },
      [x, y],
    );

  // The transmitted beam first, then three independent reflections — the picks
  // the sidecar names, which are the ones a user would choose for the same
  // reason: strong, and not collinear through the beam.
  const clicks = [truth.centre_px, ...truth.seed_spots.map((spot) => [spot.x, spot.y])];
  for (const [x, y] of clicks) {
    const point = await toScreen(x, y);
    await page.mouse.click(point.x, point.y);
  }
  await expect(page.locator('.picks__row')).toHaveCount(clicks.length);

  // The measurements are readable beside the picture rather than under it, and
  // the ratio of the two inner rings is the sqrt(3) of the hexagonal net.
  const readout = patternControl(page, '.plot__readout');
  await expect(readout).toContainText('Measured picks');

  // This exposure's phase and calibration, neither of which the practice plate's
  // are: the fields still hold the gallery's aluminium and its camera constant,
  // and indexing on those returns no solution at all — which is the honest
  // failure, and the reason the panel warns when a file is opened over a plate.
  // Three forms carry a phase picker — index, atlas, tilt — so this names the
  // one belonging to the step that indexes.
  const indexStep = page.locator('#rail-body details.step').filter({
    has: page.getByRole('button', { name: 'Index the pattern', exact: true }),
  });
  await indexStep.locator('[id^="ctl-phase-"]').selectOption('zr_hcp');
  await page
    .locator('[id^="ctl-camera_constant_mm_angstrom-"]')
    .fill(String(truth.camera_constant_mm_angstrom));
  await page.locator('[id^="ctl-pixel_size_mm-"]').fill(String(truth.pixel_size_mm));

  await page.getByRole('button', { name: 'Index the pattern', exact: true }).click();
  // Zirconium down [0001], which is what the file was built from — the answer
  // the sidecar states, reached by opening a picture and clicking on it.
  await expect(page.locator('#stage')).toContainText('Zirconium (hcp, alpha) down [0001]', {
    timeout: 30_000,
  });
  await expect(page.locator('#stage')).toContainText('3 of 3 picked spots were indexed');

  expect(browserErrors).toEqual([]);
});

/*
 * Measured pole figures, in tabs, on one scale.
 *
 * The reason to measure more than one pole figure is to compare them, and two
 * figures on two scales cannot be compared. So the assertions are that both
 * files arrive, that each gets a tab, that the tab switches the drawing, and
 * that the contour levels typed into the form are the ones the figure reports.
 */
test('XRDML pole figures open into tabs on one shared scale', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await workspaceTab(page, 'Texture').click();

  const text = readFileSync('fixtures/xrdml/synthetic_random_standard.xrdml', 'utf-8');
  await page.locator('#rail-body select[aria-label="View"]').selectOption({
    label: 'Measured pole figures',
  });
  await expect(page.locator('#stage')).toContainText('Open one or more XRDML');

  await page.setInputFiles('#texture-files input[type="file"]', [
    { name: 'ni-111.xrdml', mimeType: 'application/xml', buffer: Buffer.from(text, 'utf-8') },
    { name: 'ni-200.xrdml', mimeType: 'application/xml', buffer: Buffer.from(text, 'utf-8') },
  ]);
  await expect(page.locator('#rail-body')).toContainText('2 file(s) open');

  // Two planes, one per file, then redraw. Each index goes in its own box.
  const poles = page.locator('.indices--multi').filter({ has: page.locator('[id^="ctl-poles-"]') });
  await fillIndices(poles, [
    [1, 1, 1],
    [2, 0, 0],
  ]);
  await page.locator('[id^="ctl-contour_levels-"]').fill('0.8, 1, 1.2, 1.4');
  await page.getByRole('button', { name: 'Build texture', exact: true }).click();

  const tabs = page.locator('.figure-tab');
  await expect(tabs).toHaveCount(2, { timeout: 30_000 });
  await expect(tabs.nth(0)).toHaveText('{111}');
  await expect(tabs.nth(1)).toHaveText('{200}');
  await expect(tabs.nth(0)).toHaveAttribute('aria-selected', 'true');

  const status = page.locator('#stage .plot__status').first();
  await expect(status).toContainText('{111} from ni-111.xrdml');
  await expect(status).toContainText('one scale across every figure');
  await expect(status).toContainText('contours at 0.80, 1.00, 1.20, 1.40');

  await tabs.nth(1).click();
  await expect(status).toContainText('{200} from ni-200.xrdml');
  // The shared scale is the point: the second figure reports the same range.
  await expect(tabs.nth(1)).toHaveAttribute('aria-selected', 'true');

  // The reconstructed ODF is a further tab, after the figures it came from.
  await page.locator('[id^="ctl-reconstruct_odf-"]').check();
  await page.getByRole('button', { name: 'Build texture', exact: true }).click();
  await expect(tabs).toHaveCount(3, { timeout: 60_000 });
  await expect(tabs.nth(2)).toHaveText('ODF');
  await tabs.nth(2).click();
  await expect(status).toContainText('ODF from 2 measured figure(s)');
  // The inversion is ill-posed, and the figure has to say so where it is read.
  await expect(status).toContainText('residual');
  await expect(status).toContainText('ill-posed');

  expect(browserErrors).toEqual([]);
});

/*
 * Opening a real EBSD scan.
 *
 * The panel's practice datasets are constructions with known answers, which is
 * what makes them checkable; the point of this path is that a *file* reaches
 * exactly the same analysis. The scan below is four points, three of one
 * orientation and one 60 degrees away, so the answer — two grains — is worked
 * out by hand rather than read off the program.
 */
test('an EBSD scan file is opened and analysed like a practice map', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await workspaceTab(page, 'EBSD').click();
  await expect(page.locator('#stage svg')).toBeVisible({ timeout: 30_000 });

  const scan = [
    '# Phase 1',
    '# MaterialName  \tNickel',
    '# Symmetry              43',
    '# LatticeConstants      3.520 3.520 3.520  90.000  90.000  90.000',
    '# GRID: SqrGrid',
    '# XSTEP: 1.000000',
    '# YSTEP: 1.000000',
    '# NCOLS_ODD: 2',
    '# NCOLS_EVEN: 2',
    '# NROWS: 2',
    '#',
    '   0.00000   0.00000   0.00000      0.00000      0.00000  60.0  0.950  0  1  0.500',
    '   0.00000   0.00000   0.00000      1.00000      0.00000  55.0  0.900  0  1  0.600',
    '   0.00000   0.00000   0.00000      0.00000      1.00000  50.0  0.850  0  1  0.700',
    '   1.04720   0.00000   0.00000      1.00000      1.00000  45.0  0.800  0  1  0.800',
    '',
  ].join('\n');

  await page.setInputFiles('#rail-body input[type="file"]', {
    name: 'bicrystal.ang',
    mimeType: 'text/plain',
    buffer: Buffer.from(scan, 'utf-8'),
  });

  await expect(page.locator('#rail-body')).toContainText('bicrystal.ang');
  // The figure is of the file, and the analysis is the real one.
  await expect(page.locator('#stage')).toContainText('bicrystal.ang', { timeout: 30_000 });
  await expect(page.locator('#stage')).toContainText('2 grains');
  // And it refuses to claim the guarantee a constructed dataset carries.
  await expect(page.locator('#stage')).toContainText('measurement, not a construction');

  // Closing it goes back to the practice dataset rather than to an empty stage.
  await page.getByRole('button', { name: 'Close the scan', exact: true }).click();
  await expect(page.locator('#stage')).not.toContainText('bicrystal.ang', { timeout: 30_000 });

  expect(browserErrors).toEqual([]);
});

/*
 * Opening an HDF5 scan.
 *
 * An EDAX OIM scan is `.oh5` or `.h5` — one HDF5 container under two
 * extensions — so it cannot ride the request as text the way a `.ang` does; the
 * browser base64-encodes it into the same field instead. Building a real HDF5
 * file here would mean committing a binary blob nobody can review, so what is
 * checked is the part that is this file's to check: that the bytes survive the
 * trip. A file that is *not* HDF5, sent under an `.oh5` name, must come back
 * with the reader's own complaint about its contents — which it can only do if
 * the bytes were encoded, decoded, written to a temporary file and handed to
 * `read_oh5`. Silence, or a complaint about the extension, would mean the
 * transport broke before the reader ever saw it.
 */
test('an HDF5 scan travels as bytes and is read as one', async ({ page }) => {
  await openWorkbench(page);
  await workspaceTab(page, 'EBSD').click();
  await expect(page.locator('#stage svg')).toBeVisible({ timeout: 30_000 });

  await page.setInputFiles('#rail-body input[type="file"]', {
    name: 'not-really.oh5',
    mimeType: 'application/x-hdf5',
    buffer: Buffer.from([0x00, 0x01, 0x02, 0xff, 0xfe, 0x89, 0x48, 0x44]),
  });

  await expect(page.locator('#rail-body')).toContainText('not-really.oh5', { timeout: 30_000 });
  await expect(page.locator('#rail-body')).toContainText('could not be read as a .oh5 scan', {
    timeout: 30_000,
  });
});

/*
 * Calibrating from the image itself.
 *
 * An image that arrives by email has a scale bar and no recorded camera length.
 * Drawing a line across a length that *is* known measures the one number the
 * camera equation uses. The test does the whole loop on a practice plate whose
 * answer is known: measure beam-to-(200), say what that reflection's |g| is,
 * index, and require the zone axis back. Anything less than the full loop would
 * check that a number reached a field, not that the number is right.
 */
test('a line of known length calibrates the pattern well enough to index it', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  const toScreen = await openPlate(page);

  await patternButton(page, 'Calibrate').click();
  // The beam and the (200) reflection of the aluminium plate. Aluminium's
  // a = 4.0495 A, so d(200) = 2.02475 A and |g| = 0.493886 1/A — a number from
  // the lattice parameter, not from this program.
  for (const [x, y] of [
    [512, 512],
    [452.48, 709.3],
  ]) {
    const point = await toScreen(x, y);
    await page.mouse.click(point.x, point.y);
  }
  const strip = page.locator('.calibrate');
  await expect(strip).toContainText('px');
  await strip.locator('input').fill(String(1 / 2.02475));
  await strip.locator('select').selectOption('inv_angstrom');
  await strip.getByRole('button', { name: 'Use this scale', exact: true }).click();

  // The measured scale, in the units the form now works in.
  await expect(page.locator('#stage .plot__status').first()).toContainText('1 px =');
  const units = page.locator('#rail-body select').filter({ hasText: 'measured scale' });
  await expect(units).toHaveValue('px_scale');
  const scale = Number(
    await page.locator('[id^="ctl-reciprocal_per_px_angstrom-"]').inputValue(),
  );
  // The plate was built at 0.024 mm per pixel and a camera constant of
  // 10.0317 mm.A, so the true scale is 0.0023924 1/A per pixel. Measuring one
  // reflection recovers it to better than a percent.
  expect(Math.abs(scale - 0.0023924) / 0.0023924).toBeLessThan(0.01);

  // And it is a calibration, not a display setting: the pattern indexes on it.
  await page.getByRole('button', { name: 'Auto-pick', exact: true }).click();
  await page.getByRole('button', { name: 'Index the pattern', exact: true }).click();
  await expect(page.locator('#stage')).toContainText('[001]', { timeout: 30_000 });

  expect(browserErrors).toEqual([]);
});

/*
 * The stereogram beside the pattern.
 *
 * Two things are worth pinning. That it is *beside* the pattern rather than
 * under it — the geometry that recovered half a laptop window. And that the
 * cursor readout is right: the panel inlines the stage closed form so it can
 * answer while the pointer moves, and an inlined formula that drifts from the
 * server's is exactly the kind of defect nobody notices, because it stays
 * plausible. So the readout is compared against the server's own numbers at the
 * poles the server placed.
 */
test('the stereogram sits beside the pattern and reads out the tilt under the cursor', async ({
  page,
}) => {
  await page.setViewportSize({ width: 1600, height: 900 });
  const browserErrors = await openWorkbench(page);
  await openPanel(page, 'TEM Solver');

  const stereogram = page.locator('#stage svg[aria-label="Stereogram"]');
  await expect(stereogram).toBeVisible({ timeout: 20_000 });

  // Side by side, each taking about half the stage, not stacked.
  const geometry = await page.evaluate(() => {
    const figures = [...document.querySelectorAll('.tem-stage > .plot')];
    return figures.map((figure) => {
      const rect = figure.getBoundingClientRect();
      return { x: rect.x, y: rect.y, width: rect.width };
    });
  });
  expect(geometry).toHaveLength(2);
  expect(geometry[1].x).toBeGreaterThan(geometry[0].x + geometry[0].width - 2);
  expect(Math.abs(geometry[0].y - geometry[1].y)).toBeLessThan(2);

  // The poles the server placed, straight from the service.
  const axes = await page.evaluate(async () => {
    const response = await fetch('/api/call', {
      method: 'POST',
      headers: { 'content-type': 'application/json' },
      body: JSON.stringify({
        operation: 'tem.stereogram',
        params: { phase: { builtin: 'al_fcc' }, zone_axis: [0, 0, 1] },
      }),
    });
    const body = await response.json();
    return body.result.data.axes;
  });
  const wanted = ['[001]', '[011]', '[111]', '[112]'];
  const poles = wanted.map((label) => {
    const found = axes.find((entry) => entry.label === label);
    expect(found, `the service placed no pole ${label}`).toBeTruthy();
    return found;
  });

  for (const pole of poles) {
    const readout = await page.evaluate((entry) => {
      const svg = document.querySelector('#stage svg[aria-label="Stereogram"]');
      const box = svg.viewBox.baseVal;
      const rect = svg.getBoundingClientRect();
      const scale = Math.min(rect.width / box.width, rect.height / box.height);
      const offsetX = (rect.width - box.width * scale) / 2;
      const offsetY = (rect.height - box.height * scale) / 2;
      svg.dispatchEvent(
        new PointerEvent('pointermove', {
          bubbles: true,
          clientX: rect.left + offsetX + (entry.x - box.x) * scale,
          // The drawing's y runs up; the SVG's runs down.
          clientY: rect.top + offsetY + (-entry.y - box.y) * scale,
        }),
      );
      return svg.closest('.plot').querySelector('.plot__cursor').textContent;
    }, pole);

    const numbers = [...readout.matchAll(/-?\d+\.\d+/g)].map(Number);
    expect(numbers.length, `readout without numbers: ${readout}`).toBeGreaterThanOrEqual(3);
    // Within a tenth of a degree, which is what the readout prints.
    expect(Math.abs(numbers[0] - pole.alpha_deg)).toBeLessThan(0.11);
    expect(Math.abs(numbers[1] - pole.beta_deg)).toBeLessThan(0.11);
    expect(Math.abs(numbers[2] - pole.angle_from_beam_deg)).toBeLessThan(0.11);
    expect(readout).toContain(pole.label);
  }

  expect(browserErrors).toEqual([]);
});

/*
 * The Kikuchi bands appear with the indexing, which is what supplies an
 * orientation to draw them from.
 *
 * Three things are worth pinning in a browser rather than in Python. The toggle
 * must not exist before the pattern is indexed, because before that there is no
 * orientation and a control that draws nothing reads as a broken one; it must
 * appear as soon as there is one, without waiting for an acceptance, since
 * every candidate carries its own. And the bands must go inside the clipped
 * group with everything else: a band centre line is generated from a lattice
 * plane and runs the whole width of the pattern, so an overlay appended to the
 * outer element would paint across the margins beside the picture.
 */
test('the indexed solution draws its Kikuchi bands inside the clip', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await openPlate(page);

  const kikuchi = patternButton(page, 'Kikuchi');
  await expect(kikuchi).toBeHidden();

  await page.getByRole('button', { name: 'Auto-pick', exact: true }).click();
  await expect(page.locator('.picks__row')).toHaveCount(7, { timeout: 20_000 });
  // Picks alone are not an orientation, so there is still nothing to draw.
  await expect(kikuchi).toBeHidden();

  await page.getByRole('button', { name: 'Index the pattern', exact: true }).click();

  // It arrives with the indexing, before any acceptance.
  await expect(kikuchi).toBeVisible({ timeout: 20_000 });
  await expect(kikuchi).toHaveAttribute('aria-pressed', 'false');
  await kikuchi.click();
  await expect(kikuchi).toHaveAttribute('aria-pressed', 'true');

  const drawn = await expect
    .poll(
      async () =>
        page.evaluate(() => {
          const svg = document.querySelector('svg[aria-label="Diffraction pattern"]');
          return svg.querySelectorAll('polyline').length;
        }),
      { timeout: 20_000 },
    )
    .toBeGreaterThan(0);
  void drawn;

  const geometry = await page.evaluate(() => {
    const svg = document.querySelector('svg[aria-label="Diffraction pattern"]');
    const labels = [...svg.querySelectorAll('text')].map((node) => node.textContent);
    const box = svg.viewBox.baseVal;
    let farthest = 0;
    for (const polyline of svg.querySelectorAll('polyline')) {
      for (const pair of polyline.getAttribute('points').split(' ')) {
        const [x, y] = pair.split(',').map(Number);
        farthest = Math.max(farthest, Math.abs(x - box.width / 2), Math.abs(y - box.height / 2));
      }
    }
    return {
      labels,
      farthest,
      halfWidth: box.width / 2,
      unclipped: [...svg.children]
        .filter((child) => child.tagName !== 'defs')
        .filter((child) => !child.getAttribute('clip-path'))
        .map((child) => child.tagName),
    };
  });

  // Nothing new escaped the clip, and the band edges really do run past the
  // frame — so the case cannot pass by the overlay quietly drawing nothing.
  expect(geometry.unclipped).toEqual([]);
  expect(geometry.farthest).toBeGreaterThan(geometry.halfWidth);
  // Bands are named as planes, and the route to the target names the band to
  // follow rather than a tilt nobody can dial without the holder calibration.
  expect(geometry.labels.some((label) => /^\(\d/.test(label ?? ''))).toBe(true);
  expect(geometry.labels.some((label) => /^follow \(.*\) toward \[/.test(label ?? ''))).toBe(true);

  /*
   * Every band's name runs along that band.
   *
   * On a zone-axis plate a dozen bands cross within a few tens of pixels, and a
   * horizontal caption beside a steeply running band belongs to whichever line
   * happens to be nearest the text. The pairing here is by drawing order rather
   * than by proximity — each band appends its centre line, then its edges, then
   * its name — so the assertion is about the band the label was drawn for, not
   * about the band that ended up closest to it.
   */
  const slopes = await page.evaluate(() => {
    const svg = document.querySelector('svg[aria-label="Diffraction pattern"]');
    const fold = (angle) => (angle > 90 ? angle - 180 : angle <= -90 ? angle + 180 : angle);
    const pairs = [];
    let line = null;
    for (const node of svg.querySelectorAll('line, text')) {
      if (node.tagName === 'line') {
        line = node;
        continue;
      }
      const transform = node.getAttribute('transform') ?? '';
      const match = /rotate\(([-\d.]+)/.exec(transform);
      if (!match || !line) continue;
      const dx = Number(line.getAttribute('x2')) - Number(line.getAttribute('x1'));
      const dy = Number(line.getAttribute('y2')) - Number(line.getAttribute('y1'));
      pairs.push({
        label: Number(match[1]),
        band: fold((Math.atan2(dy, dx) * 180) / Math.PI),
      });
    }
    return pairs;
  });

  expect(slopes.length).toBeGreaterThan(0);
  for (const { label, band } of slopes) expect(Math.abs(label - band)).toBeLessThan(0.6);
  // And the plate really does carry sloping bands, so the case cannot pass by
  // every band happening to be horizontal.
  expect(slopes.some(({ band }) => Math.abs(band) > 5)).toBe(true);

  // The status line says what the overlay is and what it is not.
  const status = await patternControl(page, '.plot__status').textContent();
  expect(status).toContain('000');
  expect(status).toMatch(/not modelled/);

  // Accepting a solution settles which orientation the bands come from; it
  // neither removes the toggle nor turns it off.
  await page
    .getByRole('button', { name: 'Accept this solution', exact: true })
    .click({ timeout: 20_000 });
  await expect(kikuchi).toBeVisible();
  await expect(kikuchi).toHaveAttribute('aria-pressed', 'true');
  await expect
    .poll(() =>
      page.evaluate(() => {
        const svg = document.querySelector('svg[aria-label="Diffraction pattern"]');
        return svg.querySelectorAll('polyline').length;
      }),
    )
    .toBeGreaterThan(0);

  // Turning it off takes the bands away rather than leaving them stale.
  await kikuchi.click();
  await expect
    .poll(() =>
      page.evaluate(() => {
        const svg = document.querySelector('svg[aria-label="Diffraction pattern"]');
        return svg.querySelectorAll('polyline').length;
      }),
    )
    .toBe(0);

  expect(browserErrors).toEqual([]);
});

/*
 * About is a legal surface as much as a courtesy one: the GPL asks an
 * interactive program to display its warranty disclaimer, and a user filing a
 * bug needs the build number. The assertion is therefore on the content, not
 * merely on the drawer opening.
 */
/*
 * The third figure in the dock: the crystal's own Kikuchi map.
 *
 * Unlike the two beside it, this one does *not* turn with the camera, and that
 * is the claim being tested. The pole figure and the inverse pole figure are the
 * same view of the same crystal as the structure; the map is the atlas — fixed
 * to the crystal, centred where the user says — and what moves across it is the
 * marker showing which direction the current view has on the beam. Both halves
 * of that are asserted, because a map that quietly rotated with the camera would
 * look plausible and mean something else entirely.
 */
test('the crystal viewer maps the Kikuchi bands about an axis the user chooses', async ({
  page,
}) => {
  const browserErrors = await openWorkbench(page);
  const map = page.locator('svg[aria-label="Kikuchi map of the crystal"]');
  await expect(map).toBeVisible({ timeout: 30_000 });
  // Bands are drawn as their two edges with the plane trace between them.
  await expect.poll(() => map.locator('polyline').count()).toBeGreaterThan(10);

  const status = page.locator('.orient__cell').filter({ hasText: 'Kikuchi map' }).locator('.field__hint');
  await expect(status).toContainText('zone axes within 60° of [001]');
  // Nickel down [001]: the axes a standard cubic projection carries.
  await expect(map).toContainText('[001]');
  await expect(map).toContainText('[011]');

  // Re-centring is the user's, and it re-computes rather than re-drawing.
  await page.locator('input.orient__axis').fill('1 1 1');
  await page.locator('input.orient__axis').press('Enter');
  await expect(status).toContainText('of [111]', { timeout: 30_000 });
  await expect(map).toContainText('[111]');

  // Something that is not a direction is refused in place, without a toast and
  // without taking the map away.
  await page.locator('input.orient__axis').fill('oops');
  await page.locator('input.orient__axis').press('Enter');
  await expect(status).toContainText('Three whole numbers');
  await expect(map.locator('polyline').first()).toBeVisible();

  expect(browserErrors).toEqual([]);
});

/*
 * The map is examined, not merely displayed.
 *
 * A Kikuchi band at 200 kV is a fraction of a degree wide on a map spanning
 * sixty, so a whole-hemisphere drawing shows the *network* and not the indices:
 * the names of the bands are unreadable at the size the figure occupies in a
 * dock. Magnifying it is therefore not a convenience, it is what makes the
 * figure an atlas — and the names have to be written along their own bands, for
 * the same reason they are on the plate.
 */
test('the crystal viewer magnifies its Kikuchi map and names the bands along them', async ({
  page,
}) => {
  const browserErrors = await openWorkbench(page);
  const map = page.locator('svg[aria-label="Kikuchi map of the crystal"]');
  await expect(map).toBeVisible({ timeout: 30_000 });
  await expect.poll(() => map.locator('polyline').count()).toBeGreaterThan(10);

  // A plane is named as a plane, and its name is turned to the slope of the
  // band it belongs to. Nickel down [001] has bands at every slope, so a
  // network whose labels were all upright would fail this.
  const bandAngles = () =>
    page.evaluate(() => {
      const svg = document.querySelector('svg[aria-label="Kikuchi map of the crystal"]');
      return [...svg.querySelectorAll('text')]
        .filter((node) => /^\(/.test(node.textContent ?? ''))
        .map((node) =>
          Number(/rotate\(([-\d.]+)/.exec(node.getAttribute('transform') ?? '')?.[1]),
        );
    });
  await expect.poll(async () => (await bandAngles()).length).toBeGreaterThan(0);
  const named = await bandAngles();
  expect(named.every((angle) => Number.isFinite(angle))).toBe(true);
  expect(named.some((angle) => Math.abs(angle) > 5)).toBe(true);

  const viewWindow = () =>
    map.evaluate((node) => {
      const box = node.viewBox.baseVal;
      return { width: box.width, x: box.x, y: box.y };
    });
  const fitted = await viewWindow();

  // The wheel magnifies about the pointer: the window narrows, and it does not
  // stay centred on the middle of the map when the pointer is not there.
  const figure = page.locator('.orient__canvas--zoomable');
  const box = await figure.boundingBox();
  await page.mouse.move(box.x + box.width * 0.3, box.y + box.height * 0.35);
  for (let notch = 0; notch < 8; notch += 1) await page.mouse.wheel(0, -120);
  await expect.poll(async () => (await viewWindow()).width).toBeLessThan(fitted.width * 0.9);
  const zoomed = await viewWindow();
  expect(Math.hypot(zoomed.x - fitted.x, zoomed.y - fitted.y)).toBeGreaterThan(0);

  // Magnification buys room, and room buys names: no band loses its name for
  // being looked at more closely.
  expect((await bandAngles()).length).toBeGreaterThanOrEqual(named.length);

  // Double-click fits the whole map again, whatever state the view was left in.
  await figure.dblclick();
  await expect.poll(async () => (await viewWindow()).width).toBeCloseTo(fitted.width, 3);

  expect(browserErrors).toEqual([]);
});

/*
 * The pole figure's fly-by, and the Bunge angles of whatever is on screen.
 *
 * The inverse pole figure has had a trail since the dock was built; the pole
 * figure had none, so the one figure that shows *where the poles go* showed no
 * history of them going anywhere. And the angle fields follow the convention
 * picker, which means a reader who switched to Matthies had no Bunge triple in
 * front of them — the convention every EBSD file and every published orientation
 * is written in.
 */
test('the pole figure leaves a fly-by trail and the view states its Bunge angles', async ({
  page,
}) => {
  const browserErrors = await openWorkbench(page);
  const poleFigure = page.locator('svg[aria-label="Pole figure of the current view"]');
  await expect(poleFigure).toBeVisible({ timeout: 30_000 });

  // One moment is recorded as soon as the figure is drawn — the trail starts
  // where the crystal is — so what a fly-by adds is measured against that.
  const trailDots = poleFigure.locator('circle[r="1.5"]');
  const atRest = await trailDots.count();

  // Turn the crystal by dragging the structure, which is what a fly-by is.
  const scene = page.locator('#stage svg').first();
  const box = await scene.boundingBox();
  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  await page.mouse.down();
  for (let step = 1; step <= 12; step += 1) {
    await page.mouse.move(
      box.x + box.width / 2 + step * 8,
      box.y + box.height / 2 + step * 4,
      { steps: 2 },
    );
  }
  await page.mouse.up();

  // The poles left a trail behind them.
  await expect.poll(() => trailDots.count()).toBeGreaterThan(atRest);

  // And the readout states Bunge angles, whichever convention the picker holds.
  const readout = page.locator('.orient__bunge');
  await expect(readout).toContainText('Bunge (φ₁, Φ, φ₂)', { timeout: 30_000 });
  const inBunge = await readout.textContent();

  await page.locator('#rail-body select[aria-label="Euler-angle convention"]').selectOption('matthies');
  await expect(readout).toContainText('Bunge (φ₁, Φ, φ₂)');
  // The same orientation, so the same Bunge triple: switching how the angles are
  // *entered* must not change what the view is.
  await expect.poll(async () => readout.textContent()).toBe(inBunge);

  expect(browserErrors).toEqual([]);
});

/*
 * The EBSD workspace: six views of one scan, then two forward tools.
 *
 * The sub-tabs are not six panels with six copies of the same logic — the three
 * map tabs are one panel opened on three colourings — and they are not six
 * independent sessions either. The scan a user opens belongs to the workspace,
 * so opening a file on the summary and then going to the map must analyse *that
 * file*. Silently reverting to the practice dataset next to somebody's own data
 * is the worst answer available, and it is the one this asserts against.
 *
 * The seventh and eighth tabs are the odd ones out by design: the Kikuchi
 * simulator and ECCI workflow take no scan, because they start from a phase and
 * a stated orientation rather than from a measured map.
 */
test('the EBSD workspace shows six scan views and two forward tools', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await workspaceTab(page, 'EBSD').click();

  const subtabs = page.locator('#subtabs').getByRole('tab');
  await expect(subtabs).toHaveText([
    'IPF map',
    'GROD',
    'KAM',
    'Scan summary',
    'Distributions',
    'Pole figures',
    'Kikuchi simulator',
    'ECCI workflow',
  ]);

  // The three map tabs are the same panel opened on different colourings, and
  // each says which map it drew.
  for (const [name, title] of [
    ['GROD', 'GROD map'],
    ['KAM', 'KAM map'],
    ['IPF map', 'IPF-Z map'],
  ]) {
    await openPanel(page, name);
    await expect(page.locator('#stage')).toContainText(title, { timeout: 30_000 });
  }

  // The summary is a page of numbers, sectioned, with its headline figures out
  // in front and each of them stating what it depends on.
  await openPanel(page, 'Scan summary');
  await expect(page.locator('.summary-figure').first()).toBeVisible({ timeout: 30_000 });
  await expect(page.locator('.summary-section .card__title')).toHaveText([
    'Acquisition',
    'Indexing quality',
    'Phases',
    'Microstructure',
  ]);
  await expect(page.locator('.summary-headline')).toContainText('threshold');

  // The distribution is a histogram with its statistics under it.
  await openPanel(page, 'Distributions');
  const bars = page.locator('#stage svg rect');
  await expect.poll(() => bars.count(), { timeout: 30_000 }).toBeGreaterThan(5);
  await expect(page.locator('#stage .plot__status')).toContainText('median');

  // The discrete figure is the measurement rather than a contour of it, and
  // says so along with how much of the scan it drew.
  await openPanel(page, 'Pole figures');
  await expect(page.locator('#stage .plot__status')).toContainText('not a density estimate', {
    timeout: 30_000,
  });

  expect(browserErrors).toEqual([]);
});

test('a scan opened in one EBSD view is the scan every view analyses', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await openPanel(page, 'Scan summary');
  await expect(page.locator('.summary-figure').first()).toBeVisible({ timeout: 30_000 });

  // Opened on the summary, not on the map.
  await page.setInputFiles('#rail-body input[type="file"]', 'fixtures/ebsd/synthetic_hex_grid.ang');
  await expect(page.locator('#rail-body')).toContainText('synthetic_hex_grid.ang', {
    timeout: 30_000,
  });
  // The scan is hexagonal and tiny, unlike every practice dataset.
  await expect(page.locator('.summary-cards')).toContainText('hexagonal', { timeout: 30_000 });

  // And the map, reached afterwards, is a map of that scan rather than of the
  // practice dataset its own control still names.
  await openPanel(page, 'IPF map');
  await expect(page.locator('#rail-body')).toContainText('synthetic_hex_grid.ang');
  await expect(page.locator('#stage')).toContainText('synthetic_hex_grid.ang', {
    timeout: 30_000,
  });

  // Closing it anywhere closes it everywhere.
  await page.getByRole('button', { name: 'Close the scan', exact: true }).click();
  await expect(page.locator('#rail-body')).toContainText('No scan open', { timeout: 30_000 });
  await openPanel(page, 'Distributions');
  await expect(page.locator('#rail-body')).toContainText('No scan open');

  expect(browserErrors).toEqual([]);
});

test('the About panel states the build, the author and the licence', async ({ page }) => {
  const browserErrors = await openWorkbench(page);

  const version = await page.evaluate(async () => {
    const response = await fetch('/api/manifest');
    return (await response.json()).version;
  });

  await page.getByRole('button', { name: 'About' }).click();
  const drawer = page.locator('#about-drawer');
  await expect(drawer).toBeVisible();
  await expect(drawer).toContainText(`PyTex ${version}`);
  await expect(drawer).toContainText('Dr K V Mani Krishna');
  await expect(drawer).toContainText('Materials Group');
  await expect(drawer).toContainText('WITHOUT ANY WARRANTY');
  await expect(drawer.getByRole('link', { name: 'kvmani@barc.gov.in' })).toHaveAttribute(
    'href',
    'mailto:kvmani@barc.gov.in',
  );
  await expect(drawer.getByRole('link', { name: 'kvmani@gmail.com' })).toBeVisible();

  await page.keyboard.press('Escape');
  await expect(drawer).toBeHidden();
  expect(browserErrors).toEqual([]);
});

test('keeps all workspaces reachable in the narrow responsive layout', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  const browserErrors = await openWorkbench(page);

  for (const workspace of WORKSPACES) {
    await expect(workspaceTab(page, workspace)).toBeVisible();
  }
  await expect(page.locator('.masthead__action-label').first()).toBeHidden();
  const overflow = await page.evaluate(() => ({
    clientWidth: document.documentElement.clientWidth,
    scrollWidth: document.documentElement.scrollWidth,
  }));
  expect(overflow.scrollWidth).toBeLessThanOrEqual(overflow.clientWidth + 1);

  await workspaceTab(page, 'Calculator').click();
  await expect(page.getByRole('button', { name: 'Calculate', exact: true })).toBeVisible();
  expect(browserErrors).toEqual([]);
});

/*
 * A long table is a preview; the export is the data.
 *
 * Both halves matter and neither is visible from the other. The card must say
 * that it is showing a subset — a table that quietly truncates is worse than one
 * that is too long, because a reader counting rows gets the wrong answer and
 * never knows — and the export must be built from the whole result rather than
 * from the rows on screen. The second half is checked the only way that settles
 * it: by taking the download and counting its lines.
 */
test('the table on screen is capped and the export carries every row', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await workspaceTab(page, 'Texture').click();

  const subtitle = page.locator('.card__subtitle');
  await expect(subtitle).toContainText('Showing the first 200 of', { timeout: 30_000 });
  const caption = await subtitle.textContent();
  const total = Number(/first 200 of (\d+) rows/.exec(caption)[1]);
  expect(total).toBeGreaterThan(200);
  await expect(page.locator('.table-wrap tbody tr')).toHaveCount(200);

  const download = page.waitForEvent('download');
  await page.locator('.card__header').getByRole('button', { name: 'CSV', exact: true }).click();
  const file = await download;
  const text = readFileSync(await file.path(), 'utf-8');
  const lines = text.trim().split(/\r?\n/);
  // One header line, then every row of the result — not the two hundred drawn.
  expect(lines.length).toBe(total + 1);

  expect(browserErrors).toEqual([]);
});

/*
 * A legend button must survive being pressed.
 *
 * Both plotting panels have a legend whose buttons toggle what is drawn, and
 * both once rebuilt the whole legend as part of the redraw the click causes.
 * The button the user just pressed is then removed from the document and the
 * browser moves focus to `body`, so a keyboard user who tabs to a packet and
 * presses Enter loses their place and has to tab through the entire page to
 * reach the next one. Focus is the observable, so focus is what is asserted.
 */
test('toggling a legend entry keeps the focus on the entry', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await workspaceTab(page, 'Variants').click();

  const entry = page.locator('.legend__item').first();
  await expect(entry).toBeVisible({ timeout: 30_000 });
  const label = await entry.textContent();
  await entry.focus();
  await page.keyboard.press('Enter');

  const focused = await page.evaluate(() => ({
    tag: document.activeElement.tagName,
    text: document.activeElement.textContent,
    inLegend: Boolean(document.activeElement.closest('.legend')),
  }));
  expect(focused.tag).toBe('BUTTON');
  expect(focused.inLegend).toBe(true);
  expect(focused.text).toBe(label);
  // And the press did something, so the case cannot pass on a dead button.
  await expect(entry).toHaveAttribute('aria-pressed', 'false');

  expect(browserErrors).toEqual([]);
});

/*
 * The colour theme is the page's, not the launcher's.
 *
 * The same page runs in a browser tab and in the desktop shell's web view, so a
 * theme control in either launcher would exist in one and not the other. What
 * has to hold is that the choice is explicit, that it survives a reload, and
 * that `auto` is a real third state — a binary switch would leave a user who
 * follows their operating system with no way back to it.
 */
test('the colour theme is chosen in the page and survives a reload', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  const button = page.locator('#cycle-theme');
  const theme = () => page.evaluate(() => ({
    attribute: document.documentElement.getAttribute('data-theme'),
    stored: localStorage.getItem('pytex-theme'),
  }));

  const seen = [];
  for (let press = 0; press < 3; press += 1) {
    await button.click();
    const state = await theme();
    expect(state.attribute).toBe(state.stored);
    seen.push(state.attribute);
  }
  // Three states, cycled: light, dark, and back to following the system.
  expect(new Set(seen)).toEqual(new Set(['light', 'dark', 'auto']));

  await button.click();
  const before = await theme();
  await page.reload();
  await expect(page.locator('#tabs')).toBeVisible();
  expect(await theme()).toEqual(before);

  expect(browserErrors).toEqual([]);
});

/*
 * The crystal viewer's orientation dock.
 *
 * The claim under test is the one the dock exists to make: that the pole figure
 * beside the structure is a figure of *that* structure as it is *now*. It is
 * checked the only way it can be — by turning the crystal and watching the poles
 * move — because a screenshot of a correct-looking pole figure is exactly what a
 * dock that never updated would also produce.
 */

/** Marker positions in one of the dock's three figures, as a comparable string. */
async function figureMarkers(page, index) {
  return page.evaluate((which) => {
    const canvas = document.querySelectorAll('.orient__canvas')[which];
    return [...canvas.querySelectorAll('circle, polygon')]
      .map((node) => node.getAttribute('cx') ?? node.getAttribute('points'))
      .join('|');
  }, index);
}

/** Drag across the structure, in the small steps a hand would make. */
async function turnTheStructure(page) {
  const surface = page.locator('#stage .plot svg').first();
  const box = await surface.boundingBox();
  const startX = box.x + box.width / 2;
  const startY = box.y + box.height / 2;
  await page.mouse.move(startX, startY);
  await page.mouse.down();
  for (let step = 1; step <= 10; step += 1) {
    await page.mouse.move(startX + step * 9, startY + step * 4);
  }
  await page.mouse.up();
}

test('turns the pole figure with the crystal', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  const dock = page.locator('.orient');
  await expect(dock).toBeVisible();
  // Three figures: the pole figure, the inverse pole figure, and the Kikuchi
  // map, which arrives a moment later because it is computed rather than drawn.
  await expect(page.locator('.orient__canvas svg')).toHaveCount(3, { timeout: 30_000 });

  const before = await figureMarkers(page, 0);
  expect(before.length).toBeGreaterThan(0);
  await turnTheStructure(page);
  const after = await figureMarkers(page, 0);

  expect(after).not.toEqual(before);
  // The fly-by leaves a trail through the standard triangle, so the inverse
  // pole figure gains markers rather than merely moving one.
  const trail = await page.evaluate(
    () => document.querySelectorAll('.orient__canvas')[1].querySelectorAll('circle').length,
  );
  expect(trail).toBeGreaterThan(4);

  expect(browserErrors).toEqual([]);
});

test('sets the view from Euler angles and reads them back', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  const angles = page.locator('.orient__angles input');
  await expect(angles).toHaveCount(3);

  // Goss, {011}<100>: the sheet normal is a <011> and the rolling direction a
  // <100>. Entering its angles must put exactly that in the readout, which is
  // the whole round trip -- angles to camera in Python, camera to figures in the
  // browser, camera back to a named direction in Python.
  await angles.nth(0).fill('0');
  await angles.nth(1).fill('45');
  await angles.nth(2).fill('0');
  await page.getByRole('button', { name: 'Set view', exact: true }).click();

  const readout = page.locator('.orient__readout');
  await expect(readout).toContainText('ND ∥ [101]');
  await expect(readout).toContainText('RD ∥ [001]');

  await turnTheStructure(page);
  // The angles follow the drag once it settles, so they no longer read as Goss.
  await expect.poll(() => angles.nth(1).inputValue(), { timeout: 10_000 }).not.toEqual('45.00');

  expect(browserErrors).toEqual([]);
});

test('offers named ideal orientations and applies one', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await page.getByRole('button', { name: 'brass', exact: true }).click();

  const angles = page.locator('.orient__angles input');
  await expect.poll(() => angles.nth(1).inputValue(), { timeout: 10_000 }).toEqual('45.00');
  await expect(angles.nth(0)).toHaveValue('35.26');
  // Brass is {011}<211>.
  await expect(page.locator('.orient__readout')).toContainText('ND ∥ [101]');
  await expect(page.locator('.orient__readout')).toContainText('RD ∥ [112]');

  expect(browserErrors).toEqual([]);
});

test('keeps the structure and both orientation figures on screen together', async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  const browserErrors = await openWorkbench(page);

  // The dock is a column beside the structure at this width, and both figures
  // fit in it: the point of the layout is that nothing has to be scrolled to to
  // be watched.
  const geometry = await page.evaluate(() => {
    const dock = document.querySelector('.orient');
    const dockBox = dock.getBoundingClientRect();
    const cells = [...dock.querySelectorAll('.orient__cell')].map((cell) =>
      cell.getBoundingClientRect(),
    );
    return {
      plotWidth: document.querySelector('#stage .plot').getBoundingClientRect().width,
      dockLeft: dockBox.left,
      dockRight: dockBox.right,
      dockBottom: dockBox.bottom,
      cellBottoms: cells.map((cell) => cell.bottom),
      overflows: dock.scrollHeight > dock.clientHeight + 1,
      bodyScrollWidth: document.body.scrollWidth,
      innerWidth: window.innerWidth,
    };
  });

  expect(geometry.dockLeft).toBeGreaterThan(geometry.plotWidth);
  expect(geometry.overflows).toBe(false);
  for (const bottom of geometry.cellBottoms) expect(bottom).toBeLessThanOrEqual(geometry.dockBottom);
  expect(geometry.bodyScrollWidth).toBeLessThanOrEqual(geometry.innerWidth);

  expect(browserErrors).toEqual([]);
});

test('folds the dock away on a phone without overflowing it', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  const browserErrors = await openWorkbench(page);

  const dock = page.locator('.orient');
  await expect(dock).toHaveJSProperty('open', false);
  await dock.locator('summary').click();
  await expect(dock).toHaveJSProperty('open', true);
  await expect(page.locator('.orient__canvas svg')).toHaveCount(3, { timeout: 30_000 });

  const overflow = await page.evaluate(() => document.body.scrollWidth - window.innerWidth);
  expect(overflow).toBeLessThanOrEqual(1);

  expect(browserErrors).toEqual([]);
});
