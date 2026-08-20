import { readFileSync } from 'node:fs';

import { expect, test } from '@playwright/test';

const WORKSPACES = [
  'Crystal Viewer',
  'TEM Solver',
  'CBED',
  'Diffraction',
  'XRD',
  'EBSD',
  'Variants',
  'Texture',
  'Calculator',
];

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

  await page.getByRole('tab', { name: 'Calculator', exact: true }).click();
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

  await page.getByRole('tab', { name: 'Calculator', exact: true }).click();
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
  const overlay = patternControl(page, '.plot__overlay');
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

  /*
   * The picture is not blocked by the card that reports on it.
   *
   * The readout is pinned to the top-left corner of the drawing, which on a
   * diffraction pattern is a corner of the data. It is painted *under* the
   * figure, so the opaque part of the pattern masks it; bringing the pointer
   * onto its rectangle raises it in full, and leaving restores the clear view.
   */
  const stacking = await page.evaluate(() => {
    const card = document.querySelector('.plot__overlay');
    const canvas = document.querySelector('.plot__canvas');
    return {
      card: getComputedStyle(card).zIndex,
      canvas: getComputedStyle(canvas).zIndex,
      raised: card.classList.contains('plot__overlay--raised'),
    };
  });
  expect(Number(stacking.card)).toBeLessThan(Number(stacking.canvas));
  expect(stacking.raised).toBe(false);

  const hovered = await page.evaluate(() => {
    const card = document.querySelector('.plot__overlay');
    const svg = document.querySelector('svg[aria-label="Diffraction pattern"]');
    const box = card.getBoundingClientRect();
    const move = (clientX, clientY) =>
      svg.dispatchEvent(new PointerEvent('pointermove', { bubbles: true, clientX, clientY }));
    move(box.left + box.width / 2, box.top + box.height / 2);
    const onCard = card.classList.contains('plot__overlay--raised');
    move(box.right + 200, box.bottom + 200);
    return { onCard, afterLeaving: card.classList.contains('plot__overlay--raised') };
  });
  expect(hovered.onCard).toBe(true);
  expect(hovered.afterLeaving).toBe(false);

  // It reports the picks, so clearing them takes it away again.
  await page.getByRole('button', { name: 'Clear picks', exact: true }).click();
  await expect(overlay).toBeHidden();

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
  await page.getByRole('tab', { name: 'TEM Solver', exact: true }).click();
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
 * Measured pole figures, in tabs, on one scale.
 *
 * The reason to measure more than one pole figure is to compare them, and two
 * figures on two scales cannot be compared. So the assertions are that both
 * files arrive, that each gets a tab, that the tab switches the drawing, and
 * that the contour levels typed into the form are the ones the figure reports.
 */
test('XRDML pole figures open into tabs on one shared scale', async ({ page }) => {
  const browserErrors = await openWorkbench(page);
  await page.getByRole('tab', { name: 'Texture', exact: true }).click();

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

  // Two planes, one per file, then redraw.
  const poles = page.locator('[id^="ctl-poles-"]').first();
  await poles.fill('1 1 1\n2 0 0');
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
  await page.getByRole('tab', { name: 'EBSD', exact: true }).click();
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
  await page.getByRole('tab', { name: 'EBSD', exact: true }).click();
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
  await page.getByRole('tab', { name: 'TEM Solver', exact: true }).click();

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

/*
 * The crystal viewer's orientation dock.
 *
 * The claim under test is the one the dock exists to make: that the pole figure
 * beside the structure is a figure of *that* structure as it is *now*. It is
 * checked the only way it can be — by turning the crystal and watching the poles
 * move — because a screenshot of a correct-looking pole figure is exactly what a
 * dock that never updated would also produce.
 */

/** Marker positions in one of the dock's two figures, as a comparable string. */
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
  await expect(page.locator('.orient__canvas svg')).toHaveCount(2);

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
  await expect(page.locator('.orient__canvas svg')).toHaveCount(2);

  const overflow = await page.evaluate(() => document.body.scrollWidth - window.innerWidth);
  expect(overflow).toBeLessThanOrEqual(1);

  expect(browserErrors).toEqual([]);
});
