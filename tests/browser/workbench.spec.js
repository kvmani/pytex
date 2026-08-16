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
