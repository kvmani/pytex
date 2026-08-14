/**
 * The shell: tabs, the command palette, the help drawer, and error surfacing.
 *
 * It knows nothing about crystallography. It fetches the manifest, builds the
 * tab bar from the panels the manifest declares, mounts the matching panel
 * module, and provides the three services every panel needs — report an error,
 * show help, open the palette.
 *
 * Feature discovery lives here and is generated: the palette indexes every
 * operation, every example and every tag in the manifest, so a capability
 * registered in Python is findable the moment it exists, without anyone
 * remembering to add it to a menu.
 */

import { call, fetchManifest, ServiceCallError, fetchShell } from './core/api.js';
import { clear, el, markdown } from './core/dom.js';
import { renderHelp, setExportFormats } from './core/result.js';
import { setPhaseCatalogue } from './core/phasecontrol.js';
import * as crystal from './panels/crystal.js';
import * as tem from './panels/tem.js';
import * as diffraction from './panels/diffraction.js';
import * as xrd from './panels/xrd.js';
import * as variants from './panels/variants.js';
import * as texture from './panels/texture.js';
import * as calculator from './panels/calculator.js';

// Order is the tab order. The viewer leads because it is the panel a newcomer
// understands without being told what it is for. Variants sits beside
// diffraction because the two answer the same question from opposite ends: the
// composite pattern is what variants look like on a plate, the pole figure is
// where they point — and texture follows variants because it is the same
// pole-figure reading applied to a whole polycrystal rather than one grain.
const PANELS = [crystal, tem, diffraction, xrd, variants, texture, calculator];

const THEMES = {
  auto: { label: 'Auto', icon: '◐', description: 'follow the system' },
  light: { label: 'Light', icon: '☀', description: 'light' },
  dark: { label: 'Dark', icon: '☾', description: 'dark' },
};

const dom = {
  tabs: document.getElementById('tabs'),
  stage: document.getElementById('stage'),
  rail: document.getElementById('rail-body'),
  tagline: document.getElementById('masthead-tagline'),
  toasts: document.getElementById('toasts'),
  palette: document.getElementById('palette'),
  paletteInput: document.getElementById('palette-input'),
  paletteResults: document.getElementById('palette-results'),
  helpDrawer: document.getElementById('help-drawer'),
  helpBody: document.getElementById('help-body'),
  themeButton: document.getElementById('cycle-theme'),
  themeIcon: document.getElementById('theme-icon'),
  themeLabel: document.getElementById('theme-label'),
  activity: document.getElementById('activity'),
  activityToggle: document.getElementById('activity-toggle'),
  activityPanel: document.getElementById('activity-panel'),
  activityIndicator: document.getElementById('activity-indicator'),
  activitySummary: document.getElementById('activity-summary'),
  activityCount: document.getElementById('activity-count'),
  activityLog: document.getElementById('activity-log'),
  activityEmpty: document.getElementById('activity-empty'),
  activityClear: document.getElementById('activity-clear'),
};

const app = {
  manifest: null,
  shell: null,
  active: null,
  mounted: null,
  index: [],
  activity: { active: new Map(), history: [] },
};

wireActivity();
start();

function operationTitle(operation) {
  return app.manifest?.operations.find((entry) => entry.id === operation)?.title ?? operation;
}

function durationLabel(durationMs) {
  return durationMs < 1000 ? `${Math.round(durationMs)} ms` : `${(durationMs / 1000).toFixed(1)} s`;
}

function renderActivity() {
  const running = app.activity.active.size;
  dom.activity.classList.toggle('activity--busy', running > 0);
  dom.activityIndicator.classList.toggle('activity__indicator--busy', running > 0);
  dom.activitySummary.textContent = running
    ? `Running ${[...app.activity.active.values()].at(-1)}…`
    : app.activity.history[0]
      ? `${app.activity.history[0].title} ${app.activity.history[0].outcome}`
      : 'Ready';
  const completed = app.activity.history.length;
  dom.activityCount.textContent = running
    ? `${running} active · ${completed} recent`
    : completed
      ? `${completed} recent calculation${completed === 1 ? '' : 's'}`
      : 'No calculations yet';
  dom.activityEmpty.hidden = completed > 0;
  clear(dom.activityLog);
  for (const entry of app.activity.history) {
    dom.activityLog.append(
      el(`li.activity__entry.activity__entry--${entry.outcome}`, {}, [
        el('span.activity__entry-mark', { text: entry.outcome === 'completed' ? '✓' : '!' }),
        el('span', {}, [
          el('strong', { text: entry.title }),
          entry.message ? el('small', { text: entry.message }) : null,
        ]),
        el('time', { text: durationLabel(entry.durationMs) }),
      ]),
    );
  }
}

function wireActivity() {
  document.addEventListener('pytex:operation-start', (event) => {
    app.activity.active.set(event.detail.id, operationTitle(event.detail.operation));
    renderActivity();
  });
  document.addEventListener('pytex:operation-finish', (event) => {
    app.activity.active.delete(event.detail.id);
    app.activity.history.unshift({
      title: operationTitle(event.detail.operation),
      outcome: event.detail.outcome,
      durationMs: event.detail.durationMs,
      message: event.detail.message,
    });
    app.activity.history = app.activity.history.slice(0, 40);
    renderActivity();
  });
  dom.activityToggle.addEventListener('click', () => {
    const open = dom.activityPanel.hidden;
    dom.activityPanel.hidden = !open;
    dom.activityToggle.setAttribute('aria-expanded', String(open));
    dom.activityToggle.setAttribute(
      'aria-label',
      `${open ? 'Close' : 'Open'} calculation activity`,
    );
  });
  dom.activityClear.addEventListener('click', () => {
    app.activity.history = [];
    renderActivity();
  });
  renderActivity();
}

function savedTheme() {
  try {
    const value = localStorage.getItem('pytex-theme') ?? 'auto';
    return Object.hasOwn(THEMES, value) ? value : 'auto';
  } catch {
    return 'auto';
  }
}

function applyTheme(theme) {
  const chosen = THEMES[theme] ?? THEMES.auto;
  document.documentElement.dataset.theme = theme in THEMES ? theme : 'auto';
  dom.themeIcon.textContent = chosen.icon;
  dom.themeLabel.textContent = chosen.label;
  dom.themeButton.title = `Colour theme: ${chosen.description}`;
  dom.themeButton.setAttribute('aria-label', `Colour theme: ${chosen.description}`);
}

function cycleTheme() {
  const order = Object.keys(THEMES);
  const current = document.documentElement.dataset.theme || 'auto';
  const next = order[(order.indexOf(current) + 1) % order.length];
  try {
    localStorage.setItem('pytex-theme', next);
  } catch {
    // A locked-down intranet webview may disable storage; the theme still
    // changes for the current session, which is the useful part.
  }
  applyTheme(next);
}

applyTheme(savedTheme());

async function start() {
  try {
    app.manifest = await fetchManifest();
  } catch (error) {
    showFatal(error);
    return;
  }

  // Which shell this is decides only how a file is saved, and the answer comes
  // from Python rather than from sniffing the window, so that a shell which
  // changes how it writes files says so in one place.
  app.shell = await fetchShell();
  setExportFormats(app.manifest.export_formats);

  // The phase picker needs the catalogue before any control renders, so it is
  // fetched once here rather than lazily by each panel that shows a phase.
  try {
    const catalogue = await call('calc.catalog');
    setPhaseCatalogue({
      phases: catalogue.data.phases,
      pointGroups: catalogue.data.point_groups,
    });
  } catch (error) {
    showError(error);
  }

  buildIndex();
  buildTabs();
  wireGlobals();
  activate(PANELS[0].panel.id);
}

function buildTabs() {
  clear(dom.tabs);
  for (const module of PANELS) {
    const { id, title } = module.panel;
    dom.tabs.append(
      el('button.tab', {
        type: 'button',
        role: 'tab',
        text: title,
        'aria-selected': 'false',
        dataset: { panel: id },
        onclick: () => activate(id),
      }),
    );
  }
}

function activate(panelId) {
  const module = PANELS.find((entry) => entry.panel.id === panelId);
  if (!module) return;
  app.active = module;
  for (const tab of dom.tabs.children) {
    tab.setAttribute('aria-selected', String(tab.dataset.panel === panelId));
  }
  dom.tagline.textContent = module.panel.tagline;
  clear(dom.stage);
  clear(dom.rail);
  app.mounted = module.mount({
    stage: dom.stage,
    rail: dom.rail,
    manifest: app.manifest,
    showError,
    openHelp,
  });
}

/* ------------------------------------------------------------------ search */

function buildIndex() {
  app.index = [
    ...app.manifest.operations.map((operation) => ({
      kind: 'operation',
      id: operation.id,
      title: operation.title,
      subtitle: operation.summary,
      panel: operation.panel,
      haystack: [operation.title, operation.summary, operation.id, ...(operation.tags ?? [])]
        .join(' ')
        .toLowerCase(),
      run: () => {
        activate(operation.panel);
        openHelp(operation);
      },
    })),
    ...app.manifest.examples.map((example) => ({
      kind: 'example',
      id: example.id,
      title: example.title,
      subtitle: example.summary,
      panel: example.panel,
      haystack: [example.title, example.summary, example.teaches].join(' ').toLowerCase(),
      run: () => activate(example.panel),
    })),
  ];
}

function search(query) {
  const needle = query.trim().toLowerCase();
  if (!needle) return app.index.slice(0, 12);
  const terms = needle.split(/\s+/);
  return app.index
    .filter((entry) => terms.every((term) => entry.haystack.includes(term)))
    .slice(0, 20);
}

function openPalette() {
  dom.palette.hidden = false;
  dom.paletteInput.value = '';
  renderPalette(search(''));
  dom.paletteInput.focus();
}

function closePalette() {
  dom.palette.hidden = true;
}

function renderPalette(entries) {
  clear(dom.paletteResults);
  entries.forEach((entry, index) => {
    dom.paletteResults.append(
      el(
        'li',
        {
          role: 'option',
          'aria-selected': String(index === 0),
          dataset: { index },
          onclick: () => {
            closePalette();
            entry.run();
          },
        },
        [
          el('span.palette__kind', { text: entry.kind }),
          el('strong', { text: entry.title }),
          el('span', { text: entry.subtitle }),
        ],
      ),
    );
  });
  if (!entries.length) {
    dom.paletteResults.append(el('li', { text: 'Nothing matches that.' }));
  }
}

function moveSelection(delta) {
  const items = [...dom.paletteResults.children];
  if (!items.length) return;
  const current = items.findIndex((item) => item.getAttribute('aria-selected') === 'true');
  const next = Math.min(Math.max(current + delta, 0), items.length - 1);
  items.forEach((item, index) => item.setAttribute('aria-selected', String(index === next)));
  items[next].scrollIntoView({ block: 'nearest' });
}

function runSelection() {
  const chosen = [...dom.paletteResults.children].find(
    (item) => item.getAttribute('aria-selected') === 'true',
  );
  if (!chosen) return;
  closePalette();
  chosen.click();
}

/* -------------------------------------------------------------------- help */

function openHelp(operation) {
  const target = operation ?? app.mounted?.help?.();
  clear(dom.helpBody);
  if (target) renderHelp(dom.helpBody, target);
  else dom.helpBody.append(...markdown('Choose a calculation to see what it does.'));
  dom.helpDrawer.hidden = false;
}

function closeHelp() {
  dom.helpDrawer.hidden = true;
}

/* ------------------------------------------------------------------ errors */

/**
 * Surface a failure.
 *
 * `quiet` is for errors the panel has already placed beside the offending
 * control: repeating them as a toast would say the same thing twice, which
 * trains people to ignore both.
 */
function showError(error, { quiet = false } = {}) {
  if (quiet) return;
  const isService = error instanceof ServiceCallError;
  const toast = el('div.toast', {}, [
    el('strong', { text: isService ? error.message : 'Something went wrong.' }),
    isService && error.hint ? el('span', { text: error.hint }) : null,
    !isService ? el('span', { text: String(error?.message ?? error) }) : null,
  ]);
  dom.toasts.append(toast);
  setTimeout(() => toast.remove(), 9000);
  toast.addEventListener('click', () => toast.remove());
}

function showNotice(text) {
  const toast = el('div.toast.toast--ok', {}, [el('strong', { text })]);
  dom.toasts.append(toast);
  setTimeout(() => toast.remove(), 6000);
  toast.addEventListener('click', () => toast.remove());
}

function showFatal(error) {
  clear(dom.stage);
  dom.stage.append(
    el('section.card', {}, [
      el('div.card__body', {}, [
        el('h2.card__title', { text: 'PyTex could not start' }),
        el('p.summary', { text: error.message }),
        error.hint ? el('p.field__help', { text: error.hint }) : null,
      ]),
    ]),
  );
}

/* ----------------------------------------------------------------- globals */

function wireGlobals() {
  dom.themeButton.addEventListener('click', cycleTheme);
  document.getElementById('open-palette').addEventListener('click', openPalette);
  document.getElementById('open-help').addEventListener('click', () => openHelp());
  document.getElementById('help-close').addEventListener('click', closeHelp);

  dom.palette.addEventListener('click', (event) => {
    if (event.target === dom.palette) closePalette();
  });
  dom.helpDrawer.addEventListener('click', (event) => {
    if (event.target === dom.helpDrawer) closeHelp();
  });

  dom.paletteInput.addEventListener('input', () => renderPalette(search(dom.paletteInput.value)));
  dom.paletteInput.addEventListener('keydown', (event) => {
    if (event.key === 'ArrowDown') {
      event.preventDefault();
      moveSelection(1);
    } else if (event.key === 'ArrowUp') {
      event.preventDefault();
      moveSelection(-1);
    } else if (event.key === 'Enter') {
      event.preventDefault();
      runSelection();
    }
  });

  // Exports say where the file went. The desktop shell writes through a native
  // dialog and the web shell hands the browser a download, and neither is
  // visible from inside the page, so the one thing that must never happen is
  // pressing an export button and being told nothing at all.
  document.addEventListener('pytex:saved', (event) => {
    showNotice(event.detail?.message ?? 'Saved.');
  });

  document.addEventListener('keydown', (event) => {
    if ((event.ctrlKey || event.metaKey) && event.key.toLowerCase() === 'k') {
      event.preventDefault();
      openPalette();
    } else if (event.key === 'Escape') {
      closePalette();
      closeHelp();
    }
  });
}
