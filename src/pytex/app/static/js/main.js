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

import { call, fetchManifest, ServiceCallError } from './core/api.js';
import { clear, el, markdown } from './core/dom.js';
import { renderHelp } from './core/result.js';
import { setPhaseCatalogue } from './core/phasecontrol.js';
import * as calculator from './panels/calculator.js';

const PANELS = [calculator];

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
};

const app = {
  manifest: null,
  active: null,
  mounted: null,
  index: [],
};

start();

async function start() {
  try {
    app.manifest = await fetchManifest();
  } catch (error) {
    showFatal(error);
    return;
  }

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
