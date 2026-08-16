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

import {
  call,
  fetchManifest,
  fetchShell,
  ServiceCallError,
  setOperationTitles,
} from './core/api.js';
import { clear, el, markdown } from './core/dom.js';
import * as log from './core/logbook.js';
import { renderHelp, setExportFormats } from './core/result.js';
import { setPhaseCatalogue } from './core/phasecontrol.js';
import * as crystal from './panels/crystal.js';
import * as tem from './panels/tem.js';
import * as cbed from './panels/cbed.js';
import * as diffraction from './panels/diffraction.js';
import * as xrd from './panels/xrd.js';
import * as ebsd from './panels/ebsd.js';
import * as variants from './panels/variants.js';
import * as texture from './panels/texture.js';
import * as calculator from './panels/calculator.js';

// Order is the tab order. The viewer leads because it is the panel a newcomer
// understands without being told what it is for. Variants sits beside
// diffraction because the two answer the same question from opposite ends: the
// composite pattern is what variants look like on a plate, the pole figure is
// where they point — and texture follows variants because it is the same
// pole-figure reading applied to a whole polycrystal rather than one grain.
// CBED sits beside the TEM solver because it is the same specimen under a
// converged probe: the solver indexes a pattern of spots, and a CBED disc is
// what one of those spots becomes when the beam stops being parallel.
// EBSD follows XRD because it is the other way a texture is measured rather
// than modelled, and texture follows both because it is what those measurements
// are read as.
const PANELS = [crystal, tem, cbed, diffraction, xrd, ebsd, variants, texture, calculator];

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
  aboutDrawer: document.getElementById('about-drawer'),
  aboutBody: document.getElementById('about-body'),
  themeButton: document.getElementById('cycle-theme'),
  themeIcon: document.getElementById('theme-icon'),
  themeLabel: document.getElementById('theme-label'),
  console: document.getElementById('console'),
};

const app = {
  manifest: null,
  shell: null,
  active: null,
  mounted: null,
  index: [],
  logConsole: null,
};

app.logConsole = log.mountConsole(dom.console);
start();

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
  setOperationTitles(app.manifest.operations);
  log.notice(
    `PyTex ${app.manifest.version ?? ''} ready in the ${app.shell.shell} shell: ` +
      `${app.manifest.operations.length} operations across ${PANELS.length} workspaces.`,
    { source: 'app', detail: { shell: app.shell.shell } },
  );

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
  log.info(`Opened the ${module.panel.title} workspace.`, { source: panelId });
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

/* ------------------------------------------------------------------- about */

/**
 * Render the identity panel: what this program is, which build is running, who
 * wrote it, and the licence it is distributed under.
 *
 * Every fact comes from `manifest.about`, which Python builds in
 * `pytex/app/about.py` from the same version literal the package metadata reads.
 * Nothing here is written twice, so the version on screen is the version that
 * answered the request.
 */
function renderAbout(body, about) {
  if (!about) {
    body.append(...markdown('This build did not publish an About document.'));
    return;
  }
  const { author = {}, license = {}, links = [] } = about;
  body.append(
    el('h2', { text: `${about.name} ${about.version ?? ''}`.trim() }),
    about.tagline ? el('p.summary', { text: about.tagline }) : null,

    el('h3', { text: 'What it is' }),
    ...markdown(about.description ?? ''),
    about.audience ? el('p', { text: about.audience }) : null,

    el('h3', { text: 'Developed by' }),
    el('p.about__author', {}, [
      el('strong', { text: author.name ?? '' }),
      author.affiliation ? el('span', { text: author.affiliation }) : null,
      ...(author.emails ?? []).map((address) =>
        el('a', { href: `mailto:${address}`, text: address }),
      ),
    ]),

    el('h3', { text: 'Licence' }),
    el('p', {}, [
      el('strong', { text: license.name ?? '' }),
      license.spdx ? ' · ' : null,
      license.spdx ? el('code', { text: license.spdx }) : null,
    ]),
    license.notice ? el('p.about__notice', { text: license.notice }) : null,

    links.length ? el('h3', { text: 'Links' }) : null,
    links.length
      ? el(
          'ul.about__links',
          {},
          links.map((link) =>
            el('li', {}, [
              el('a', { href: link.url, target: '_blank', rel: 'noreferrer', text: link.label }),
            ]),
          ),
        )
      : null,

    el('h3', { text: 'This session' }),
    el('p.field__help', {
      text:
        `${app.manifest?.operations?.length ?? 0} operations across ${PANELS.length} ` +
        `workspaces, running in the ${app.shell?.shell ?? 'unknown'} shell.`,
    }),
  );
}

function openAbout() {
  clear(dom.aboutBody);
  renderAbout(dom.aboutBody, app.manifest?.about);
  dom.aboutDrawer.hidden = false;
  log.info('Opened the About panel.', {
    source: 'app',
    detail: { version: app.manifest?.about?.version ?? null },
  });
}

function closeAbout() {
  dom.aboutDrawer.hidden = true;
}

/* ------------------------------------------------------------------ errors */

/**
 * Surface a failure.
 *
 * `quiet` is for errors the panel has already placed beside the offending
 * control: repeating them as a toast would say the same thing twice, which
 * trains people to ignore both. It does *not* suppress the log entry — the
 * console is the record of the session, and an error the user saw and dismissed
 * is exactly the kind of thing they later need to find again.
 */
function showError(error, { quiet = false } = {}) {
  const isService = error instanceof ServiceCallError;
  // A ServiceCallError already produced a Python-side record on the way here;
  // logging it again would double every rejected input. Anything else is a
  // frontend fault that nothing else has reported.
  if (!isService) {
    log.error(String(error?.message ?? error), {
      source: app.active?.panel.id ?? 'app',
      detail: { kind: error?.name ?? 'Error' },
    });
  }
  if (quiet) return;
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
  log.critical(`PyTex could not start: ${error.message}`, { source: 'app' });
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
  document.getElementById('open-about').addEventListener('click', openAbout);
  document.getElementById('about-close').addEventListener('click', closeAbout);

  dom.palette.addEventListener('click', (event) => {
    if (event.target === dom.palette) closePalette();
  });
  dom.helpDrawer.addEventListener('click', (event) => {
    if (event.target === dom.helpDrawer) closeHelp();
  });
  dom.aboutDrawer.addEventListener('click', (event) => {
    if (event.target === dom.aboutDrawer) closeAbout();
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
    const message = event.detail?.message ?? 'Saved.';
    showNotice(message);
    // A toast lasts six seconds; where a file went is worth longer than that.
    log.notice(message, { source: app.active?.panel.id ?? 'app', detail: event.detail ?? {} });
  });

  document.addEventListener('keydown', (event) => {
    if ((event.ctrlKey || event.metaKey) && event.key.toLowerCase() === 'k') {
      event.preventDefault();
      openPalette();
    } else if (event.key === 'Escape') {
      closePalette();
      closeHelp();
      closeAbout();
    }
  });
}
