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
  fetchExperience,
  fetchManifest,
  fetchShell,
  ServiceCallError,
  setOperationTitles,
} from './core/api.js';
import { clear, el, markdown } from './core/dom.js';
import { mountFeedback } from './core/feedback.js';
import * as log from './core/logbook.js';
import { createTour } from './core/tour.js';
import { renderHelp, setExportFormats } from './core/result.js';
import { setPhaseCatalogue } from './core/phasecontrol.js';
import * as crystal from './panels/crystal.js';
import * as saedsim from './panels/saedsim.js';
import * as tem from './panels/tem.js';
import * as cbed from './panels/cbed.js';
import * as diffraction from './panels/diffraction.js';
import * as xrd from './panels/xrd.js';
import * as ebsd from './panels/ebsd.js';
import * as ebsdDistribution from './panels/ebsddistribution.js';
import * as ebsdFigures from './panels/ebsdfigures.js';
import * as ebsdSummary from './panels/ebsdsummary.js';
import * as ebsdKikuchi from './panels/ebsdkikuchi.js';
import * as ecciWorkflow from './panels/ecciWorkflow.js';
import * as variants from './panels/variants.js';
import * as texture from './panels/texture.js';
import * as calculator from './panels/calculator.js';

// Order is the tab order. The viewer leads because it is the panel a newcomer
// understands without being told what it is for. Variants sits beside
// diffraction because the two answer the same question from opposite ends: the
// composite pattern is what variants look like on a plate, the pole figure is
// where they point — and texture follows variants because it is the same
// pole-figure reading applied to a whole polycrystal rather than one grain.
// EBSD follows XRD because it is the other way a texture is measured rather
// than modelled, and texture follows both because it is what those measurements
// are read as.
//
// A *workspace* is a top-level tab. Most are one panel, and read exactly as
// they did when the tab bar was a flat list of panels. One is a group: every
// transmission-electron surface is one technique seen four ways — a pattern
// simulated, a pattern indexed, a disc under a converged probe, and the
// composite a set of variants makes — and four top-level tabs said they were
// four subjects. Grouped, the tab bar names the subject and the sub-tab bar
// names the view, which is what they are.
const TEM_ANALYSIS = {
  id: 'tem-analysis',
  title: 'TEM Analysis',
  tagline: 'Simulate, index, and interpret transmission-electron diffraction.',
  panels: [saedsim, tem, cbed, diffraction],
};

// EBSD is grouped for the same reason: one scan, seen six ways — followed by
// two forward/experiment-planning tools. The Kikuchi simulator and ECCI
// workflow take no scan; both start from a standalone phase and orientation,
// and belong here because their geometry is the EBSD geometry and nowhere
// else's. The scan itself is shared between the first six (see
// `core/ebsdscan.js`), because a file opened in the summary is open for the map
// — anything else would silently analyse the practice dataset next to a user's
// own data.
const EBSD_ANALYSIS = {
  id: 'ebsd',
  title: 'EBSD',
  tagline: 'One orientation scan: its map, its statistics, and where it points.',
  panels: [
    ebsd,
    ebsd.grodPanel,
    ebsd.kamPanel,
    ebsdSummary,
    ebsdDistribution,
    ebsdFigures,
    ebsdKikuchi,
    ecciWorkflow,
  ],
};

const WORKSPACES = [
  solo(crystal),
  TEM_ANALYSIS,
  solo(xrd),
  EBSD_ANALYSIS,
  solo(variants),
  solo(texture),
  solo(calculator),
];

/** A workspace that is one panel: the tab is the panel, as it always was. */
function solo(module) {
  return {
    id: module.panel.id,
    title: module.panel.title,
    tagline: module.panel.tagline,
    panels: [module],
  };
}

/** Every panel, in tab order, for the counts the log and the About panel state. */
const PANELS = WORKSPACES.flatMap((workspace) => workspace.panels);

/** The workspace a panel id belongs to, and the panel itself. */
function locate(panelId) {
  for (const workspace of WORKSPACES) {
    const module = workspace.panels.find((entry) => entry.panel.id === panelId);
    if (module) return { workspace, module };
  }
  return null;
}

const THEMES = {
  auto: { label: 'Auto', icon: '◐', description: 'follow the system' },
  light: { label: 'Light', icon: '☀', description: 'light' },
  dark: { label: 'Dark', icon: '☾', description: 'dark' },
};

const dom = {
  tabs: document.getElementById('tabs'),
  subtabs: document.getElementById('subtabs'),
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
  feedbackDrawer: document.getElementById('feedback-drawer'),
  feedbackBody: document.getElementById('feedback-body'),
  feedbackButton: document.getElementById('open-feedback'),
  themeButton: document.getElementById('cycle-theme'),
  themeIcon: document.getElementById('theme-icon'),
  themeLabel: document.getElementById('theme-label'),
  console: document.getElementById('console'),
};

const app = {
  manifest: null,
  shell: null,
  experience: null, // what this deployment decided about greeting and feedback
  active: null, // the open workspace
  activePanel: null, // the open panel inside it
  mounted: null,
  index: [],
  logConsole: null,
  feedback: null,
  tour: null,
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
  // Fetched before anything is drawn so that the Feedback button is either
  // there from the start or never appears, rather than arriving late.
  app.experience = await fetchExperience();
  setExportFormats(app.manifest.export_formats);
  setOperationTitles(app.manifest.operations);
  log.notice(
    `PyTex ${app.manifest.version ?? ''} ready in the ${app.shell.shell} shell: ` +
      `${app.manifest.operations.length} operations across ${PANELS.length} panels ` +
      `in ${WORKSPACES.length} workspaces.`,
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
  activate(WORKSPACES[0].id);

  // After the first panel is mounted, not before: the tour points at the tab
  // bar, the control rail and the console, and a welcome that arrives while
  // the page is still empty is talking about furniture that is not there yet.
  wireExperience();
}

/**
 * Turn on the two things the deployment decides: the feedback form and the
 * welcome tour.
 *
 * Both are optional, and a server that does not publish `/api/experience` at
 * all — an older build, or a proxy that drops the route — leaves both off.
 * That is the safe failure: an unexplained modal on startup is worse than no
 * tour, and a feedback button that posts into the void is worse than none.
 */
function wireExperience() {
  const feedback = app.experience?.feedback ?? { enabled: false };
  if (feedback.enabled) {
    app.feedback = mountFeedback({
      drawer: dom.feedbackDrawer,
      body: dom.feedbackBody,
      config: feedback,
      // Gathered at submit time rather than at open time: people open the form,
      // go back to look at the thing they are describing, and come back.
      context: () => ({
        workspace: app.active?.title ?? null,
        panel: app.activePanel?.panel.id ?? null,
        shell: app.shell?.shell ?? null,
        viewport: `${window.innerWidth}x${window.innerHeight}`,
        user_agent: navigator.userAgent,
      }),
    });
    dom.feedbackButton.addEventListener('click', () => app.feedback.open());
    document.getElementById('feedback-close').addEventListener('click', () => app.feedback.close());
    dom.feedbackDrawer.addEventListener('click', (event) => {
      if (event.target === dom.feedbackDrawer) app.feedback.close();
    });
  } else {
    dom.feedbackButton.hidden = true;
  }

  app.tour = createTour({
    config: app.experience?.tour ?? { enabled: false },
    steps: tourSteps(),
    version: app.manifest?.version ?? 'unknown',
  });
  if (app.tour.shouldGreet()) app.tour.start();
}

/**
 * What the tour says, in order.
 *
 * Every step points at something on the page, and a step whose target is not
 * present is dropped rather than shown pointing at nothing — which is what
 * makes it safe to name the feedback button here even in a deployment that
 * turned the form off.
 */
function tourSteps() {
  const names = WORKSPACES.map((workspace) => workspace.title);
  return [
    {
      title: 'Welcome to the PyTex Workbench',
      body: [
        'This is the PyTex library with a face on it — the same code a script would call, ' +
          'so anything you work out here you can reproduce in a few lines of Python.',
        'Two minutes now will save you hunting later. You can leave at any point.',
      ],
      note: 'Every result states the convention it was computed under. Nothing here is a black box.',
    },
    {
      title: 'Seven workspaces, one data model',
      target: '#tabs',
      body: [
        `Along the top: ${names.join(', ')}. Crystal structures, electron and X-ray ` +
          'diffraction, EBSD maps, variant and orientation-relationship analysis, and texture.',
        'They share one description of frames, symmetry and orientations, so a phase you set up ' +
          'in one is the same phase in all of them.',
      ],
    },
    {
      title: 'The controls come from the code',
      target: '#rail',
      body: [
        'Every control on the right is generated from the same declaration Python validates ' +
          'against, down to its help text — press the ? beside any field.',
        'Miller indices get one box per index, so 1 1 0 can never be read as something else.',
      ],
    },
    {
      title: 'Results you can take away',
      target: '#stage',
      body: [
        'Figures are drawn here, and every result can be exported as CSV, XLSX or JSON — the ' +
          'numbers you are looking at, not a second calculation of them.',
      ],
    },
    {
      title: 'Find anything by name',
      target: '#open-palette',
      body: [
        'Press Ctrl+K, or use this button, to search every operation and every worked example ' +
          'the build carries. It is generated from the manifest, so it is never out of date.',
      ],
      note: 'Start with an example if you are not sure what a panel does — they all run.',
    },
    {
      title: 'The session keeps a record',
      target: '#console',
      body: [
        'This strip is the message log: what ran, what it was given, and what it warned about. ' +
          'It is the place to look when a number is not what you expected.',
      ],
    },
    {
      title: 'Tell us what is missing',
      target: '#open-feedback',
      body: [
        'PyTex grows from what its users say they need. This button opens a short form — a ' +
          'feature you want, something that reads wrongly, or a paper you wish were implemented.',
        'Nothing is too small, and you need not leave your name.',
      ],
      note: 'That is the tour. Choose a workspace and have a look around.',
    },
  ];
}

function buildTabs() {
  clear(dom.tabs);
  for (const workspace of WORKSPACES) {
    dom.tabs.append(
      el('button.tab', {
        type: 'button',
        role: 'tab',
        text: workspace.title,
        'aria-selected': 'false',
        dataset: { workspace: workspace.id },
        onclick: () => activate(workspace.id),
      }),
    );
  }
}

/**
 * Draw the sub-tab strip for the open workspace.
 *
 * Hidden entirely when the workspace is one panel, rather than shown with a
 * single tab in it: a navigation control offering one destination is furniture,
 * and it would cost every ungrouped workspace a strip of screen to say nothing.
 */
function buildSubtabs(workspace) {
  clear(dom.subtabs);
  if (workspace.panels.length < 2) {
    dom.subtabs.hidden = true;
    return;
  }
  dom.subtabs.hidden = false;
  for (const module of workspace.panels) {
    const { id, title, tagline } = module.panel;
    dom.subtabs.append(
      el('button.subtab', {
        type: 'button',
        role: 'tab',
        text: title,
        title: tagline,
        'aria-selected': String(app.activePanel?.panel.id === id),
        dataset: { panel: id },
        onclick: () => activate(workspace.id, id),
      }),
    );
  }
}

/**
 * Open a workspace, and one panel inside it.
 *
 * `target` names the workspace or, for a grouped one, the panel to land on;
 * omitting it keeps the panel already open in that workspace if there is one,
 * and otherwise takes the first. The palette calls this with an operation's
 * panel id, which may be a sub-panel — so both are accepted, and the workspace
 * is looked up from whichever was given.
 */
function activate(target, panelId = null) {
  const workspace =
    WORKSPACES.find((entry) => entry.id === target) ?? locate(target)?.workspace ?? null;
  if (!workspace) return;
  const wanted = panelId ?? (workspace.id === target ? null : target);
  const module =
    workspace.panels.find((entry) => entry.panel.id === wanted) ??
    (workspace.panels.includes(app.activePanel) ? app.activePanel : null) ??
    workspace.panels[0];

  app.active = workspace;
  app.activePanel = module;
  for (const tab of dom.tabs.children) {
    tab.setAttribute('aria-selected', String(tab.dataset.workspace === workspace.id));
  }
  buildSubtabs(workspace);
  dom.tagline.textContent = module.panel.tagline;
  log.info(`Opened the ${module.panel.title} workspace.`, { source: module.panel.id });
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
  // The tour is remembered as seen, which would make it unrepeatable; offering
  // it from the help panel is what makes remembering safe rather than final.
  if (app.tour) {
    dom.helpBody.append(
      el('p.field__help', {}, [
        el('button.button', {
          type: 'button',
          text: 'Show the welcome tour again',
          onclick: () => {
            closeHelp();
            app.tour.start();
          },
        }),
      ]),
    );
  }
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
        `panels in ${WORKSPACES.length} workspaces, running in the ` +
        `${app.shell?.shell ?? 'unknown'} shell.`,
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
      source: app.activePanel?.panel.id ?? 'app',
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
    log.notice(message, { source: app.activePanel?.panel.id ?? 'app', detail: event.detail ?? {} });
  });

  document.addEventListener('keydown', (event) => {
    if ((event.ctrlKey || event.metaKey) && event.key.toLowerCase() === 'k') {
      event.preventDefault();
      openPalette();
    } else if (event.key === 'Escape') {
      closePalette();
      closeHelp();
      closeAbout();
      app.feedback?.close();
    }
  });
}
