/**
 * The measured-texture panel: the inputs once, every reading of them in tabs.
 *
 * Quantitative texture analysis asks one set of inputs — the pole-figure files,
 * the plane of each, the crystal, the sample symmetry — for several answers:
 * the figures themselves, the orientation distribution, the fractions of named
 * components. A panel that rebuilt its form for each answer made the user state
 * the inputs again every time, and made it possible to read an ODF of one
 * symmetry beside figures of another. Here the rail holds the inputs, one
 * request answers everything, and the tabs only choose what to look at. The
 * service remembers the inversion, so choosing a different section returns
 * without solving again.
 *
 * Three rules from the model-texture panel carry over unchanged: every density
 * is in m.r.d. with the colour scale anchored at 1; the measured and
 * recalculated figures share one scale, because a comparison on two scales is
 * not one; and the difference is drawn on its own diverging scale centred on
 * zero, because "too high" and "too low" are the two things it is read for.
 */

import { el, formatNumber, svg } from '../core/dom.js';
import { explainer } from '../core/explainer.js';
import { buildForm } from '../core/controls.js';
import { plotFrame } from '../core/plotframe.js';
import { download, renderResult } from '../core/result.js';
import { call } from '../core/api.js';
import {
  closeFiles,
  describeOpenedFiles,
  onFilesChanged,
  openFiles,
  openedFiles,
} from '../core/texturefiles.js';
import {
  VIEW,
  contourLevels,
  contourPath,
  defaultContourStyle,
  discFrame,
  drawContourGrid,
  interpolatePoleFigure,
  paletteColor,
} from './texture.js';

export const panel = {
  id: 'texture_analysis',
  title: 'Measured texture',
  tagline:
    'Open pole figures and choose the symmetries once: measured, recalculated and difference '
    + 'figures, ODF sections and volume fractions follow.',
};

const OPERATION = 'texture.analysis';

const TABS = [
  { id: 'figures', label: 'Pole figures' },
  { id: 'odf', label: 'ODF sections' },
  { id: 'fractions', label: 'Volume fractions' },
];

/** What a pole-figure view can show. `plate` is every figure, every reading. */
const FIGURE_VIEWS = [
  { id: 'plate', label: 'All figures', title: 'Measured, recalculated and difference, for every figure' },
  { id: 'measured', label: 'Measured', title: 'As recorded, normalised to m.r.d.' },
  { id: 'symmetrized', label: 'Symmetrized', title: 'With the sample symmetry imposed: what was inverted' },
  { id: 'recalculated', label: 'Recalculated', title: 'Projected back from the ODF onto the measured directions' },
  { id: 'difference', label: 'Difference', title: 'Recalculated minus measured' },
];

const POINT_COLUMNS = [
  { key: 'measured', label: 'Measured', units: 'm.r.d.', numeric: true, digits: 3 },
  { key: 'symmetrized', label: 'Symmetrized', units: 'm.r.d.', numeric: true, digits: 3 },
  { key: 'recalculated', label: 'Recalculated', units: 'm.r.d.', numeric: true, digits: 3 },
  { key: 'difference', label: 'Difference', units: 'm.r.d.', numeric: true, digits: 3 },
  { key: 'polar_deg', label: 'From ND', units: '°', numeric: true, digits: 1 },
  { key: 'azimuth_deg', label: 'Azimuth', units: '°', numeric: true, digits: 1 },
];

const DIFFERENCE_NEGATIVE = [49, 84, 140];
const DIFFERENCE_NEUTRAL = [240, 240, 236];
const DIFFERENCE_POSITIVE = [154, 32, 32];

export function mount(context) {
  const operation = context.manifest.operations.find((entry) => entry.id === OPERATION);
  const examples = context.manifest.examples.filter((entry) => entry.panel === panel.id);
  const symbol = (name, fallback) => context.manifest.symbols?.[name]?.text ?? fallback;
  const state = {
    form: null,
    result: null,
    teaches: null,
    tab: 'figures',
    view: 'plate',
    figureIndex: 0,
    grids: new Map(),
    contour: defaultContourStyle(),
    plotNode: null,
    busy: false,
  };

  const frame = plotFrame({
    title: 'Texture analysis',
    toolbar: [
      el('button.button', {
        type: 'button',
        text: 'SVG',
        title: 'Save the figure on screen as SVG',
        onclick: () => {
          if (!state.plotNode) return;
          const markup = new XMLSerializer().serializeToString(state.plotNode);
          download(`pytex-texture-${state.tab}.svg`, markup, 'image/svg+xml');
        },
      }),
    ],
  });
  const mainTabs = el('div.figure-tabs.analysis-tabs', { role: 'tablist', 'aria-label': 'Reading' });
  const subBar = el('div.figure-tabs', { role: 'tablist', 'aria-label': 'View' });
  const legend = el('div.legend');
  const details = el('div');
  const formHost = el('div');
  const fileHost = el('div.group__body');
  const fileSummary = el('p.field__help');

  const runButton = el('button.button.button--primary.button--block', {
    type: 'button',
    text: 'Analyse texture',
    onclick: () => run(),
  });

  context.rail.append(
    el('details.group', { open: true, id: 'texture-analysis-files' }, [
      el('summary', { text: 'Open pole figures' }),
      fileHost,
    ]),
    formHost,
    runButton,
    el('details.group', { open: examples.length > 0 }, [
      el('summary', { text: 'Try an example' }),
      el('div.group__body', {}, [
        explainer(
          'Each example analyses demonstration figures measured from a known texture, so every '
            + 'number has an answer to be read against. Open your own files to analyse them '
            + 'instead; the examples always use the demonstration.',
          { label: 'What these examples show' },
        ),
        el(
          'div.examples',
          {},
          examples.map((example) =>
            el('button.example', { type: 'button', onclick: () => loadExample(example) }, [
              el('strong', { text: example.title }),
              el('span', { text: example.summary }),
            ]),
          ),
        ),
      ]),
    ]),
  );

  /* ----------------------------------------------------------------- rail */

  function renderFileControls() {
    const input = el('input', {
      type: 'file',
      accept: '.xrdml',
      multiple: true,
      'aria-label': 'Open XRDML pole-figure files',
      onchange: async (event) => {
        const chosen = [...(event.target.files ?? [])];
        if (!chosen.length) return;
        try {
          await openFiles(chosen);
          await run();
        } catch (error) {
          context.showError(error);
        }
      },
    });
    fileSummary.textContent = describeOpenedFiles(
      'No file open: the demonstration figures are analysed. Choose one or more .xrdml files.',
    );
    fileHost.replaceChildren(
      explainer(
        'Panalytical XRDML pole-figure files, one per reflection. Give the plane of each in '
          + '"Plane of each file", in the order they are listed here: the file records the '
          + 'diffraction angle, not the reflection. Files opened here are open in the Kearns '
          + 'parameter panel too.',
        { label: 'Which files, in which order' },
      ),
      input,
      fileSummary,
      openedFiles().length
        ? el('button.button', {
            type: 'button',
            text: 'Close them (use the demonstration)',
            onclick: () => {
              closeFiles();
              run();
            },
          })
        : null,
    );
  }

  function renderControls(initial = {}) {
    state.form = buildForm(operation, { initial });
    formHost.replaceChildren(state.form.element);
    // The opened files travel beside the form, from the control above.
    for (const field of state.form.element.querySelectorAll('.field')) {
      if (field.querySelector('[id^="ctl-files-"]')) field.hidden = true;
    }
  }

  function loadExample(example) {
    state.teaches = example.teaches;
    state.view = 'plate';
    state.figureIndex = 0;
    renderControls(example.request);
    run({ demonstration: true });
  }

  /* -------------------------------------------------------------- request */

  async function run({ demonstration = false } = {}) {
    if (state.busy) return;
    state.busy = true;
    runButton.disabled = true;
    runButton.textContent = 'Analysing…';
    frame.setStatus('Analysing: reading the figures, imposing the sample symmetry, inverting…');
    state.form.clearErrors();
    try {
      const request = { ...state.form.values() };
      if (!demonstration && openedFiles().length) request.files = { items: openedFiles() };
      else delete request.files;
      const result = await call(OPERATION, request);
      state.result = result;
      state.grids = new Map();
      state.figureIndex = Math.min(state.figureIndex, (result.data.figures?.length ?? 1) - 1);
      draw();
      renderResult(details, result, { teaches: state.teaches });
    } catch (error) {
      if (!state.form.showError(error)) context.showError(error);
      else context.showError(error, { quiet: true });
    } finally {
      state.busy = false;
      runButton.disabled = false;
      runButton.textContent = 'Analyse texture';
    }
  }

  /* ----------------------------------------------------------------- tabs */

  function renderMainTabs() {
    mainTabs.replaceChildren(
      ...TABS.map((tab) =>
        el('button.figure-tab', {
          type: 'button',
          role: 'tab',
          text: tab.label,
          'data-tab': tab.id,
          'aria-selected': String(tab.id === state.tab),
          onclick: () => {
            state.tab = tab.id;
            draw();
          },
        }),
      ),
    );
  }

  function tabButton(label, selected, onclick, title = '') {
    return el('button.figure-tab', {
      type: 'button',
      role: 'tab',
      text: label,
      title,
      'aria-selected': String(selected),
      onclick,
    });
  }

  /** The row under the main tabs: what exactly the chosen reading shows. */
  function renderSubBar() {
    const data = state.result?.data;
    if (!data) {
      subBar.replaceChildren();
      return;
    }
    if (state.tab === 'figures') {
      const symmetric = data.sample_symmetry !== 'triclinic';
      const views = FIGURE_VIEWS.filter((view) => symmetric || view.id !== 'symmetrized');
      const nodes = views.map((view) =>
        tabButton(view.label, view.id === state.view, () => {
          state.view = view.id;
          draw();
        }, view.title),
      );
      if (state.view !== 'plate') {
        nodes.push(el('span.analysis-tabs__gap', { 'aria-hidden': 'true', text: '·' }));
        data.figures.forEach((figure, index) => {
          nodes.push(
            tabButton(figure.label, index === state.figureIndex, () => {
              state.figureIndex = index;
              draw();
            }, figure.file),
          );
        });
      }
      subBar.replaceChildren(...nodes);
      return;
    }
    if (state.tab === 'odf') {
      const values = state.form.values();
      const kindSelect = el(
        'select',
        {
          'aria-label': 'Section kind',
          onchange: (event) => {
            state.form.setValues({ section_kind: event.currentTarget.value });
            run();
          },
        },
        choiceOptions('section_kind', values.section_kind),
      );
      const presetSelect = el(
        'select',
        {
          'aria-label': 'Which sections',
          onchange: (event) => {
            state.form.setValues({ section_preset: event.currentTarget.value });
            run();
          },
        },
        choiceOptions('section_preset', values.section_preset),
      );
      subBar.replaceChildren(
        el('label.analysis-tabs__control', {}, [el('span', { text: 'Sections' }), kindSelect]),
        el('label.analysis-tabs__control', {}, [el('span', { text: 'Which' }), presetSelect]),
      );
      return;
    }
    subBar.replaceChildren(
      el('span.field__help', {
        text: `Tolerance ${formatNumber(data.volume_fractions.tolerance_deg, 1)}° — change it `
          + 'under Volume fractions in the rail.',
      }),
    );
  }

  function choiceOptions(name, current) {
    const parameter = operation.parameters.find((entry) => entry.name === name);
    return (parameter?.options ?? []).map((option) => {
      const [value, label] = Array.isArray(option) ? option : [option.value, option.label];
      return el('option', { value, text: label, selected: value === current });
    });
  }

  /* -------------------------------------------------------------- drawing */

  function draw() {
    renderMainTabs();
    renderSubBar();
    const data = state.result?.data;
    if (!data) return;
    if (state.tab === 'figures') drawFigures(data);
    else if (state.tab === 'odf') drawSections(data);
    else drawFractions(data);
  }

  function grid(index, key) {
    const cacheKey = `${index}:${key}:${state.contour.gridSize}`;
    if (!state.grids.has(cacheKey)) {
      const points = state.result.data.figures[index].points.map((point) => ({
        x: point.x,
        y: point.y,
        mrd: point[key],
      }));
      state.grids.set(cacheKey, interpolatePoleFigure(points, state.contour.gridSize));
    }
    return state.grids.get(cacheKey);
  }

  function intensityStyle(data) {
    return {
      ...state.contour,
      customLevels: data.levels.join(' '),
      scaleMax: data.intensity_scale.maximum,
    };
  }

  /** One disc — an intensity or a difference — into `parent`, centred at the origin. */
  function drawDisc(parent, index, key, data) {
    const clipId = `ta-clip-${index}-${key}-${Math.random().toString(36).slice(2, 7)}`;
    const defs = svg('defs');
    const clip = svg('clipPath', { id: clipId });
    clip.append(svg('circle', { cx: 0, cy: 0, r: VIEW }));
    defs.append(clip);
    parent.append(defs);
    const body = svg('g', { 'clip-path': `url(#${clipId})` });
    const mapPoint = (x, y) => ({ x: x * VIEW, y: -y * VIEW });
    if (key === 'difference') {
      drawDifferenceGrid(body, grid(index, key), data.difference_scale.maximum, state.contour, mapPoint);
    } else {
      drawContourGrid(body, grid(index, key), intensityStyle(data), data.intensity_scale.maximum, mapPoint);
    }
    parent.append(body);
    discFrame(parent, data.specimen_axes ?? ['RD', 'TD']);
    return body;
  }

  function drawFigures(data) {
    const figures = data.figures;
    if (state.view === 'plate') {
      frame.configure({ toData: null, formatCursor: null });
      state.plotNode = renderPlate(data);
      frame.setContent(state.plotNode);
      renderLegends(data, true);
      frame.setStatus(
        `${figures.length} figure(s) · ${data.sample_symmetry_label} sample symmetry · RP `
          + `${formatNumber(data.mean_fit_rp_percent, 1)}% to the inverted figures, `
          + `${formatNumber(data.mean_rp_percent, 1)}% to the measured · one intensity scale, `
          + `difference on its own ±${formatNumber(data.difference_scale.maximum, 2)} m.r.d. scale`,
      );
      return;
    }
    const index = Math.min(state.figureIndex, figures.length - 1);
    const figure = figures[index];
    const key = state.view;
    frame.configure({
      toData: (x, y) => ({ x: x / VIEW, y: -y / VIEW }),
      formatCursor: (point) => {
        const radius = Math.hypot(point.x, point.y);
        if (radius > 1.0001) return 'outside the projection';
        const polar = data.projection === 'equal_area'
          ? 2 * Math.asin(Math.min(radius / Math.SQRT2, 1))
          : 2 * Math.atan(radius);
        const azimuth = ((Math.atan2(point.y, point.x) * 180) / Math.PI + 360) % 360;
        return `${formatNumber((polar * 180) / Math.PI, 1)}° from ND · ${formatNumber(azimuth, 1)}° azimuth`;
      },
    });
    const root = svg('svg', {
      viewBox: `${-VIEW * 1.16} ${-VIEW * 1.16} ${2.32 * VIEW} ${2.32 * VIEW}`,
      preserveAspectRatio: 'xMidYMid meet',
      'aria-label': `${figure.label} ${key} pole figure`,
      'data-figure-view': key,
    });
    const body = drawDisc(root, index, key, data);
    const hitRadius = (2 * VIEW) / Math.sqrt(Math.max(figure.points.length, 1));
    for (const point of figure.points) {
      const node = svg('circle', {
        cx: point.x * VIEW,
        cy: -point.y * VIEW,
        r: hitRadius,
        fill: 'transparent',
        'pointer-events': 'all',
      });
      body.append(node);
      frame.hoverable(node, point, POINT_COLUMNS);
    }
    state.plotNode = root;
    frame.setContent(root);
    renderLegends(data, key === 'difference' ? 'difference' : 'intensity');
    const range = figure.ranges[key];
    const residual = figure.residual;
    frame.setStatus(
      `${figure.label} from ${figure.file} · ${FIGURE_VIEWS.find((view) => view.id === key).title} · `
        + `${formatNumber(range.minimum, 2)} to ${formatNumber(range.maximum, 2)} m.r.d. · RP `
        + `${formatNumber(figure.fit_residual.rp_percent, 1)}% to inverted, `
        + `${formatNumber(residual.rp_percent, 1)}% to measured · hover for every reading at a point`,
    );
  }

  /** Every figure as a row: measured, (symmetrized), recalculated, difference. */
  function renderPlate(data) {
    const columns = [
      ['measured', 'Measured'],
      ...(data.sample_symmetry !== 'triclinic' ? [['symmetrized', 'Symmetrized']] : []),
      ['recalculated', 'Recalculated'],
      ['difference', 'Difference'],
    ];
    const cellWidth = 2.5 * VIEW;
    const cellHeight = 2.55 * VIEW;
    const top = 0.35 * VIEW;
    const left = 0.9 * VIEW;
    const width = left + columns.length * cellWidth;
    const height = top + data.figures.length * cellHeight;
    const root = svg('svg', {
      viewBox: `0 0 ${width} ${height}`,
      preserveAspectRatio: 'xMidYMid meet',
      'aria-label': `${data.figures.length} pole figures: measured, recalculated and difference`,
      'data-plate-rows': String(data.figures.length),
      'data-plate-columns': String(columns.length),
    });
    columns.forEach(([, label], column) => {
      root.append(svg('text', {
        x: left + (column + 0.5) * cellWidth,
        y: top * 0.7,
        'text-anchor': 'middle',
        'font-size': 14,
        'font-weight': 600,
        fill: 'currentColor',
        text: label,
      }));
    });
    data.figures.forEach((figure, row) => {
      root.append(svg('text', {
        x: left * 0.45,
        y: top + (row + 0.5) * cellHeight,
        'text-anchor': 'middle',
        'font-size': 14,
        'font-weight': 600,
        fill: 'currentColor',
        text: figure.label,
      }));
      columns.forEach(([key], column) => {
        const cell = svg('g', {
          transform: `translate(${left + (column + 0.5) * cellWidth} ${top + (row + 0.5) * cellHeight})`,
          'data-plate-cell': `${figure.label}:${key}`,
        });
        drawDisc(cell, row, key, data);
        root.append(cell);
      });
    });
    return root;
  }

  function drawSections(data) {
    const odf = data.odf;
    frame.configure({ toData: null, formatCursor: null });
    state.plotNode = renderSectionPlate(odf);
    frame.setContent(state.plotNode);
    legend.replaceChildren(renderRamp(odf.max_mrd, { ...state.contour, customLevels: '', scaleMax: 0 }));
    const across = symbol(odf.horizontal_coordinate, odf.horizontal_coordinate === 'phi2' ? 'φ₂' : 'φ₁');
    frame.setStatus(
      `${odf.sections.length} ${odf.section_kind_label.toLowerCase()} section(s) · ${across} `
        + `across, Φ down · peak ${formatNumber(odf.max_mrd, 2)} m.r.d. · ${odf.method_label} · `
        + `range φ₁ 0–${odf.ranges.phi1_max_deg}°, Φ 0–${odf.ranges.big_phi_max_deg}°, `
        + `φ₂ 0–${odf.ranges.phi2_max_deg}° for this crystal and sample symmetry`,
    );
  }

  function renderSectionPlate(odf) {
    const sections = odf.sections;
    const count = sections.length;
    const columns = count <= 4 ? count : Math.ceil(Math.sqrt(count * 1.4));
    const rows = Math.ceil(count / columns);
    const across = sections[0].across_deg;
    const bigPhi = sections[0].big_phi_deg;
    const aspect = Math.min(Math.max(across.at(-1) / Math.max(bigPhi.at(-1), 1), 0.5), 4);
    const cellHeight = 100;
    const cellWidth = cellHeight * aspect;
    const gapX = 16;
    const gapY = 26;
    const style = { ...state.contour, customLevels: '', scaleMax: 0 };
    const width = columns * cellWidth + (columns - 1) * gapX;
    const height = rows * cellHeight + rows * gapY;
    const root = svg('svg', {
      viewBox: `-8 -4 ${width + 16} ${height + 8}`,
      preserveAspectRatio: 'xMidYMid meet',
      'aria-label': `${count} ODF sections`,
      'data-odf-sections': String(count),
      'data-section-kind': odf.section_kind,
    });
    const constant = {
      phi2: symbol('phi2', 'φ₂'),
      phi1: symbol('phi1', 'φ₁'),
      sigma: 'σ',
    }[odf.section_kind];
    const readHover = count <= 6;
    const hoverColumns = [
      { key: 'section', label: `${constant}`, units: '°', numeric: true, digits: 1 },
      { key: 'across', label: odf.horizontal_coordinate === 'phi2' ? 'φ₂' : 'φ₁', units: '°', numeric: true, digits: 1 },
      { key: 'big_phi', label: 'Φ', units: '°', numeric: true, digits: 1 },
      { key: 'mrd', label: 'Density', units: 'm.r.d.', numeric: true, digits: 3 },
    ];
    sections.forEach((section, index) => {
      const column = index % columns;
      const row = Math.floor(index / columns);
      const originX = column * (cellWidth + gapX);
      const originY = row * (cellHeight + gapY) + gapY - 6;
      const cell = svg('g', { 'data-section-value': String(section.value_deg) });
      const values = section.densities;
      drawContourGrid(
        cell,
        {
          xValues: across.map((_, i) => i / Math.max(across.length - 1, 1)),
          yValues: bigPhi.map((_, i) => i / Math.max(bigPhi.length - 1, 1)),
          values,
        },
        style,
        odf.max_mrd,
        (x, y) => ({ x: originX + x * cellWidth, y: originY + y * cellHeight }),
      );
      cell.append(
        svg('rect', {
          x: originX, y: originY, width: cellWidth, height: cellHeight,
          fill: 'none', stroke: 'currentColor', 'stroke-opacity': 0.5, 'stroke-width': 0.6,
        }),
        svg('text', {
          x: originX + cellWidth / 2,
          y: originY - 5,
          'text-anchor': 'middle',
          'font-size': 9,
          fill: 'currentColor',
          text: `${constant} = ${formatNumber(section.value_deg, 1)}° · max ${formatNumber(section.max_mrd, 1)}`,
        }),
      );
      if (readHover) {
        const stepX = cellWidth / Math.max(across.length - 1, 1);
        const stepY = cellHeight / Math.max(bigPhi.length - 1, 1);
        bigPhi.forEach((phiValue, i) => {
          across.forEach((acrossValue, j) => {
            const node = svg('rect', {
              x: originX + j * stepX - stepX / 2,
              y: originY + i * stepY - stepY / 2,
              width: stepX,
              height: stepY,
              fill: 'transparent',
              'pointer-events': 'all',
            });
            cell.append(node);
            frame.hoverable(node, {
              section: section.value_deg, across: acrossValue, big_phi: phiValue, mrd: values[i][j],
            }, hoverColumns);
          });
        });
      }
      root.append(cell);
    });
    return root;
  }

  function drawFractions(data) {
    const rows = data.volume_fractions.rows;
    frame.configure({ toData: null, formatCursor: null });
    legend.replaceChildren(
      el('span.legend__item', {}, [
        el('span.legend__swatch', { style: 'background:var(--accent)' }),
        el('span', { text: 'volume within the tolerance' }),
      ]),
      el('span.legend__item', {}, [
        el('span', { text: '┆ what a random texture holds in the same ball' }),
      ]),
    );
    if (!rows.length) {
      state.plotNode = null;
      frame.setContent(el('div.stage__placeholder', {
        text: 'No ideal-orientation catalogue is defined for this crystal system.',
      }));
      frame.setStatus('');
      return;
    }
    state.plotNode = renderFractionChart(rows, data.volume_fractions.tolerance_deg);
    frame.setContent(state.plotNode);
    const top = [...rows].sort((a, b) => b.fraction - a.fraction)[0];
    frame.setStatus(
      `${rows.length} ideal orientations within ${formatNumber(data.volume_fractions.tolerance_deg, 1)}° · `
        + `strongest ${top.component} at ${formatNumber(top.percent, 1)}%, `
        + `${formatNumber(top.times_random, 2)}× random · components are independent, so the `
        + 'fractions need not sum to one',
    );
  }

  function renderFractionChart(rows, tolerance) {
    const width = 520;
    const rowHeight = 30;
    const left = 190;
    const right = 150;
    const top = 26;
    const height = top + rows.length * rowHeight + 10;
    const maximum = Math.max(...rows.map((row) => Math.max(row.percent, row.random_percent)), 1);
    const x = (value) => left + (value / maximum) * (width - left - right);
    const root = svg('svg', {
      viewBox: `0 0 ${width} ${height}`,
      preserveAspectRatio: 'xMidYMid meet',
      'aria-label': `Volume fractions within ${tolerance} degrees`,
      'data-fraction-rows': String(rows.length),
    });
    root.append(svg('text', {
      x: left, y: 14, 'font-size': 10, fill: 'currentColor', 'fill-opacity': 0.7,
      text: `Volume within ${formatNumber(tolerance, 1)}° (%)`,
    }));
    rows.forEach((row, index) => {
      const y = top + index * rowHeight + rowHeight / 2;
      root.append(
        svg('text', {
          x: left - 8, y: y - 2, 'text-anchor': 'end', 'font-size': 11, 'font-weight': 600,
          fill: 'currentColor', text: row.component,
        }),
        svg('text', {
          x: left - 8, y: y + 10, 'text-anchor': 'end', 'font-size': 8, fill: 'currentColor',
          'fill-opacity': 0.65, text: row.miller,
        }),
        svg('rect', {
          x: left, y: y - 7, width: Math.max(x(row.percent) - left, 0.5), height: 14,
          fill: 'var(--accent)', 'fill-opacity': 0.8,
        }, [svg('title', { text: `${row.component}: ${formatNumber(row.percent, 2)}%` })]),
        svg('line', {
          x1: x(row.random_percent), x2: x(row.random_percent), y1: y - 11, y2: y + 11,
          stroke: 'currentColor', 'stroke-width': 1.2, 'stroke-dasharray': '2 2',
        }),
        svg('text', {
          x: x(row.percent) + 6, y: y + 4, 'font-size': 10, fill: 'currentColor',
          text: `${formatNumber(row.percent, 1)}% · ${formatNumber(row.times_random, 2)}× random`,
        }),
      );
    });
    return root;
  }

  /* -------------------------------------------------------------- legends */

  function renderRamp(maximum, style) {
    const stops = contourLevels(maximum, style);
    return el('span.legend__group', {}, [
      el('span.legend__item', {}, [el('span', { text: 'm.r.d.' })]),
      ...stops.map((value) =>
        el('span.legend__item', {}, [
          el('span.legend__swatch', { style: `background:${paletteColor(value, style, maximum)}` }),
          el('span', { text: value === 1 ? '1 (random)' : formatNumber(value, 2) }),
        ]),
      ),
    ]);
  }

  function renderDifferenceKey(maximum) {
    const levels = differenceLevels(maximum);
    return el('span.legend__group', {}, [
      el('span.legend__item', {}, [el('span', { text: 'recalculated − measured' })]),
      ...[-maximum, ...levels, maximum].map((value) =>
        el('span.legend__item', {}, [
          el('span.legend__swatch', { style: `background:${divergingColor(value, maximum)}` }),
          el('span', { text: formatNumber(value, 2) }),
        ]),
      ),
    ]);
  }

  function renderLegends(data, which) {
    const nodes = [];
    if (which === true || which === 'intensity') {
      nodes.push(renderRamp(data.intensity_scale.maximum, intensityStyle(data)));
    }
    if (which === true || which === 'difference') {
      nodes.push(renderDifferenceKey(data.difference_scale.maximum));
    }
    legend.replaceChildren(...nodes);
  }

  renderFileControls();
  renderControls();
  renderMainTabs();
  frame.setControls(mainTabs, subBar, legend);
  frame.setContent(el('div.stage__placeholder', { text: 'Press Analyse texture, or try an example.' }));
  context.stage.append(frame.element, details);
  const unsubscribe = onFilesChanged(() => renderFileControls());
  if (examples.length) loadExample(examples[0]);

  return {
    help: () => operation,
    unmount: () => unsubscribe(),
  };
}

/* ------------------------------------------------------------------ helpers */

function divergingColor(value, maximum) {
  const t = Math.max(-1, Math.min(1, value / (maximum || 1)));
  const target = t < 0 ? DIFFERENCE_NEGATIVE : DIFFERENCE_POSITIVE;
  const amount = Math.abs(t);
  const mix = DIFFERENCE_NEUTRAL.map((channel, index) =>
    Math.round(channel + (target[index] - channel) * amount));
  return `rgb(${mix[0]} ${mix[1]} ${mix[2]})`;
}

/** Symmetric, round contour levels for a difference figure, excluding zero. */
function differenceLevels(maximum) {
  if (!(maximum > 0)) return [];
  const raw = maximum / 4;
  const magnitude = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((factor) => factor * magnitude).find((value) => value >= raw);
  const levels = [];
  for (let value = step; value < maximum - 1e-12; value += step) levels.push(-value, value);
  return levels.sort((a, b) => a - b);
}

/**
 * A difference figure: diverging fill centred on zero, dashed lines below it.
 *
 * Its own scale on purpose. Drawn on the m.r.d. ramp a difference would read as
 * a texture — and a residual of +0.3 m.r.d. is not "a weak pole", it is the ODF
 * over-predicting by 0.3 at that direction.
 */
function drawDifferenceGrid(parent, grid, maximum, style, mapPoint) {
  const fill = svg('g', { 'fill-opacity': style.fillOpacity });
  for (let row = 0; row < grid.yValues.length - 1; row += 1) {
    for (let column = 0; column < grid.xValues.length - 1; column += 1) {
      const samples = [
        grid.values[row][column], grid.values[row][column + 1],
        grid.values[row + 1][column], grid.values[row + 1][column + 1],
      ].filter(Number.isFinite);
      if (!samples.length) continue;
      const value = samples.reduce((sum, sample) => sum + sample, 0) / samples.length;
      const a = mapPoint(grid.xValues[column], grid.yValues[row]);
      const b = mapPoint(grid.xValues[column + 1], grid.yValues[row + 1]);
      fill.append(svg('rect', {
        x: Math.min(a.x, b.x),
        y: Math.min(a.y, b.y),
        width: Math.abs(b.x - a.x) + 0.35,
        height: Math.abs(b.y - a.y) + 0.35,
        fill: divergingColor(value, maximum),
      }));
    }
  }
  parent.append(fill);
  const lines = svg('g', { 'aria-label': 'Difference contour lines' });
  for (const level of differenceLevels(maximum)) {
    const d = contourPath(grid, level, mapPoint);
    if (!d) continue;
    lines.append(svg('path', {
      d,
      fill: 'none',
      stroke: style.lineColor,
      'stroke-width': style.lineWidth,
      'stroke-dasharray': level < 0 ? '3 2' : null,
      'vector-effect': 'non-scaling-stroke',
      'data-contour-level': level,
    }));
  }
  parent.append(lines);
}
