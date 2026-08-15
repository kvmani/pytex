/**
 * Shared controls for presentation-only plot choices.
 *
 * Scientific parameters belong in the operation manifest and go to Python.
 * Marker shape, colour and visual scale do not change a diffraction result, so
 * they stay in the shared frontend and redraw the rows Python already made.
 */

import { el } from './dom.js';

export const DEFAULT_MARKER_STYLE = Object.freeze({
  shape: 'circle',
  fill: 'filled',
  variantEncoding: 'color_shape_size',
  scale: 1,
  sizeMode: 'perceptual',
  parentColor: '#64748b',
  palette: 'distinct',
  childColor: '#2563eb',
});

const SHAPES = [
  ['circle', 'Circle'],
  ['square', 'Square'],
  ['triangle', 'Triangle'],
  ['diamond', 'Diamond'],
  ['star', 'Star'],
  ['cross', 'Cross'],
];

const FILLS = [
  ['filled', 'Filled'],
  ['outline', 'Outline (unfilled)'],
];

const VARIANT_ENCODINGS = [
  ['color', 'Colour'],
  ['shape', 'Shape'],
  ['size', 'Size'],
  ['color_shape', 'Shape + colour'],
  ['color_size', 'Size + colour'],
  ['shape_size', 'Shape + size'],
  ['color_shape_size', 'Shape + size + colour (24 distinct)'],
];

const SIZE_MODES = [
  ['perceptual', 'Perceptual intensity'],
  ['area', 'Area follows intensity'],
  ['constant', 'Constant size'],
];

const PALETTES = [
  ['distinct', 'Distinct variants'],
  ['accessible', 'Colourblind-safe cycle'],
  ['single', 'Single product colour'],
];

function selectField(label, help, options, value, oninput) {
  const select = el('select', { oninput });
  select.append(
    ...options.map(([optionValue, optionLabel]) =>
      el('option', { value: optionValue, text: optionLabel, selected: optionValue === value }),
    ),
  );
  return el('label.field', {}, [
    el('span.field__label', { text: label }),
    select,
    el('span.field__hint', { text: help }),
  ]);
}

function colorField(label, help, value, oninput) {
  const input = el('input', { type: 'color', value, oninput });
  return {
    input,
    element: el('label.field', {}, [
      el('span.field__label', { text: label }),
      el('span.color-control', {}, [input, el('output', { text: value.toUpperCase() })]),
      el('span.field__hint', { text: help }),
    ]),
  };
}

/**
 * Build marker controls and return their live style object.
 *
 * @param {object} options
 * @param {(style: object) => void} options.onChange redraw callback
 */
export function markerStyleControl({ onChange }) {
  const style = { ...DEFAULT_MARKER_STYLE };
  const scaleOutput = el('output', { text: '1.00×', for: 'marker-scale' });
  const scaleInput = el('input', {
    id: 'marker-scale',
    type: 'range',
    min: 0.5,
    max: 2.5,
    step: 0.05,
    value: style.scale,
    oninput: (event) => {
      style.scale = Number(event.currentTarget.value);
      scaleOutput.textContent = `${style.scale.toFixed(2)}×`;
      onChange(style);
    },
  });

  const parent = colorField(
    'Parent colour',
    'The achromatic reference lattice. Choose a colour with enough contrast for the current theme.',
    style.parentColor,
    (event) => {
      style.parentColor = event.currentTarget.value;
      event.currentTarget.nextElementSibling.textContent = style.parentColor.toUpperCase();
      onChange(style);
    },
  );
  const child = colorField(
    'Product colour',
    'Used when colour is not a variant channel or the palette uses one product colour.',
    style.childColor,
    (event) => {
      style.childColor = event.currentTarget.value;
      event.currentTarget.nextElementSibling.textContent = style.childColor.toUpperCase();
      onChange(style);
    },
  );
  const usesPaletteColor = () =>
    style.variantEncoding.includes('color') && style.palette !== 'single';
  child.input.disabled = usesPaletteColor();

  const root = el('details.group.appearance', {}, [
    el('summary', { text: 'Appearance' }),
    el('div.group__body', {}, [
      el('p.field__help', {
        text: 'These controls change only the drawing. Spot coordinates, indices, intensities and exports stay unchanged.',
      }),
      selectField(
        'Spot shape',
        'Use shape as a second visual channel or to match a publication convention.',
        SHAPES,
        style.shape,
        (event) => {
          style.shape = event.currentTarget.value;
          onChange(style);
        },
      ),
      selectField(
        'Spot fill',
        'Outline leaves spots unfilled. Double-diffraction spots retain a dashed edge in either mode.',
        FILLS,
        style.fill,
        (event) => {
          style.fill = event.currentTarget.value;
          onChange(style);
        },
      ),
      selectField(
        'Variant encoding',
        'Choose shape, size, colour, or a combination. Shape + size + colour guarantees 24 distinct symbols even if a palette repeats.',
        VARIANT_ENCODINGS,
        style.variantEncoding,
        (event) => {
          style.variantEncoding = event.currentTarget.value;
          child.input.disabled = usesPaletteColor();
          onChange(style);
        },
      ),
      el('label.field', {}, [
        el('span.field__label', { text: 'Spot-size scale' }),
        el('span.range-control', {}, [scaleInput, scaleOutput]),
        el('span.field__hint', { text: 'Scales every marker without changing detector coordinates.' }),
      ]),
      selectField(
        'Intensity sizing',
        'Perceptual keeps weak reflections visible; area is quantitatively literal; constant compares geometry alone.',
        SIZE_MODES,
        style.sizeMode,
        (event) => {
          style.sizeMode = event.currentTarget.value;
          onChange(style);
        },
      ),
      parent.element,
      selectField(
        'Variant palette',
        'Distinct maximises separation; the safe cycle avoids red-green ambiguity; single colour emphasises phase over variant.',
        PALETTES,
        style.palette,
        (event) => {
          style.palette = event.currentTarget.value;
          child.input.disabled = usesPaletteColor();
          onChange(style);
        },
      ),
      child.element,
      el('button.button', {
        type: 'button',
        text: 'Reset appearance',
        onclick: () => {
          Object.assign(style, DEFAULT_MARKER_STYLE);
          root.replaceWith(markerStyleControl({ onChange }).element);
          onChange(style);
        },
      }),
    ]),
  ]);

  return { element: root, style };
}

const ACCESSIBLE_COLORS = [
  '#0072b2', '#e69f00', '#009e73', '#cc79a7',
  '#d55e00', '#56b4e9', '#f0e442', '#7c3aed',
];

/** Return the colour for one product source under a marker style. */
export function productColor(index, style) {
  if (!style.variantEncoding.includes('color') || style.palette === 'single') {
    return style.childColor;
  }
  if (style.palette === 'accessible') return ACCESSIBLE_COLORS[index % ACCESSIBLE_COLORS.length];
  const hue = (index * 137.508 + 212) % 360;
  return `hsl(${hue.toFixed(1)} 72% 54%)`;
}

/** Shape and size channels for one product variant. */
export function variantMarkerStyle(index, style) {
  const useShape = style.variantEncoding.includes('shape');
  const useSize = style.variantEncoding.includes('size');
  const shape = useShape ? SHAPES[index % SHAPES.length][0] : style.shape;
  const sizeCycle = [0.82, 0.94, 1.06, 1.18];
  const sizeIndex = Math.floor(index / SHAPES.length) % sizeCycle.length;
  return { shape, scale: useSize ? sizeCycle[sizeIndex] : 1 };
}

/** Convert relative intensity to marker radius in plot-view units. */
export function markerRadius(intensity, style) {
  const safe = Math.max(Number(intensity) || 0, 0);
  let radius = 2.2;
  if (style.sizeMode === 'perceptual') radius = 0.7 + 2.6 * Math.pow(safe, 0.25);
  else if (style.sizeMode === 'area') radius = 0.7 + 2.6 * Math.sqrt(safe);
  return radius * style.scale;
}

/** Build an SVG marker centred on x/y. */
export function markerNode(svg, {
  x, y, radius, shape, color, hollow = false, dashed = false,
}) {
  const paint = {
    fill: hollow || shape === 'cross' ? 'none' : color,
    stroke: hollow || shape === 'cross' ? color : 'none',
    'stroke-width': shape === 'cross' ? Math.max(radius * 0.38, 0.5) : 0.4,
    'stroke-dasharray': dashed && shape !== 'cross' ? '1 1' : null,
  };
  if (shape === 'square') {
    return svg('rect', { x: x - radius, y: y - radius, width: radius * 2, height: radius * 2, ...paint });
  }
  if (shape === 'diamond') {
    return svg('polygon', {
      points: `${x},${y - radius} ${x + radius},${y} ${x},${y + radius} ${x - radius},${y}`,
      ...paint,
    });
  }
  if (shape === 'triangle') {
    const halfWidth = radius * 0.92;
    return svg('polygon', {
      points: `${x},${y - radius} ${x + halfWidth},${y + radius * 0.72} ` +
        `${x - halfWidth},${y + radius * 0.72}`,
      ...paint,
    });
  }
  if (shape === 'star') {
    const points = Array.from({ length: 10 }, (_, index) => {
      const angle = -Math.PI / 2 + index * Math.PI / 5;
      const pointRadius = index % 2 === 0 ? radius : radius * 0.42;
      return `${x + Math.cos(angle) * pointRadius},${y + Math.sin(angle) * pointRadius}`;
    }).join(' ');
    return svg('polygon', { points, ...paint });
  }
  if (shape === 'cross') {
    return svg('path', {
      d: `M ${x - radius} ${y} L ${x + radius} ${y} M ${x} ${y - radius} L ${x} ${y + radius}`,
      ...paint,
    });
  }
  return svg('circle', { cx: x, cy: y, r: radius, ...paint });
}
