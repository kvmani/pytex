/**
 * The phase picker, which is the one control every panel shares.
 *
 * A phase is chosen far more often than anything else in this application, so
 * it gets a purpose-built control rather than a JSON box: a catalogue dropdown,
 * a CIF-file choice, and six cell parameters plus a point group for everything
 * else. Picking a catalogue entry fills the parameter fields in, so "start from
 * zirconium and stretch c" is two clicks and one edit rather than a retyped
 * structure. A CIF remains file-backed until Python normalizes it into the same
 * canonical phase contract; no crystallography is reimplemented in JavaScript.
 *
 * The control emits exactly what `PhaseSpec.from_json` accepts: `{builtin: id}`
 * when nothing was edited, `{cif: {name, text}}` for a file, and a full
 * description once anything was. That asymmetry is deliberate — an unedited
 * catalogue choice should stay a *reference*, so a later correction to a
 * published lattice parameter reaches saved analyses.
 */

import { el } from './dom.js';
import { explainer } from './explainer.js';

const CATALOGUE = { phases: [], pointGroups: [] };

/** Load the catalogue once, at start-up, so the control renders synchronously. */
export function setPhaseCatalogue({ phases, pointGroups }) {
  CATALOGUE.phases = phases ?? [];
  CATALOGUE.pointGroups = pointGroups ?? [];
}

export function phaseControl(parameter, value, onChange, id) {
  const state = { spec: normalise(value), edited: false, cif: cifOf(value) };

  const cifStatus = el('p.field__help', {
    text: state.cif ? `${state.cif.name} loaded; it replaces the catalogue phase.` : '',
    hidden: !state.cif,
  });
  const cifInput = el('input', {
    type: 'file',
    accept: '.cif,chemical/x-cif',
    'aria-label': 'Load a CIF crystal structure',
    onchange: (event) => loadCif(event.target.files?.[0]),
  });
  const removeCif = el('button.button', {
    type: 'button',
    text: 'Remove CIF',
    hidden: !state.cif,
    onclick: () => {
      state.cif = null;
      cifInput.value = '';
      cifStatus.hidden = true;
      removeCif.hidden = true;
      catalogue.disabled = false;
      details.hidden = false;
      writeFields(state.spec);
      onChange();
    },
  });

  const catalogue = el(
    'select',
    {
      id,
      onchange: () => {
        const chosen = CATALOGUE.phases.find((entry) => entry.id === catalogue.value);
        if (!chosen) return;
        clearCif();
        state.spec = { ...chosen };
        state.edited = false;
        writeFields(state.spec);
        onChange();
      },
    },
    [
      el('option', { value: '', text: 'Custom phase…' }),
      ...CATALOGUE.phases.map((entry) =>
        el('option', {
          value: entry.id,
          text: `${entry.name} · ${entry.point_group}`,
          title: entry.source ?? '',
        }),
      ),
    ],
  );

  const cellInputs = {};
  const cellRow = (names, step) =>
    el(
      'div.cell-row',
      {},
      names.map((name) => {
        const input = el('input', {
          type: 'number',
          step,
          'aria-label': name,
          placeholder: name,
          oninput: () => {
            state.edited = true;
            onChange();
          },
        });
        cellInputs[name] = input;
        return el('label.cell-cell', {}, [el('span', { text: name }), input]);
      }),
    );

  const pointGroup = el(
    'select',
    {
      'aria-label': 'Point group',
      onchange: () => {
        state.edited = true;
        // The space group came from the catalogue entry this edit started from,
        // and a space group belongs to exactly one crystal system. Left in
        // place across a symmetry change it decides which reflections are
        // absent for a phase it does not describe -- a cubic F number quietly
        // deleting most of a tetragonal pattern. Clearing it computes without
        // systematic absences, which is the honest default for a phase whose
        // space group is no longer known.
        if (crystalSystemOf(pointGroup.value) !== crystalSystemOf(state.spec?.point_group)) {
          spaceGroup.value = '';
          spaceGroupNumber.value = '';
        }
        onChange();
      },
    },
    CATALOGUE.pointGroups.map((group) =>
      el('option', {
        value: group.symbol,
        text: `${group.symbol} · ${group.crystal_system}`,
      }),
    ),
  );

  const spaceGroup = el('input', {
    type: 'text',
    placeholder: 'Fm-3m',
    'aria-label': 'Space-group symbol',
    oninput: () => {
      state.edited = true;
      onChange();
    },
  });
  const spaceGroupNumber = el('input', {
    type: 'number',
    min: 1,
    max: 230,
    placeholder: '225',
    'aria-label': 'Space-group number',
    oninput: () => {
      state.edited = true;
      onChange();
    },
  });

  const details = el('details.group', {}, [
    el('summary', { text: 'Cell parameters and symmetry' }),
    el('div.group__body', {}, [
      cellRow(['a', 'b', 'c'], 'any'),
      cellRow(['alpha', 'beta', 'gamma'], 'any'),
      el('label.field__label', { text: 'Point group' }),
      pointGroup,
      el('label.field__label', { text: 'Space group (sets systematic absences)' }),
      el('div.cell-row', {}, [spaceGroup, spaceGroupNumber]),
    ]),
  ]);

  const cifGroup = el('div.phase-cif', {}, [
    // A label, then the detail behind it. The sentence this replaces named all
    // four things a CIF brings and ran to three lines in a rail this narrow,
    // every time, for a control most sessions never touch.
    el('span.field__label', { text: 'Or load a .cif file' }),
    explainer(
      'Its lattice, symmetry, space group and atomic sites replace the catalogue ' +
        'phase entirely. The cell and symmetry controls above are disabled while one is ' +
        'loaded, because they would then describe something the calculation is not using.',
      { label: 'What a CIF replaces' },
    ),
    cifInput,
    cifStatus,
    removeCif,
  ]);
  const element = el('div.phase-control', {}, [catalogue, cifGroup, details]);

  async function loadCif(file) {
    if (!file) return;
    state.cif = { name: file.name, text: '' };
    cifStatus.hidden = false;
    cifStatus.textContent = `Reading ${file.name}…`;
    removeCif.hidden = true;
    try {
      state.cif.text = await file.text();
      catalogue.value = '';
      catalogue.disabled = true;
      details.hidden = true;
      cifStatus.textContent =
        `${file.name} loaded. Run the analysis to validate and normalize its crystal structure.`;
      removeCif.hidden = false;
      onChange();
    } catch {
      state.cif = null;
      cifInput.value = '';
      cifStatus.textContent = `${file.name} could not be read in the browser.`;
      removeCif.hidden = true;
    }
  }

  function clearCif() {
    state.cif = null;
    cifInput.value = '';
    cifStatus.hidden = true;
    removeCif.hidden = true;
    catalogue.disabled = false;
    details.hidden = false;
  }

  function writeFields(spec) {
    if (!spec) return;
    for (const name of ['a', 'b', 'c', 'alpha', 'beta', 'gamma']) {
      if (cellInputs[name]) cellInputs[name].value = spec[name] ?? '';
    }
    if (spec.point_group) pointGroup.value = spec.point_group;
    spaceGroup.value = spec.space_group_symbol ?? '';
    spaceGroupNumber.value = spec.space_group_number ?? '';
    catalogue.value = spec.id ?? '';
  }

  writeFields(state.spec);
  if (state.cif) {
    catalogue.value = '';
    catalogue.disabled = true;
    details.hidden = true;
  }

  return {
    id,
    element,
    read() {
      if (state.cif) return { cif: state.cif };
      const chosen = catalogue.value;
      if (chosen && !state.edited) return { builtin: chosen };
      const spec = {
        // An edited catalogue phase is not that phase. Titles, summaries and
        // the provenance sheet of every export take this name, so a 4 x 4 x 6
        // tetragonal cell must not go out labelled "Nickel (fcc)" because the
        // edit happened to start there.
        name: state.spec?.name ? `${state.spec.name} (edited)` : 'Custom phase',
        a: numberOf(cellInputs.a),
        b: numberOf(cellInputs.b),
        c: numberOf(cellInputs.c),
        alpha: numberOf(cellInputs.alpha),
        beta: numberOf(cellInputs.beta),
        gamma: numberOf(cellInputs.gamma),
        point_group: pointGroup.value,
      };
      if (spaceGroup.value.trim()) {
        spec.space_group_symbol = spaceGroup.value.trim();
        spec.space_group_number = numberOf(spaceGroupNumber);
      }
      // Carry the atomic basis of the catalogue entry the edit started from:
      // changing a lattice parameter should not silently discard the structure
      // and turn every intensity into an unavailable result.
      if (chosen && state.spec?.sites) spec.sites = state.spec.sites;
      return spec;
    },
    write(next) {
      state.cif = cifOf(next);
      state.spec = normalise(next);
      state.edited = Boolean(next && !next.builtin && !next.id && !next.cif);
      if (state.cif) {
        cifStatus.textContent = `${state.cif.name} loaded; it replaces the catalogue phase.`;
        cifStatus.hidden = false;
        removeCif.hidden = false;
        catalogue.disabled = true;
        details.hidden = true;
      } else {
        clearCif();
      }
      writeFields(state.spec);
    },
  };
}

function numberOf(input) {
  if (!input || input.value === '') return null;
  return Number(input.value);
}

function crystalSystemOf(symbol) {
  if (!symbol) return null;
  return CATALOGUE.pointGroups.find((group) => group.symbol === symbol)?.crystal_system ?? null;
}

function normalise(value) {
  if (!value) return CATALOGUE.phases[0] ?? null;
  if (typeof value === 'string') {
    return CATALOGUE.phases.find((entry) => entry.id === value) ?? null;
  }
  if (value.builtin) {
    return CATALOGUE.phases.find((entry) => entry.id === value.builtin) ?? null;
  }
  if (value.cif) return null;
  return value;
}

function cifOf(value) {
  const cif = value && typeof value === 'object' ? value.cif : null;
  if (!cif || typeof cif !== 'object') return null;
  return { name: String(cif.name ?? ''), text: String(cif.text ?? '') };
}
