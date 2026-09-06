/**
 * The candidate-phase list, for the one operation that is told several phases
 * rather than one.
 *
 * Every other control on the workbench asks "which phase is this?" and expects
 * an answer. Identification asks "which of these is it?", so the control has to
 * hold a *set* — and the set is the whole scientific input, not a
 * convenience: a ranking against two candidates and a ranking against six are
 * different measurements, and the one thing a user must be able to see at a
 * glance is exactly what was in the running.
 *
 * Three consequences shape the control:
 *
 * **Several CIFs open at once.** A user comparing candidates has a folder of
 * them, and choosing them one file dialogue at a time is a worse experience for
 * no gain, so the file input takes `multiple` and every chosen file becomes a
 * row. Rows are read in the browser and travel as text in the ordinary JSON
 * request, exactly as the single-phase control's CIF does; no crystallography
 * happens here.
 *
 * **A row can be removed but the list cannot be emptied silently.** The list is
 * the question being asked. An empty one is a question with no content, so the
 * control says so in place of the rows rather than submitting nothing and
 * letting the server produce the complaint.
 *
 * **A catalogue row stays a reference.** As in the single-phase control, an
 * unedited catalogue choice emits `{builtin: id}` rather than a copied cell, so
 * a later correction to a published lattice parameter reaches saved analyses.
 * A file row emits `{cif: {name, text}}`.
 *
 * The emitted value is `{phases: [{label, phase}, ...]}`, which is what the
 * `xrd.phase_identification` handler reads.
 */

import { el } from './dom.js';
import { explainer } from './explainer.js';

const CATALOGUE = { phases: [] };

/** Share the phase catalogue, loaded once at start-up alongside the phase control's. */
export function setCandidateCatalogue(phases) {
  CATALOGUE.phases = phases ?? [];
}

/** A row is either a catalogue reference or a file; nothing else is offered. */
function rowsOf(value) {
  const entries = value && typeof value === 'object' ? value.phases : null;
  if (!Array.isArray(entries)) return [];
  return entries
    .map((entry) => {
      const phase = entry && typeof entry === 'object' ? entry.phase : null;
      if (!phase || typeof phase !== 'object') return null;
      if (phase.cif) {
        return {
          kind: 'cif',
          name: String(phase.cif.name ?? ''),
          text: String(phase.cif.text ?? ''),
          label: String(entry.label ?? phase.cif.name ?? ''),
        };
      }
      const builtin = String(phase.builtin ?? '');
      if (!builtin) return null;
      return { kind: 'builtin', builtin, label: String(entry.label ?? '') };
    })
    .filter(Boolean);
}

function catalogueName(id) {
  const found = CATALOGUE.phases.find((entry) => entry.id === id);
  return found ? found.name : id;
}

export function phaseCandidatesControl(parameter, value, onChange, id) {
  const state = { rows: rowsOf(value) };

  const list = el('div.candidates');
  const status = el('p.field__help');

  const picker = el('select', {
    id,
    'aria-label': 'Add a built-in phase to the candidates',
    onchange: (event) => {
      const chosen = event.currentTarget.value;
      if (!chosen) return;
      event.currentTarget.value = '';
      // A duplicate row would be ranked twice and score identically, which
      // reads as a bug rather than as a choice. Adding one already present is
      // therefore a no-op with a word about why.
      if (state.rows.some((row) => row.kind === 'builtin' && row.builtin === chosen)) {
        status.textContent = `${catalogueName(chosen)} is already a candidate.`;
        return;
      }
      state.rows.push({ kind: 'builtin', builtin: chosen, label: '' });
      render();
      onChange();
    },
  }, [
    el('option', { value: '', text: 'Add a built-in phase…' }),
    ...CATALOGUE.phases.map((entry) =>
      el('option', { value: entry.id, text: `${entry.name} · ${entry.point_group}` }),
    ),
  ]);

  const fileInput = el('input', {
    type: 'file',
    accept: '.cif,chemical/x-cif',
    multiple: true,
    'aria-label': 'Open one or more CIF files as candidates',
    onchange: (event) => loadFiles(Array.from(event.currentTarget.files ?? [])),
  });

  async function loadFiles(files) {
    if (!files.length) return;
    status.textContent = `Reading ${files.length} file${files.length === 1 ? '' : 's'}…`;
    const failed = [];
    for (const file of files) {
      try {
        const text = await file.text();
        state.rows.push({ kind: 'cif', name: file.name, text, label: file.name });
      } catch {
        failed.push(file.name);
      }
    }
    fileInput.value = '';
    render();
    status.textContent = failed.length
      ? `Could not read ${failed.join(', ')} in the browser.`
      : `${files.length} candidate${files.length === 1 ? '' : 's'} added. `
        + 'Run the identification to validate and normalize their structures.';
    onChange();
  }

  function remove(index) {
    state.rows.splice(index, 1);
    render();
    onChange();
  }

  function render() {
    if (!state.rows.length) {
      list.replaceChildren(
        el('p.candidates__empty', {
          text: 'No candidates yet. Add built-in phases or open .cif files; two or more '
            + 'make this an identification rather than a check on one phase.',
        }),
      );
      return;
    }
    list.replaceChildren(
      ...state.rows.map((row, index) =>
        el('div.candidates__row', {}, [
          el('span.candidates__name', {
            text: row.kind === 'builtin' ? catalogueName(row.builtin) : row.name,
            title: row.kind === 'builtin' ? 'Built-in catalogue phase' : row.name,
          }),
          el('span.candidates__origin', {
            text: row.kind === 'builtin' ? 'catalogue' : 'CIF',
          }),
          el('button.button.candidates__remove', {
            type: 'button',
            text: 'Remove',
            'aria-label': `Remove ${row.kind === 'builtin' ? catalogueName(row.builtin) : row.name}`,
            onclick: () => remove(index),
          }),
        ]),
      ),
    );
  }

  render();

  const element = el('div.candidates__control', {}, [
    list,
    el('div.candidates__actions', {}, [picker, fileInput]),
    status,
    explainer(
      'Two or more candidates make this a choice. One is a check on that phase, which is a '
        + 'useful thing to run but is not an identification, and the result says so.',
      { label: 'How many candidates' },
    ),
  ]);

  return {
    id,
    element,
    read: () => ({
      phases: state.rows.map((row) =>
        row.kind === 'builtin'
          ? { label: catalogueName(row.builtin), phase: { builtin: row.builtin } }
          : { label: row.label || row.name, phase: { cif: { name: row.name, text: row.text } } },
      ),
    }),
    write: (next) => {
      state.rows = rowsOf(next);
      render();
    },
  };
}
