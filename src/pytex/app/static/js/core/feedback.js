/**
 * The feedback and feature-request form.
 *
 * What it is for
 * --------------
 * PyTex is written by two or three people for a community much larger than
 * that, and almost everything worth adding next is something a user noticed
 * and nobody asked them about. This is the place they can say it without
 * leaving the workbench, finding an address, or knowing who to write to.
 *
 * The tone is deliberate. A form headed "Report a bug" collects bugs; a form
 * that says it wants to hear what would save you an afternoon collects the
 * thing that is actually scarce. The invitation text comes from the server, so
 * a group running its own copy can name its own maintainer.
 *
 * What it asks, and what it does not
 * ----------------------------------
 * Only the message is required. Name, e-mail and group are optional on
 * purpose: a note worth reading is worth reading anonymously, and requiring a
 * name suppresses exactly the frank criticism that is most useful.
 *
 * The page attaches *context* — which workspace and panel were open — because
 * "the legend is unreadable" is a different report depending on which figure
 * it was about, and nobody should have to remember to say. That context is
 * recorded as claimed rather than trusted; the server stamps what it knows
 * itself.
 */

import { sendFeedback } from './api.js';
import { clear, el, markdown } from './dom.js';
import * as log from './logbook.js';

/** Where the draft is kept between openings, so a half-written note survives. */
const DRAFT_KEY = 'pytex-feedback-draft';

/**
 * Mount the feedback drawer.
 *
 * @param {object} options
 * @param {HTMLElement} options.drawer - The drawer element to fill.
 * @param {HTMLElement} options.body - The drawer's body.
 * @param {object} options.config - The `feedback` block from `/api/experience`.
 * @param {Function} options.context - Called at submit time for `{panel, workspace}`.
 * @returns {{open: Function, close: Function}}
 */
export function mountFeedback({ drawer, body, config, context = () => ({}) }) {
  let state = 'form';

  function close() {
    drawer.hidden = true;
  }

  function open() {
    state = 'form';
    render();
    drawer.hidden = false;
    body.querySelector('textarea')?.focus();
    log.info('Opened the feedback form.', { source: 'feedback' });
  }

  function render() {
    clear(body);
    if (state === 'sent') return;
    body.append(...formNodes());
  }

  function formNodes() {
    const categories = config.categories ?? [{ value: 'feedback', label: 'General feedback' }];
    const draft = readDraft();

    const category = el(
      'select',
      { id: 'feedback-category' },
      categories.map((entry) => el('option', { value: entry.value, text: entry.label })),
    );
    if (draft.category) category.value = draft.category;

    const message = el('textarea', {
      id: 'feedback-message',
      rows: 7,
      // Deliberately not `required`: the browser's own bubble would pre-empt
      // the form's submit handler, and "Please fill out this field" is a worse
      // thing to say to somebody about to give you their time than the
      // sentence this form says instead.
      maxlength: config.max_message_characters ?? 8000,
      placeholder:
        'What would you change, add, or explain differently? One line is plenty — ' +
        'say what you expected and what happened instead.',
      oninput: () => {
        saveDraft();
        counter.textContent = remaining();
        error.hidden = true;
      },
    });
    message.value = draft.message ?? '';

    const name = shortField('feedback-name', 'Your name', draft.name, 'Entirely optional');
    const email = shortField(
      'feedback-email',
      'E-mail',
      draft.email,
      'Only if you would like an answer',
      'email',
    );
    const organisation = shortField(
      'feedback-organisation',
      'Group or laboratory',
      draft.organisation,
      'Helps us know who PyTex is reaching',
    );

    const consent = el('input', { type: 'checkbox', id: 'feedback-consent' });
    consent.checked = draft.contact_consent ?? true;

    const rating = el(
      'select',
      { id: 'feedback-rating' },
      [
        el('option', { value: '', text: 'Rather not say' }),
        ...[5, 4, 3, 2, 1].map((value) =>
          el('option', { value: String(value), text: `${value} — ${RATING_WORDS[value]}` }),
        ),
      ],
    );
    if (draft.rating) rating.value = String(draft.rating);

    const counter = el('span.feedback__counter', { text: '' });
    const error = el('p.field__error', { hidden: true });
    const submit = el('button.button.button--primary', { type: 'submit', text: 'Send it' });

    function remaining() {
      const limit = config.max_message_characters ?? 8000;
      const left = limit - message.value.length;
      return left > limit / 4 ? '' : `${left} characters left`;
    }

    function saveDraft() {
      writeDraft({
        category: category.value,
        message: message.value,
        name: name.input.value,
        email: email.input.value,
        organisation: organisation.input.value,
        rating: rating.value,
        contact_consent: consent.checked,
      });
    }

    const form = el(
      'form.feedback',
      {
        onsubmit: async (event) => {
          event.preventDefault();
          if (!message.value.trim()) {
            error.hidden = false;
            error.textContent = 'There is nothing in the note yet — even one line is useful.';
            message.focus();
            return;
          }
          submit.disabled = true;
          submit.textContent = 'Sending…';
          try {
            const payload = await sendFeedback({
              category: category.value,
              message: message.value,
              name: name.input.value,
              email: email.input.value,
              organisation: organisation.input.value,
              rating: rating.value === '' ? null : Number(rating.value),
              contact_consent: consent.checked,
              context: context(),
            });
            clearDraft();
            state = 'sent';
            showReceipt(payload);
          } catch (failure) {
            error.hidden = false;
            error.textContent = [failure.message, failure.hint].filter(Boolean).join(' ');
            log.error(`Feedback could not be sent: ${failure.message}`, { source: 'feedback' });
          } finally {
            submit.disabled = false;
            submit.textContent = 'Send it';
          }
        },
      },
      [
        el('h2', { text: 'Tell us what would make PyTex better' }),
        ...markdown(config.invitation ?? ''),
        el('p.feedback__privacy', {
          text: config.relayed
            ? 'Your note is filed with this installation and e-mailed to the maintainer.'
            : 'Your note is filed with this installation, where the person who runs it will read it.',
        }),

        field('feedback-category', 'What is this about?', category),
        field('feedback-message', 'Your note', message),
        counter,

        el('details.feedback__about-you', { open: true }, [
          el('summary', { text: 'About you (all optional)' }),
          el('div.feedback__grid', {}, [name.element, email.element, organisation.element]),
          field('feedback-rating', 'How is PyTex working out so far?', rating),
          el('label.checkbox', {}, [
            consent,
            el('span', { text: 'You may write back to me about this' }),
          ]),
        ]),

        error,
        el('div.feedback__actions', {}, [
          submit,
          el('button.button', { type: 'button', text: 'Not now', onclick: close }),
        ]),
      ],
    );

    counter.textContent = remaining();
    for (const control of [category, rating, consent, name.input, email.input, organisation.input]) {
      control.addEventListener('change', saveDraft);
    }
    return [form];
  }

  /**
   * What became of the note, said plainly.
   *
   * The receipt distinguishes "filed here" from "filed and e-mailed", and it
   * says so rather than showing one cheerful message for both. A user in a
   * deployment with no relay is entitled to know that what happens next
   * depends on somebody reading a local file.
   */
  function showReceipt(payload) {
    const receipt = payload.receipt ?? {};
    clear(body);
    body.append(
      el('div.feedback__receipt', {}, [
        el('h2', { text: receipt.stored ? 'Thank you — it is on its way' : 'That did not save' }),
        ...markdown(payload.acknowledgement ?? ''),
        el('p.field__help', {
          text: receipt.delivery_detail ?? '',
        }),
        el('div.feedback__actions', {}, [
          el('button.button.button--primary', {
            type: 'button',
            text: 'Close',
            onclick: close,
          }),
          el('button.button', {
            type: 'button',
            text: 'Send another',
            onclick: () => {
              state = 'form';
              render();
              body.querySelector('textarea')?.focus();
            },
          }),
        ]),
      ]),
    );
    log.notice(
      receipt.stored
        ? `Feedback recorded${receipt.delivered ? ' and e-mailed' : ''}.`
        : 'Feedback could not be recorded.',
      { source: 'feedback' },
    );
  }

  return { open, close };
}

const RATING_WORDS = {
  5: 'it is doing what I need',
  4: 'good, with rough edges',
  3: 'mixed',
  2: 'more friction than help',
  1: 'it is getting in my way',
};

function field(id, label, control) {
  return el('div.field', {}, [el('label.field__label', { for: id }, [label]), control]);
}

function shortField(id, label, value, hint, type = 'text') {
  const input = el('input', { id, type, autocomplete: 'off' });
  input.value = value ?? '';
  const element = el('div.field', {}, [
    el('label.field__label', { for: id }, [label]),
    input,
    hint ? el('p.field__help.field__help--inline', { text: hint }) : null,
  ]);
  return { element, input };
}

/* ------------------------------------------------------------------ drafts */

/**
 * A half-written note survives closing the drawer, and a browser crash.
 *
 * People start a note, go back to check the thing they are complaining about,
 * and come back. Losing the text at that point loses the report — they do not
 * type it twice.
 */
function readDraft() {
  try {
    return JSON.parse(localStorage.getItem(DRAFT_KEY) ?? '{}') ?? {};
  } catch {
    return {};
  }
}

function writeDraft(draft) {
  try {
    localStorage.setItem(DRAFT_KEY, JSON.stringify(draft));
  } catch {
    // A locked-down webview may disable storage. The form still works; only
    // the convenience of a surviving draft is lost.
  }
}

function clearDraft() {
  try {
    localStorage.removeItem(DRAFT_KEY);
  } catch {
    // As above.
  }
}
