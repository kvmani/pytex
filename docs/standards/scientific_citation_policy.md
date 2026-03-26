# Scientific Citation Policy

PyTex is a scientific software project and should be documented like one.

## Citation Requirement

Major scientific docs must contain explicit citations for:

- crystallographic conventions
- mathematical formulations
- algorithm families
- standard data-exchange assumptions
- comparison or parity claims against external tools

Web-facing concept, tutorial, and workflow pages should also cite scientific sources when they explain conventions, mathematics, or algorithm choices, rather than leaving citations only to the PDF layer.

## Source Hierarchy

Use this hierarchy when fixing conventions or semantics:

1. IUCr and International Tables for Crystallography
2. formal data standards and dictionaries such as CIF-related standards
3. canonical textbooks used in the field
4. peer-reviewed papers
5. maintained project documentation for tools such as MTEX or ORIX
6. vendor or educational notes

If PyTex follows a tool convention rather than a formal standard, the docs must say so explicitly.

## Normative Vs Informative References

- `Normative`: a source that fixes behavior, notation, or semantics that PyTex intentionally adopts
- `Informative`: a source that explains, motivates, or compares behavior without fixing it

Every major scientific note should separate these categories.

## Repository Rule

Do not introduce a crystallographic convention, algorithm description, or notation mapping into the stable docs without at least one explicit reference.
