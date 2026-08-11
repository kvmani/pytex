# Euler Convention Handling

PyTex now exposes explicit Euler convention handling for the MTEX-relevant conventions used in the parity suite.

## Supported Conventions

- Bunge $(\phi_1,\Phi,\phi_2)$ with a $Z X Z$ axis sequence
- Matthies $(\alpha,\beta,\gamma)$ with a $Z Y Z$ axis sequence
- ABG as an alias for the Matthies $Z Y Z$ convention

## Policy

PyTex keeps Bunge as the canonical convenience API while also exposing a convention-aware Euler conversion entrypoint so parity tests can target MTEX-facing semantics directly.

## Current Limits

- Only the conventions required for the current parity scope are exposed.
- The parity policy is documented against MTEX 6.1.1 conventions, not every historical toolbox variant.
