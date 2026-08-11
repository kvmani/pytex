# Hexagonal and Trigonal Conventions in PyTex

## Purpose

Hexagonal and trigonal systems are a common source of avoidable confusion because internal coordinate storage, direct-space notation, reciprocal-space notation, and teaching conventions are often mixed without warning.

## Canonical Policy

PyTex fixes the following rules:

- Internal direct-space vectors use a three-component direct basis.
- Internal reciprocal-space planes use a three-component reciprocal basis.
- Four-index notation is an interface and documentation convention, not an internal storage requirement.

## Direction Mapping

For three-index direct-basis directions $[u\,v\,w]$, the four-index Weber form follows

\begin{align*}
U &= \frac{2u - v}{3}, \\
  V &= \frac{2v - u}{3}, \\
  T &= -\frac{u + v}{3}, \\
  W &= w,
\end{align*}

with $U + V + T = 0$.

## Plane Mapping

For reciprocal-space planes, the Miller–Bravais redundancy is

$$
i = -(h + k)
$$

so the four-index plane $(h\,k\,i\,l)$ is normalized internally to the three-component reciprocal-basis form $(h\,k\,l)$.

## Normative References

- M. De Graef, *Introduction to Conventional Transmission Electron Microscopy*, Cambridge University Press, 2003. DOI: <https://doi.org/10.1017/CBO9780511615092>.
- Th. Hahn (ed.), *International Tables for Crystallography, Volume A: Space-Group Symmetry*, IUCr / Springer. DOI: <https://doi.org/10.1107/97809553602060000100>.
