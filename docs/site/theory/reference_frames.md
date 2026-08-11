# Reference Frames and Conventions in PyTex

## Frame Domains

PyTex distinguishes at least the following frame domains:

- crystal frame,
- specimen frame,
- map frame,
- detector frame,
- laboratory frame,
- reciprocal frame.

## Transform Semantics

A frame transform is represented as a proper rigid transform

$$
\mathbf{x}_{\text{target}} = \mathbf{R}\mathbf{x}_{\text{source}} + \mathbf{t}
$$

with $\mathbf{R} \in SO(3)$ and explicit source and target frames.

## Why This Matters

The same numerical vector can mean different things in different software environments. PyTex aims to force those meanings to remain inspectable.
