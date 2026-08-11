# Diffraction Geometry And Bragg Rings

PyTex now includes a minimal diffraction foundation that evaluates detector coordinates, outgoing directions, scattering vectors, and Bragg ring radii from explicit geometry.

## Detector Coordinates

Detector coordinates are centered at the explicit pattern center and converted from pixels to millimeters using the detector pixel size. Zero tilt aligns the detector normal with the beam direction.

## Outgoing Directions

Given a detector point $\mathbf{r}$ in the laboratory frame, the outgoing unit direction is

$$
\hat{\mathbf{k}}_{\mathrm{out}} = \frac{\mathbf{r}}{\lVert \mathbf{r} \rVert}
$$

## Scattering Vector

With electron wavelength $\lambda$, the incident and outgoing wavevectors are

$$
\mathbf{k}_{\mathrm{in}} = \frac{\hat{\mathbf{k}}_{\mathrm{in}}}{\lambda},
  \qquad
  \mathbf{k}_{\mathrm{out}} = \frac{\hat{\mathbf{k}}_{\mathrm{out}}}{\lambda}
$$

and the scattering vector is

$$
\mathbf{q} = \mathbf{k}_{\mathrm{out}} - \mathbf{k}_{\mathrm{in}}
$$

## Bragg Angle And Ring Radius

For spacing $d$, PyTex currently uses the geometric Bragg relation

$$
2 d \sin \theta = \lambda
$$

with detector ring radius

$$
R = L \tan(2\theta)
$$

where $L$ is the camera length.

## Current Limits

- The current implementation is a geometry foundation, not a full diffraction simulation stack.
- Detector tilt is represented as an intrinsic detector-axis rotation in $(u,v,n)$ order.
- Intensity modeling, calibration refinement, and adapter-backed simulation workflows remain ahead.
