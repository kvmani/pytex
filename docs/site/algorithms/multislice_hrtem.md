# Multislice HRTEM Simulation

**Surface:** `pytex.diffraction.multislice` — `multislice`, `MultisliceExitWave`
(`image`, `focal_series`, `defocus_thickness_map`, `diffraction_pattern`, `beam_intensities`),
`hrtem_image`, `simulate_multislice_hrem`, `periodic_slab`, `zone_axis_cell`, `SlicedPotential`,
`MultisliceGrid` — with `pytex.adapters.abtem.simulate_hrem(engine=...)` as the unified entry point
and the workbench operations `tem.simulate_hrem` and `tem.hrtem_series`.

This page states how an HRTEM image is computed from atoms: the steps a reader could reimplement,
the constraints on each, what it costs, and how it fails. The derivations are in
{doc}`../theory/multislice_hrtem`; the objective-lens optics in
{doc}`../theory/hrem_multislice_and_ctf`; the verified numbers in
{doc}`the worked examples </examples/generated/multislice-hrtem>`; and a guided tour of the cases in
{doc}`the multislice tutorial </tutorials/notebooks/36_multislice_hrtem>`.

```{figure} ../../figures/multislice_hrtem_algorithm.svg
:alt: Three-lane flow sheet. Lane 1 builds an exactly periodic orthogonal slab with the zone axis
  along the beam, cuts it into slices, and forms each slice's projected potential from
  parametrized electron scattering factors. Lane 2 alternates band-limited transmission with
  Fresnel propagation and stores exit waves at the requested depths. Lane 3 applies the objective
  lens with coherence envelopes, averages intensities over frozen phonons, and repeats only the
  lens step for every image of a focal or thickness series.
:width: 100%

The multislice pipeline. Only lane 3 is repeated for a focal series.
```

## 1. The algorithm

**Input.** An `AtomicSnapshot` in an orthogonal box $[0, L_x)\times[0, L_y)\times[0, L_z)$ with the
beam along $+z$; the beam energy; a pixel size $\Delta x$ (or an explicit `MultisliceGrid`); a slice
thickness $\Delta z$; a scattering-factor parametrization; optionally a beam tilt, exit depths,
a static Debye–Waller displacement and a frozen-phonon ensemble.

1. **Specimen.** For a crystal, `periodic_slab(phase, [uvw], (m1, m2), thickness_angstrom=t)`
   builds the box from `zone_axis_cell`: $\mathbf t_3 \parallel [uvw]$ primitive,
   $\mathbf t_1 \perp \mathbf t_3$ shortest, $\mathbf t_2 \perp$ both shortest, from an integer
   search up to index 6. It tiles $m_1\times m_2$ laterally and
   $m_3 = \lceil t/|\mathbf t_3|\rceil$ along the beam, keeps atoms with fractional box coordinates
   in $[0,1)$, and checks the count.
2. **Grid.** $N_x = \mathrm{next\_fast\_len}(\lceil L_x/\Delta x\rceil)$, likewise $N_y$, so the
   delivered pixel is never coarser than requested. The band limit is
   $g_{\max} = \tfrac23\cdot\tfrac{1}{2\max(\Delta x,\Delta y)}$, with a raised-cosine edge of width
   $0.01/\max(\Delta x, \Delta y)$.
3. **Slicing.** $n = \lceil L_z/\Delta z\rceil$ equal slices of $L_z/n$. Each atom goes to the slice
   holding its centre (with a $10^{-6}$ tolerance so an atom on a boundary goes to the slice that
   begins there). Each slice gets a key — a hash of its sorted $(Z, x, y)$ — so slices with the same
   projection are recognised.
4. **Projected potential**, per distinct slice: for each species, the structure factor
   $S(\mathbf g) = \sum_j e^{-2\pi i g_x x_j}\,e^{-2\pi i g_y y_j}$ as one
   $(N_y\times N_{\mathrm{at}})(N_{\mathrm{at}}\times N_x)$ product; multiply by $f_e(g)$ (and the
   Debye–Waller factor), sum over species, scale by $47.878/A$, and inverse FFT.
5. **Transmission**: $\tau_n = \mathcal F^{-1}\{\mathcal F\{e^{i\sigma_e v_n}\}\cdot
   A_{\mathrm{aa}}\}$. Transmission functions of repeated slices are cached.
6. **Propagation**: starting from $\psi_0 = 1$, for $n = 0 \dots N-1$:
   $\psi \leftarrow \mathcal F^{-1}\{\mathcal F\{\psi\,\tau_n\}\cdot\mathcal P_{\Delta z}\}$, with
   $\mathcal P_{\Delta z}$ band-limited and carrying the tilt. Transmit, then propagate — abTEM's
   order — so the exit plane is the bottom of the box. After slice $n$, if its lower edge is the
   slice boundary nearest a requested depth, store $\psi$.
7. **Frozen phonons** (optional): repeat 3–6 for each of $N_{\mathrm{fp}}$ configurations with
   independent Gaussian displacements (seeded), keeping every configuration's exit waves.
8. **Imaging**, for each requested depth and defocus: $H = A\,E_s\,E_c\,e^{-i\chi}$ on the Fourier
   grid; $I = \langle|\mathcal F^{-1}\{\mathcal F\{\psi^{(c)}\}H\}|^2\rangle_c$. With focal
   integration, $E_c$ is replaced by a Gauss–Hermite sum over defocus whose node count
   $n = \lceil a^2/4\rceil + 16$ (24 to 512) resolves the largest defocus phase rate
   $a = \sqrt2\pi\lambda\Delta g^2$ the aperture and grid pass.
9. **Diagnostics**: the retained intensity $\langle|\psi|^2\rangle$ at each stored depth, the band
   limit in Å⁻¹ and mrad, slice and distinct-slice counts; all in `describe()`.

**Output.** `MultisliceExitWave` holding waves of shape
`(configurations, depths, Ny, Nx)`; from it, images, `FocalSeries`, `DefocusThicknessMap`,
diffraction patterns and beam intensities.

## 2. Constraints and their calibration

| Constraint | Value | Why, and what it is checked against |
| --- | --- | --- |
| Band limit | $\tfrac23$ Nyquist | The largest limit for which $\psi\tau$ cannot alias into the band; abTEM's value; worked example *band limit*. |
| Taper | $0.01/\Delta x$ | abTEM's default; the abTEM parity tests use it. |
| Slice assignment tolerance | $10^{-6}$ (in slice units) | An atom exactly on a boundary must not be assigned by rounding noise; without it the abTEM exit-wave parity fell from $6\times10^{-4}$ to $2\times10^{-2}$. |
| Zone-axis search | index $\le 6$ | Enough for every low-index axis of the cubic and hexagonal systems; refusal beyond it rather than a strained box. |
| Focal-integration nodes | $n = \lceil a^2/4\rceil + 16$, at most 512 | Integrates $e^{iax}$ to $10^{-8}$ for $a \le 2\sqrt{n-16}$, checked numerically; SciPy's `roots_hermite` stays finite where NumPy's `hermgauss` overflows above 370 nodes. Worked example *Frank envelope*. |
| Grid ceiling | $4096^2$ points | Several complex arrays of this size are live at once (256 MB each). |

## 3. Cost

Per configuration, with $N = N_xN_y$ grid points, $N_s$ slices and $N_u$ distinct slices:

- potential: $O(N_u\,N\,N_{\mathrm{at}}/N_u)$ — each atom enters one structure-factor product — plus
  $N_u$ FFTs;
- propagation: $2N_s$ FFTs, $O(N_s N\log N)$ — the dominant term;
- imaging: 2 FFTs per image, or $2n$ with $n$-node focal integration, independent of the specimen
  thickness.

A perfect crystal whose period along the beam is divided by $\Delta z$ has $N_u$ equal to the
slices per period, however thick the foil. A focal series of $F$ images costs one propagation and
$2F$ FFTs. Frozen phonons multiply the propagation cost by $N_{\mathrm{fp}}$ and defeat the slice
cache.

## 4. Failure modes

- **Non-periodic box.** Built by bounding-box padding (as `AtomicSnapshot.from_phase` does) instead
  of `periodic_slab`: seams at the edges, beams off the grid, forbidden reflections. Always build
  crystals with `periodic_slab`; the workbench does.
- **Too coarse a sampling.** High-angle scattering is cut and the retained intensity falls; below
  0.95 `describe()` says so. Refine $\Delta x$.
- **Too thick a slice.** The splitting error grows with $\Delta z$ and with heavy atoms; check by
  halving $\Delta z$.
- **Quasi-coherent envelope on a strong object.** Heavy columns at large focal spread: compare with
  `temporal_coherence="focal_integration"`.
- **Frozen phonons with too few configurations.** The TDS background is noisy; the elastic wave is
  still correct.
- **Low-symmetry zone axis.** `zone_axis_cell` refuses when no orthogonal lattice pair exists
  within the search range.

## 5. Validation

The validation table and its sources are in the theory note,
{ref}`Validation <multislice-validation>`. The independent checks are closed forms (potential
integral, mean inner potential, Kirkland's real-space projection, propagator identities, Frank's
envelope), a second method (Bloch waves on the same potential) and a second code (abTEM, when
installed). The corresponding rows of the {doc}`diffraction validation matrix
<../validation/diffraction_validation_matrix>` point at the tests.

## References

- Kirkland, E. J. (2010). *Advanced Computing in Electron Microscopy*, 2nd ed., ch. 6. Springer.
- Madsen, J. & Susi, T. (2021). The abTEM code: transmission electron microscopy from first
  principles. *Open Research Europe* **1**, 24.
- Cowley, J. M. & Moodie, A. F. (1957). *Acta Cryst.* **10**, 609–619.
