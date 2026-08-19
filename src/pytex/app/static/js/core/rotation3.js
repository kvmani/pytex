/**
 * Row-major 3x3 rotation arithmetic, shared by everything that holds a camera.
 *
 * These nine-number matrices are the browser's whole share of orientation work.
 * Every crystallographic decision — which poles exist, which symmetry operators
 * apply, where the fundamental sector is — is made in Python and arrives as
 * finished vectors; what happens here is multiplication.
 *
 * The functions live in one module because the crystal viewer and its
 * orientation figures must compose *the same* camera. Two copies of `multiply`
 * that disagree in transpose convention would put the pole figure a transpose
 * away from the structure it claims to describe, and the picture would look
 * plausible the whole time.
 */

/** The identity rotation. */
export function identity() {
  return [1, 0, 0, 0, 1, 0, 0, 0, 1];
}

/** Matrix product `a b`, both row-major. */
export function multiply(a, b) {
  const out = new Array(9).fill(0);
  for (let row = 0; row < 3; row += 1) {
    for (let column = 0; column < 3; column += 1) {
      let sum = 0;
      for (let k = 0; k < 3; k += 1) sum += a[row * 3 + k] * b[k * 3 + column];
      out[row * 3 + column] = sum;
    }
  }
  return out;
}

/** Rotation by `angle` radians about the x axis. */
export function rotationX(angle) {
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  return [1, 0, 0, 0, c, -s, 0, s, c];
}

/** Rotation by `angle` radians about the y axis. */
export function rotationY(angle) {
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  return [c, 0, s, 0, 1, 0, -s, 0, c];
}

/** `m v`. */
export function applyMatrix(m, v) {
  return [
    m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
    m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
    m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
  ];
}

/** `mᵀ v` — for a rotation, the inverse map. */
export function applyTranspose(m, v) {
  return [
    m[0] * v[0] + m[3] * v[1] + m[6] * v[2],
    m[1] * v[0] + m[4] * v[1] + m[7] * v[2],
    m[2] * v[0] + m[5] * v[1] + m[8] * v[2],
  ];
}

/** Vector cross product. */
export function cross(a, b) {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

/** Unit vector in the direction of `v`; the zero vector is returned unchanged in direction. */
export function normalise(v) {
  const length = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / length, v[1] / length, v[2] / length];
}

/**
 * A camera looking down a crystal direction.
 *
 * The third row of a camera is the viewing direction in crystal coordinates, so
 * "look along c" is the rotation whose third row is c. The other two rows are
 * any orthonormal pair completing it; the reference vector is switched near the
 * poles so the cross product never collapses.
 */
export function lookAlong(vector) {
  const forward = normalise(vector);
  const reference = Math.abs(forward[1]) > 0.95 ? [0, 0, 1] : [0, 1, 0];
  const right = normalise(cross(reference, forward));
  const up = cross(forward, right);
  return [right[0], right[1], right[2], up[0], up[1], up[2], forward[0], forward[1], forward[2]];
}
