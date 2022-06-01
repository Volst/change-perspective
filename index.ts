export type QuadPoints = [
  number,
  number,
  number,
  number,
  number,
  number,
  number,
  number
];

export type SinglePoint = [number, number];

function dim(x: number[][]) {
  let y: number[];
  if (typeof x === 'object') {
    y = x[0];
    if (typeof y === 'object') {
      return [x.length, y.length];
    }
    return [x.length];
  }
  return [];
}

function _foreach2(x: any, s: number[], k: number, f: (x: any) => void) {
  if (k === s.length - 1) {
    return f(x);
  }
  let i;
  const n = s[k];
  const ret = Array(n);
  for (i = n - 1; i >= 0; i--) {
    ret[i] = _foreach2(x[i], s, k + 1, f);
  }
  return ret;
}

function cloneV(x: number[]) {
  const _n = x.length;
  let i: number;
  const ret = Array(_n);

  for (i = _n - 1; i !== -1; --i) {
    ret[i] = x[i];
  }
  return ret;
}

function clone(x: any) {
  if (typeof x !== 'object') return x;
  const V = cloneV;
  const s = dim(x);
  return _foreach2(x, s, 0, V);
}

function diag(d: number[]) {
  let i: number;
  let i1: number;
  let j: number;
  const n = d.length;
  const A = Array(n);
  let Ai;
  for (i = n - 1; i >= 0; i--) {
    Ai = Array(n);
    i1 = i + 2;
    for (j = n - 1; j >= i1; j -= 2) {
      Ai[j] = 0;
      Ai[j - 1] = 0;
    }
    if (j > i) {
      Ai[j] = 0;
    }
    Ai[i] = d[i];
    for (j = i - 1; j >= 1; j -= 2) {
      Ai[j] = 0;
      Ai[j - 1] = 0;
    }
    if (j === 0) {
      Ai[0] = 0;
    }
    A[i] = Ai;
  }
  return A;
}

function rep(s: number[], v: number, k?: number) {
  if (typeof k === 'undefined') {
    k = 0;
  }
  const n = s[k];
  const ret = Array(n);
  let i;
  if (k === s.length - 1) {
    for (i = n - 2; i >= 0; i -= 2) {
      ret[i + 1] = v;
      ret[i] = v;
    }
    if (i === -1) {
      ret[0] = v;
    }
    return ret;
  }
  for (i = n - 1; i >= 0; i--) {
    ret[i] = rep(s, v, k + 1);
  }
  return ret;
}

function identity(n: number) {
  return diag(rep([n], 1));
}

function inv(a: number[][]) {
  const s = dim(a);
  const abs = Math.abs;
  const m = s[0];
  const n = s[1];
  const A = clone(a);
  let Ai: number[];
  let Aj: number[];
  const I = identity(m);
  let Ii;
  let Ij;
  var i, j, k, x;
  for (j = 0; j < n; ++j) {
    var i0 = -1;
    var v0 = -1;
    for (i = j; i !== m; ++i) {
      k = abs(A[i][j]);
      if (k > v0) {
        i0 = i;
        v0 = k;
      }
    }
    Aj = A[i0];
    A[i0] = A[j];
    A[j] = Aj;
    Ij = I[i0];
    I[i0] = I[j];
    I[j] = Ij;
    x = Aj[j];
    for (k = j; k !== n; ++k) Aj[k] /= x;
    for (k = n - 1; k !== -1; --k) Ij[k] /= x;
    for (i = m - 1; i !== -1; --i) {
      if (i !== j) {
        Ai = A[i];
        Ii = I[i];
        x = Ai[j];
        for (k = j + 1; k !== n; ++k) Ai[k] -= Aj[k] * x;
        for (k = n - 1; k > 0; --k) {
          Ii[k] -= Ij[k] * x;
          --k;
          Ii[k] -= Ij[k] * x;
        }
        if (k === 0) Ii[0] -= Ij[0] * x;
      }
    }
  }
  return I;
}

function dotMMsmall(x: number[][], y: number[][]) {
  var i, j, k, p, q, r, ret, foo, bar, woo, i0;
  p = x.length;
  q = y.length;
  r = y[0].length;
  ret = Array(p);
  for (i = p - 1; i >= 0; i--) {
    foo = Array(r);
    bar = x[i];
    for (k = r - 1; k >= 0; k--) {
      woo = bar[q - 1] * y[q - 1][k];
      for (j = q - 2; j >= 1; j -= 2) {
        i0 = j - 1;
        woo += bar[j] * y[j][k] + bar[i0] * y[i0][k];
      }
      if (j === 0) {
        woo += bar[0] * y[0][k];
      }
      foo[k] = woo;
    }
    ret[i] = foo;
  }
  return ret;
}

function dotMV(x: number[][], y: number[]) {
  const p = x.length;
  let i: number;
  const ret = Array(p);
  for (i = p - 1; i >= 0; i--) {
    ret[i] = dotVV(x[i], y);
  }
  return ret;
}

function dotVV(x: number[], y: number[]) {
  let i;
  const n = x.length;
  let i1;
  let ret = x[n - 1] * y[n - 1];
  for (i = n - 2; i >= 1; i -= 2) {
    i1 = i - 1;
    ret += x[i] * y[i] + x[i1] * y[i1];
  }
  if (i === 0) {
    ret += x[0] * y[0];
  }
  return ret;
}

function transpose(x: number[][]) {
  let i;
  let j;
  const m = x.length;
  const n = x[0].length;
  const ret = Array(n);
  let A0;
  let A1;
  let Bj;
  for (j = 0; j < n; j++) ret[j] = Array(m);
  for (i = m - 1; i >= 1; i -= 2) {
    A1 = x[i];
    A0 = x[i - 1];
    for (j = n - 1; j >= 1; --j) {
      Bj = ret[j];
      Bj[i] = A1[j];
      Bj[i - 1] = A0[j];
      --j;
      Bj = ret[j];
      Bj[i] = A1[j];
      Bj[i - 1] = A0[j];
    }
    if (j === 0) {
      Bj = ret[0];
      Bj[i] = A1[0];
      Bj[i - 1] = A0[0];
    }
  }
  if (i === 0) {
    A0 = x[0];
    for (j = n - 1; j >= 1; --j) {
      ret[j][0] = A0[j];
      --j;
      ret[j][0] = A0[j];
    }
    if (j === 0) {
      ret[0][0] = A0[0];
    }
  }
  return ret;
}

function round(num: number) {
  return Math.round(num * 10000000000) / 10000000000;
}

export function getNormalizationCoefficients(
  srcPts: QuadPoints,
  dstPts: QuadPoints,
  isInverse?: boolean
): QuadPoints {
  if (isInverse) {
    const tmp = dstPts;
    dstPts = srcPts;
    srcPts = tmp;
  }
  const r1 = [
    srcPts[0],
    srcPts[1],
    1,
    0,
    0,
    0,
    -1 * dstPts[0] * srcPts[0],
    -1 * dstPts[0] * srcPts[1],
  ];
  const r2 = [
    0,
    0,
    0,
    srcPts[0],
    srcPts[1],
    1,
    -1 * dstPts[1] * srcPts[0],
    -1 * dstPts[1] * srcPts[1],
  ];
  const r3 = [
    srcPts[2],
    srcPts[3],
    1,
    0,
    0,
    0,
    -1 * dstPts[2] * srcPts[2],
    -1 * dstPts[2] * srcPts[3],
  ];
  const r4 = [
    0,
    0,
    0,
    srcPts[2],
    srcPts[3],
    1,
    -1 * dstPts[3] * srcPts[2],
    -1 * dstPts[3] * srcPts[3],
  ];
  const r5 = [
    srcPts[4],
    srcPts[5],
    1,
    0,
    0,
    0,
    -1 * dstPts[4] * srcPts[4],
    -1 * dstPts[4] * srcPts[5],
  ];
  const r6 = [
    0,
    0,
    0,
    srcPts[4],
    srcPts[5],
    1,
    -1 * dstPts[5] * srcPts[4],
    -1 * dstPts[5] * srcPts[5],
  ];
  const r7 = [
    srcPts[6],
    srcPts[7],
    1,
    0,
    0,
    0,
    -1 * dstPts[6] * srcPts[6],
    -1 * dstPts[6] * srcPts[7],
  ];
  const r8 = [
    0,
    0,
    0,
    srcPts[6],
    srcPts[7],
    1,
    -1 * dstPts[7] * srcPts[6],
    -1 * dstPts[7] * srcPts[7],
  ];

  const matA = [r1, r2, r3, r4, r5, r6, r7, r8];
  const matB = dstPts;
  let matC;

  try {
    matC = inv(dotMMsmall(transpose(matA), matA));
  } catch (e) {
    return [1, 0, 0, 0, 1, 0, 0, 0];
  }

  const matD = dotMMsmall(matC, transpose(matA));
  const matX = dotMV(matD, matB);
  for (var i = 0; i < matX.length; i++) {
    matX[i] = round(matX[i]);
  }
  matX[8] = 1;

  return matX as QuadPoints;
}

function applyTransform(coeffs: QuadPoints, x: number, y: number): SinglePoint {
  const coordinates = [];
  coordinates[0] =
    (coeffs[0] * x + coeffs[1] * y + coeffs[2]) /
    (coeffs[6] * x + coeffs[7] * y + 1);
  coordinates[1] =
    (coeffs[3] * x + coeffs[4] * y + coeffs[5]) /
    (coeffs[6] * x + coeffs[7] * y + 1);
  return coordinates as SinglePoint;
}

export default function fixPerspective(srcPts: QuadPoints, dstPts: QuadPoints) {
  const coeffs = getNormalizationCoefficients(srcPts, dstPts, false);
  return (x: number, y: number) => applyTransform(coeffs, x, y);
}
