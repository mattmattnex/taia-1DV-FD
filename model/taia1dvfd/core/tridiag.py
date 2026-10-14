import torch

def _shift(P: torch.Tensor, s: int, pad: torch.Tensor, up: bool) -> torch.Tensor:
    """Row i of the result is row i-s (up=False) or i+s (up=True) of P, padded with `pad`."""
    if up:
        return torch.cat((P[..., s:, :], pad[..., :s, :]), dim=-2)
    return torch.cat((pad[..., :s, :], P[..., :-s, :]), dim=-2)


def pcr(a: torch.Tensor, b: torch.Tensor, c: torch.Tensor, d: torch.Tensor) -> torch.Tensor:
    """
    Batched parallel cyclic reduction for tridiagonal systems  L x = d.

      a : (*batch, n-1)  sub-diagonal    a[i] = L[i+1, i]
      b : (*batch, n)    diagonal        b[i] = L[i, i]
      c : (*batch, n-1)  super-diagonal  c[i] = L[i, i+1]
      d : (*batch, n, k) right-hand side

    ceil(log2 n) fully vectorized steps instead of a loop over n: far fewer
    (larger) torch ops, so much less kernel-launch / dispatch overhead.
    Costs O(n log n) flops instead of O(n), which is irrelevant when the
    solve is overhead-bound. Same stability as Thomas (diagonal dominance).
    """
    n, k = d.shape[-2], d.shape[-1]
    dtype = torch.promote_types(torch.promote_types(a.dtype, b.dtype), d.dtype)
    batch = torch.broadcast_shapes(a.shape[:-1], b.shape[:-1], c.shape[:-1], d.shape[:-2])

    # Full-length diagonals: a[0] = 0, c[n-1] = 0. Pack [a, b, c, d] -> (*batch, n, 3+k)
    z = torch.zeros(*batch, 1, dtype=dtype, device=d.device)
    a = torch.cat((z, a.to(dtype).expand(*batch, n - 1)), dim=-1)
    c = torch.cat((c.to(dtype).expand(*batch, n - 1), z), dim=-1)
    b = b.to(dtype).expand(*batch, n)
    P = torch.cat((a.unsqueeze(-1), b.unsqueeze(-1), c.unsqueeze(-1),
                   d.to(dtype).expand(*batch, n, k)), dim=-1)

    # Padding rows are the trivial equation 1 * x = 0 (only b = 1)
    pad = torch.zeros(*batch, n, 3 + k, dtype=dtype, device=d.device)
    pad[..., 1] = 1

    s = 1
    while s < n:
        Pm = _shift(P, s, pad, up=False)    # rows i-s
        Pp = _shift(P, s, pad, up=True)     # rows i+s

        a,  b,  c,  d  = P[..., 0:1],  P[..., 1:2],  P[..., 2:3],  P[..., 3:]
        am, bm, cm, dm = Pm[..., 0:1], Pm[..., 1:2], Pm[..., 2:3], Pm[..., 3:]
        ap, bp, cp, dp = Pp[..., 0:1], Pp[..., 1:2], Pp[..., 2:3], Pp[..., 3:]

        k1 = a / bm                          # eliminates x_{i-s}
        k2 = c / bp                          # eliminates x_{i+s}
        P = torch.cat((-k1 * am,
                       b - k1 * cm - k2 * ap,
                       -k2 * cp,
                       d - k1 * dm - k2 * dp), dim=-1)
        s *= 2

    return P[..., 3:] / P[..., 1:2]


def pcr_solve(Lop: torch.Tensor, rhs: torch.Tensor) -> torch.Tensor:
    """Drop-in replacement for torch.linalg.solve(Lop, rhs) when Lop is tridiagonal (PCR).
    Kept for the equations not yet converted to banded form."""
    a = torch.diagonal(Lop, offset=-1, dim1=-2, dim2=-1)
    b = torch.diagonal(Lop, offset=0,  dim1=-2, dim2=-1)
    c = torch.diagonal(Lop, offset=1,  dim1=-2, dim2=-1)
    return pcr(a, b, c, rhs)


# =====================================================================
#  Banded helpers (no dense (ne, n, n) matrices)
# =====================================================================
def tridiag_apply(lo: torch.Tensor, dg: torch.Tensor, up: torch.Tensor,
                  x: torch.Tensor) -> torch.Tensor:
    """
    y = T x for a tridiagonal T given by its diagonals (same convention as pcr):
      lo : (ne, n-1)  lo[i] = T[i+1, i]
      dg : (ne, n)    dg[i] = T[i, i]
      up : (ne, n-1)  up[i] = T[i, i+1]
      x  : (ne, n)
    """
    y = dg * x
    y[:, :-1] = y[:, :-1] + up * x[:, 1:]
    y[:, 1:]  = y[:, 1:]  + lo * x[:, :-1]
    return y


def ddz(x: torch.Tensor, dz) -> torch.Tensor:
    """
    True d/dz of x (ne, n) on a grid whose index increases downward (z[i+1] = z[i] - dz).
    Centred in the interior, second-order one-sided at both ends.
    Same result as the old dense operators d_dz_c / d_dz_w.
    """
    d = torch.empty_like(x)
    d[:, 1:-1] = (x[:, :-2] - x[:, 2:]) / (2 * dz)
    d[:, 0]    = ( 1.5 * x[:, 0]  - 2 * x[:, 1]  + 0.5 * x[:, 2])  / dz
    d[:, -1]   = (-1.5 * x[:, -1] + 2 * x[:, -2] - 0.5 * x[:, -3]) / dz
    return d