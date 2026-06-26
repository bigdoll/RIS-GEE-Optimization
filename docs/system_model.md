# System model & notation

Uplink of a multi-user RIS-aided network (paper Sec. II).

## Entities

| Symbol | Code | Meaning |
|--------|------|---------|
| `K` | `K` | number of single-antenna users |
| `N_R` | `NR` | number of BS receive antennas |
| `N` | `N` | number of RIS reflecting elements |
| `h_k` | `H(:,k)` | `N×1` channel user `k` → RIS |
| `G` | `G` | `N_R×N` channel RIS → BS |
| `H_k = diag(h_k)` | `diag(H(:,k))` | diagonal user channel |
| `A_k = G H_k` | `G*diag(H(:,k))` | cascaded user→RIS→BS channel |
| `gamma` | `gamma` | `N×1` RIS reflection vector |
| `X = gamma gamma^H` | `X` | lifted RIS variable (Approach 2) |
| `p_k` | `p(k)` | transmit power of user `k` |
| `c_k` | `C(:,k)` | linear receive filter for user `k` |
| `B` | `BW` / `BW0` | bandwidth (Hz / MHz) |
| `mu_k` | `mu` | inverse transmit-amplifier efficiency (≥ 1) |

## Key quantities

- **Colored-noise covariance** `W = sigma^2 I_{N_R} + sigma_RIS^2 G diag(gamma gamma^H) G^H`.
- **SINR** (Eq. 2): `SINR_k = p_k |c_k^H A_k gamma|^2 / ( c_k^H W c_k + sum_{m≠k} p_m |c_k^H A_m gamma|^2 )`.
- **Matrix R** (Eq. 6): `R = sum_k p_k H_k^H H_k + sigma_RIS^2 I_N`  → `func_R.m`.
- **Global reflection constraint**:
  - active (Eq. 7–8): `tr(R) ≤ tr(R gamma gamma^H) ≤ PR_max + tr(R)`;
  - nearly-passive (Eq. 14b): `tr(R gamma gamma^H) ≤ tr(R)`.
- **Total power / GEE denominator** (Eqs. 10–11): static circuit power `Pc`
  plus `sum_k mu_k p_k` plus the RIS amplification term. Implemented in
  `func_ptot*`, `func_gtot*`, `func_Xtot*`.
- **LMMSE receiver** (Eq. 44): `c_k = sqrt(p_k) M_k^{-1} A_k gamma`,
  `M_k = sum_{m≠k} p_m A_m gamma gamma^H A_m^H + W`.

## Nearly-passive = special case of active

Set `PR = 1` and `sigma_RIS = 0` (no amplification, no RIS noise). Then
Problem (12) reduces to Problem (14); the passive code mirrors the active code
with the `sigma_sqris` term removed.

## Default numerical parameters (paper Sec. V)

| Parameter | Value |
|-----------|-------|
| `K`, `N_R`, `N` | 4, 4, 100 |
| Bandwidth `B` | 20 MHz |
| `P0` (all nodes except RIS) | 40 dBm |
| `P0,RIS` (active / passive) | 30 / 20 dBm |
| `Pc,n` per element (active / passive) | 20 / 0 dBm |
| `PR_max` | 10 dBW |
| Noise PSD / figure | −174 dBm/Hz / 10 dB |
| Cell radius / BS–RIS distance | 100 m / 50 m |
| Heights: user / RIS / BS | [0,5] m / 15 m / 10 m |
| Path-loss exponent | 4 |
| Rician factors (RIS–BS / user–RIS) | 4 / 2 |
| Convergence tolerance `eps` | 1e-3 … 1e-6 |

> The original research scripts shipped with slightly different defaults
> (e.g. noise figure 2 dB, larger heights). The example scripts in `examples/`
> use the **paper** values above.
