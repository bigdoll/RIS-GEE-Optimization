# Equation & algorithm mapping

## Example scripts → paper figures

| Script | Paper item |
|--------|-----------|
| `examples/run_active_ris.m`, `examples/run_passive_ris.m` | Figs. 2-3 (GEE / sum-rate vs. transmit power) |
| `examples/figures/fig4_gee_vs_pcn.m` | Fig. 4 (GEE vs. active per-element static power) |
| `examples/figures/fig5_gee_vs_N.m` | Fig. 5 (GEE vs. number of RIS elements) |
| `examples/figures/fig6_global_vs_local.m` | Fig. 6 (global vs. local reflection constraint) |
| `examples/figures/fig7_gee_vs_position.m` | Fig. 7 (GEE vs. RIS-BS distance) |
| `examples/figures/tables_convergence.m` | Tables I-II (iterations & convergence time) |
| `examples/figures/gee_active.m`, `gee_passive.m` | per-realization evaluation helpers |

## Function index

Every MATLAB function in this repository, the paper equation/algorithm it implements, and a one-line description. Paths are relative to the repository root.

| File | Paper ref. | Purpose |
|------|-----------|---------|
| `src/active_ris/LMMSE_receiver_active.m` | Eq. (44) | Closed-form linear MMSE receive filters c_k (active RIS). |
| `src/active_ris/SINR_active.m` | Eq. (2) | Per-user SINR with an active RIS (includes RIS noise amplification). |
| `src/active_ris/approach1_alternating/altopt1_active.m` | Algorithm 3, Problem (12) | ALGORITHM 3 (active): alternating optimization of power p, RIS vector gamma and MMSE filters C. |
| `src/active_ris/approach1_alternating/data_rate_active.m` | Eq. (23) | Per-user and sum rate for the active branch (bandwidth in Mbit/s). |
| `src/active_ris/approach1_alternating/parameters_active.m` | Eqs. (34)-(36) | Sequential-approximation constants A_bar,B_bar,D_bar,E_bar,F_bar for the RIS subproblem. |
| `src/active_ris/approach1_alternating/power/SRp_active.m` | Eq. (40) | MMSE sum-rate written as g1(p) - g2(p). |
| `src/active_ris/approach1_alternating/power/SRpapprox_active.m` | Eq. (41) | Concave surrogate of the sum-rate in p:  g1 - g2 - grad g2'(p-p_bar). |
| `src/active_ris/approach1_alternating/power/SRpapprox_active_cvx.m` | Eq. (41) | CVX-expression form of the power surrogate. |
| `src/active_ris/approach1_alternating/power/func_gtot_active.m` | Eq. (11) | Total consumed power as a function of gamma (RIS subproblem). |
| `src/active_ris/approach1_alternating/power/func_gtot_active_cvx.m` | Eq. (11) | CVX-expression form of func_gtot_active (quad_form, convex-certified). |
| `src/active_ris/approach1_alternating/power/func_ptot_active.m` | Eqs. (10)-(11),(39) | Total consumed power as a function of p (power subproblem). |
| `src/active_ris/approach1_alternating/power/funcg1_active.m` | Eq. (40) | Concave term g1(p) of the sum-rate decomposition. |
| `src/active_ris/approach1_alternating/power/funcg1_active_cvx.m` | Eq. (40) | CVX-expression form of g1(p). |
| `src/active_ris/approach1_alternating/power/funcg2_active.m` | Eq. (40) | Convex term g2(p) of the sum-rate decomposition. |
| `src/active_ris/approach1_alternating/power/gradg2_active.m` | Eq. (42) | Gradient of f2(p) (numerator of g2) used for the linearization. |
| `src/active_ris/approach1_alternating/power/p_cvxopt_active.m` | Algorithm 2, Problem (43) | ALGORITHM 2 (active): sequential fractional-programming power optimization. |
| `src/active_ris/approach1_alternating/ris_reflection/SRgapprox_active.m` | Eqs. (35)-(36) | Concave surrogate of the MMSE sum-rate in gamma (true value). |
| `src/active_ris/approach1_alternating/ris_reflection/SRgapprox_active_cvx.m` | Eqs. (35)-(36) | CVX-expression form of the gamma surrogate (objective inside cvx_begin). |
| `src/active_ris/approach1_alternating/ris_reflection/gamma_cvxopt_active.m` | Algorithm 1, Problem (38) | ALGORITHM 1 (active): sequential fractional-programming RIS optimization. |
| `src/active_ris/approach2_mmse/altopt2_active.m` | Algorithm 6, Problem (47) | ALGORITHM 6 (active): alternating optimization of X = gamma*gamma^H and p (MMSE embedded). |
| `src/active_ris/approach2_mmse/power/SRpapprox_active2.m` | Eq. (55) | Concave surrogate of the sum-rate in p (MMSE-embedded branch). |
| `src/active_ris/approach2_mmse/power/SRpapprox_active_cvx2.m` | Eq. (55) | CVX-expression form of the power surrogate (MMSE branch). |
| `src/active_ris/approach2_mmse/power/func_ptot_active2.m` | Eqs. (10)-(11) | Total consumed power as a function of p (MMSE branch). |
| `src/active_ris/approach2_mmse/power/funcg1_active2.m` | Eq. (55) | Concave term g1(p) of the MMSE-embedded sum-rate. |
| `src/active_ris/approach2_mmse/power/funcg1_active_cvx2.m` | Eq. (55) | CVX-expression form of g1(p) (MMSE branch). |
| `src/active_ris/approach2_mmse/power/funcg2_active2.m` | Eq. (55) | Convex term g2(p) of the MMSE-embedded sum-rate. |
| `src/active_ris/approach2_mmse/power/gradg2_active2.m` | Eq. (55) | Gradient of f2(p) (MMSE branch) for the linearization. |
| `src/active_ris/approach2_mmse/power/p_cvxopt_active2.m` | Algorithm 5, Problem (56) | ALGORITHM 5 (active): sequential power optimization (MMSE branch). |
| `src/active_ris/approach2_mmse/ris_reflection/G1func_active2.m` | Eq. (50) | Concave term F1(X). |
| `src/active_ris/approach2_mmse/ris_reflection/G1func_active_cvx2.m` | Eq. (50) | CVX-expression form of F1(X). |
| `src/active_ris/approach2_mmse/ris_reflection/G2func_active2.m` | Eq. (50) | Concave term F2(X). |
| `src/active_ris/approach2_mmse/ris_reflection/G2grad_active2.m` | Eq. (53) | Gradient grad F2(X) for the linearization. |
| `src/active_ris/approach2_mmse/ris_reflection/SRX_active2.m` | Eq. (50) | MMSE sum-rate in X as F1(X) - F2(X). |
| `src/active_ris/approach2_mmse/ris_reflection/SRXapprox_active2.m` | Eq. (52) | Concave surrogate of the sum-rate in X (F2 linearized). |
| `src/active_ris/approach2_mmse/ris_reflection/X_cvxopt_active2.m` | Algorithm 4, Problem (54) | ALGORITHM 4 (active): semidefinite-relaxation RIS optimization in X. |
| `src/active_ris/approach2_mmse/ris_reflection/func_Xtot_active2.m` | Eq. (51) | Total consumed power as a function of X:  tr(R X) + ... . |
| `src/common/available_power.m` | System model, Sec. V | Build the transmit-power sweep: convert a dB gain range to linear Watts. |
| `src/common/data_rate.m` | Eq. (23)/(123) | Per-user and sum rate  R_k = B*log2(1+SINR_k). |
| `src/common/func_R.m` | Eq. (6) | Positive-definite matrix  R = sum_k p_k H_k^H H_k + sigma_RIS^2 I_N. |
| `src/common/generate_channels.m` | System model, Sec. II / V | Geometry-based Rician channels for the RIS-aided uplink (G: RIS->BS, H: users->RIS). |
| `src/common/noise_power.m` | System model, Sec. II | Thermal-noise power (variance) from PSD, noise figure and bandwidth. |
| `src/common/static_power_consumption.m` | Eq. (10) | Static circuit power  Pc = P0 + N*Pc,n + P0,RIS. |
| `src/passive_ris/LMMSE_receiver_passive.m` | Eq. (44) | Closed-form linear MMSE receive filters c_k (nearly-passive RIS). |
| `src/passive_ris/SINR_passive.m` | Eq. (2), sigma_RIS=0 | Per-user SINR with a nearly-passive RIS (no RIS noise amplification). |
| `src/passive_ris/approach1_alternating/alt_opt1.m` | Algorithm 3, Problem (14) | ALGORITHM 3 (nearly-passive): alternating optimization of p, gamma and MMSE filters C. |
| `src/passive_ris/approach1_alternating/parameters_passive.m` | Eqs. (34)-(36) | Sequential-approximation constants A_bar..F_bar for the RIS subproblem. |
| `src/passive_ris/approach1_alternating/power/SRpapprox_passive.m` | Eq. (41) | Concave surrogate of the sum-rate in p:  g1 - g2 - grad g2'(p-p_bar). |
| `src/passive_ris/approach1_alternating/power/SRpapprox_passive_cvx.m` | Eq. (41) | CVX-expression form of the power surrogate. |
| `src/passive_ris/approach1_alternating/power/func_g1.m` | Eq. (40) | Concave term g1(p) of the sum-rate decomposition. |
| `src/passive_ris/approach1_alternating/power/func_g1_cvx.m` | Eq. (40) | CVX-expression form of g1(p). |
| `src/passive_ris/approach1_alternating/power/func_g2.m` | Eq. (40) | Convex term g2(p) of the sum-rate decomposition. |
| `src/passive_ris/approach1_alternating/power/gradient_g2.m` | Eq. (42) | Gradient of f2(p) (numerator of g2) for the linearization. |
| `src/passive_ris/approach1_alternating/power/gradient_g2_cvx.m` | Eq. (42) | CVX-expression form of grad f2(p). |
| `src/passive_ris/approach1_alternating/power/p_cvxopt.m` | Algorithm 2, Problem (43) | ALGORITHM 2 (nearly-passive): sequential power optimization. |
| `src/passive_ris/approach1_alternating/ris_reflection/SRgapprox_passive.m` | Eqs. (35)-(36) | Concave surrogate of the MMSE sum-rate in gamma (true value). |
| `src/passive_ris/approach1_alternating/ris_reflection/SRgapprox_passive_cvx.m` | Eqs. (35)-(36) | CVX-expression form of the gamma surrogate. |
| `src/passive_ris/approach1_alternating/ris_reflection/gamma_cvxopt.m` | Algorithm 1, Problem (38) | ALGORITHM 1 (nearly-passive): sequential fractional-programming RIS optimization. |
| `src/passive_ris/approach2_mmse/alt_opt2.m` | Algorithm 6, Problem (47) | ALGORITHM 6 (nearly-passive): alternating optimization of X and p (MMSE embedded). |
| `src/passive_ris/approach2_mmse/power/SRpapprox2_passive.m` | Eq. (55) | Concave surrogate of the sum-rate in p (MMSE-embedded branch). |
| `src/passive_ris/approach2_mmse/power/SRpapprox2_passive_cvx.m` | Eq. (55) | CVX-expression form of the power surrogate (MMSE branch). |
| `src/passive_ris/approach2_mmse/power/func2_g1.m` | Eq. (55) | Concave term g1(p) of the MMSE-embedded sum-rate. |
| `src/passive_ris/approach2_mmse/power/func2_g2.m` | Eq. (55) | Convex term g2(p) of the MMSE-embedded sum-rate. |
| `src/passive_ris/approach2_mmse/power/func2_gradient.m` | Eq. (55) | Gradient of f2(p) (MMSE branch) for the linearization. |
| `src/passive_ris/approach2_mmse/power/power2_cvxopt.m` | Algorithm 5, Problem (56) | ALGORITHM 5 (nearly-passive): sequential power optimization (MMSE branch). |
| `src/passive_ris/approach2_mmse/ris_reflection/cvxopt_X.m` | Algorithm 4, Problem (54) | ALGORITHM 4 (nearly-passive): semidefinite-relaxation RIS optimization in X. |
| `src/passive_ris/approach2_mmse/ris_reflection/funcG1_cvx.m` | Eq. (50) | CVX-expression form of F1(X). |
| `src/passive_ris/approach2_mmse/ris_reflection/funcG2_grad.m` | Eq. (53) | Gradient grad F2(X) for the linearization. |
| `src/passive_ris/approach2_mmse/ris_reflection/funcG_inner.m` | Eq. (50) | Inner per-user log-det term of the MMSE sum-rate. |
| `src/passive_ris/approach2_mmse/ris_reflection/func_G.m` | Eq. (50) | Concave term F1(X) of the MMSE sum-rate. |
| `src/passive_ris/approach2_mmse/ris_reflection/func_SRmmse.m` | Eq. (50) | MMSE sum-rate in X as F1(X) - F2(X). |
| `src/passive_ris/approach2_mmse/ris_reflection/func_SRmmse_approx.m` | Eq. (52) | Concave surrogate of the sum-rate in X (F2 linearized). |
| `src/passive_ris/func_ptot.m` | Eqs. (10),(14) | Total consumed power  sum_k mu*p_k + Pc (nearly-passive RIS). |
