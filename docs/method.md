# Method

This note derives what the library computes and explains the numerical choices. Code references point to `include/risgeo/`.

## 1. Network model

The typical user sits at the origin and is served by a base station (BS) at distance $r$.

* **Interfering BSs** form a Poisson point process (PPP) $\Phi_B$ of density $\lambda_B$ outside the disk of radius $r$ (the serving BS is the nearest one). All BSs transmit with power $P$ through antennas of gain $G$.
* **Fading.** BS links are Rayleigh faded: the power gain $h \sim \mathrm{Exp}(1)$, independent across links. The pathloss is $\ell(d) = (1+d)^{-4}$ (NLoS).
* **RISs** form a PPP of density $\lambda_R$ in the ring $[R_1, R_2]$ around the user, a Matérn-cluster-like model with a Poisson number of RISs per user. Each RIS has $N$ elements. By the central limit theorem its amplitude is Gaussian, $A \sim \mathcal{N}(N\mu_\zeta, N\sigma_\zeta^2)$, where $\mu_\zeta, \sigma_\zeta^2$ are the mean and variance of the product of two Rician gains. The reflected power is $P G^2 A^2\, \ell_{\mathrm{LoS}}(\text{BS} \to \text{RIS})\, \ell_{\mathrm{LoS}}(\text{RIS} \to \text{user})$ with $\ell_{\mathrm{LoS}}(d) = (1+d)^{-3}$.

With noise power $N_0$, the signal-to-interference-plus-noise ratio is

$$
\mathrm{SINR} = \frac{S + R}{I + N_0},\qquad S = PGh\,\ell(r),\quad I = \sum_{x\in\Phi_B} PGh_x\,\ell(|x|),\quad R = \sum_{y \in \Phi_R} R_y .
$$

The same machinery accepts other fading laws (Gamma for Nakagami-m, `laplace/Gamma.hpp`), binomial point processes (`laplace/BinomialPP.hpp`) and other pathloss models (`utils/pathloss.hpp`).

## 2. Laplace transforms

**Point processes** (`laplace/PoissonPP.hpp`). For a PPP of density $\lambda$ on a region $D$ with independent marks, the probability generating functional gives

$$
\mathcal{L}_{\sum_x f(x)}(s) = \exp\Bigl(-\lambda \int_D \bigl(1 - \mathcal{L}_{f(x)}(s)\bigr)\,dx\Bigr).
$$

For the interference, in polar coordinates, $\mathcal{L}_I(s) = \exp\bigl(-2\pi\lambda_B \int_r^\infty x\,(1 - \tfrac{1}{1 + sPG\ell(x)})\,dx\bigr)$, a 1D integral (`PoissonPP2D_Reduced`). For the RIS cluster it is a 2D integral over $(y, \varphi)$ in the ring (`PoissonPP2D`). The exponent per unit density, $\partial \log \mathcal{L} / \partial \lambda$, is exposed as `eval_exponent`. The research code uses it to differentiate performance with respect to the BS and RIS densities.

**Fading distributions** (`laplace/`). Exponential: $\mathcal{L}(s) = (1+vs)^{-1}$. Gamma: $(1+vs)^{-k}$. Scaled non-central chi-square with one degree of freedom, i.e. $(m + \sqrt{v}Z)^2$: $\mathcal{L}(s) = \exp\bigl(-\tfrac{s\,\lambda v}{1+2vs}\bigr)(1+2vs)^{-1/2}$ with $\lambda = m^2/v$. The `*_D1` functions return $d\mathcal{L}/ds$, which gives means and sensitivities.

## 3. Coverage with sign-indefinite effective interference

Normalise by the mean signal power and write coverage as $\mathbb{P}(h > X)$ with

$$
X = \frac{\theta (I + N_0) - R}{PG\,\ell(r)}, \qquad
\mathcal{L}_X(s) = \mathcal{L}_I\Bigl(\frac{s\theta}{PG\ell(r)}\Bigr)\, e^{-s\theta N_0/(PG\ell(r))}\, \mathcal{L}_R\Bigl(\frac{-s}{PG\ell(r)}\Bigr).
$$

**Without RISs**, $X \ge 0$ and $\mathbb{P}(h > X) = \mathbb{E}[e^{-X}] = \mathcal{L}_X(1)$: one 1D integral, about 2 µs (`CoverageProbability::calculate_rayleigh`).

**With RISs**, $X$ takes both signs, so $\mathbb{P}(h > X) = \mathbb{E}[e^{-X}\mathbf{1}\{X > 0\}] + \mathbb{P}(X \le 0)$. Gil-Pelaez inversion gives the second term from the characteristic function $\varphi_X(u) = \mathcal{L}_X(-iu)$:

$$
\mathbb{P}(X > 0) = \frac12 + \frac1\pi \int_0^\infty \frac{\operatorname{Im} \mathcal{L}_X(-iu)}{u}\,du .
$$

For the first term, apply the same formula under the tilted measure $d\tilde{\mathbb{P}} = e^{-X} d\mathbb{P} / \mathcal{L}_X(1)$, whose characteristic function is $\mathcal{L}_X(1 - iu) / \mathcal{L}_X(1)$:

$$
\mathbb{E}[e^{-X}\mathbf{1}\{X > 0\}] = \frac{\mathcal{L}_X(1)}{2} + \frac1\pi \int_0^\infty \frac{\operatorname{Im} \mathcal{L}_X(1 - iu)}{u}\,du .
$$

Adding the two,

$$
\boxed{\;\mathbb{P}(h > X) = \tfrac12\bigl(1 + \mathcal{L}_X(1)\bigr) + \frac1\pi \int_0^\infty \frac{\operatorname{Im}\bigl[\mathcal{L}_X(1 - iu) - \mathcal{L}_X(-iu)\bigr]}{u}\,du\;}
$$

(`inversion/GilPelaez.hpp`, `exponential_exceeds`). Both oscillating terms sit in one integrand, so they partially cancel before integration. The integrand is finite at $u \to 0$ because $\operatorname{Im}\mathcal{L}_X(a) = 0$ for real $a$. The half-line integral is evaluated with exp-sinh double-exponential quadrature, which handles the slow algebraic/oscillatory decay.

## 4. Ergodic rate

$\mathbb{E}[\ln(1 + \mathrm{SINR})] = \int_0^\infty \mathbb{P}(\mathrm{SINR} > \theta)\,\frac{d\theta}{1+\theta}$. Without RISs this is a 1D integral of a closed form (Gauss-Kronrod). With RISs, the inversion integral is nested inside the $\theta$ integral. Substituting $u = e^{w}$ turns $du/u$ into $dw$ and removes the $1/u$ singularity, leaving a smooth, bounded integrand on $[0, \theta_{\max}] \times [-20, 20]$ for the adaptive cubature (`GilPelaez::integrated_exponential_exceeds`). The truncations ($\theta_{\max} = 10^4$, $|w| \le 20$) are in `geometry/Parameter.hpp`.

## 5. Numerical integration

**Adaptive cubature** (`integration/cuhre/`). This is a from-scratch implementation of the globally adaptive strategy of Berntsen, Espelid and Genz (the algorithm behind Cuba's Cuhre):

* each subregion is integrated with the fully symmetric degree-7 rule of Genz and Malik ($2^N + 2N^2 + 4N + 1$ points); three embedded null rules of degrees 5, 5 and 3 give a robust error estimate;
* the subregion with the largest error is bisected along the axis with the largest fourth divided difference;
* subregions live in a max-heap keyed on their error, and the totals are maintained incrementally (exact sums are refreshed periodically and before accepting convergence), so a bisection costs $O(\log n)$;
* with OpenMP, $3\times$ the number of threads subregions are bisected per sweep;
* integrands may return real or complex values (templates, no virtual calls in the inner loop).

`integration/genz_malik/GM2D.hpp` is a 2D variant with precomputed tables for `float`, `double` and `long double`, selected at compile time as Boost.Math does. One caveat from experience: like every Genz-Malik implementation, it can be unreliable when the integrand varies at very different rates along the two axes, e.g. $e^{-xy^4}$.

**1D quadrature** uses Boost.Math: 15-point Gauss-Kronrod for finite and semi-infinite smooth integrands, exp-sinh for $[0, \infty)$ with endpoint behaviour.

## 6. Validation

`validation/monte_carlo.cpp` simulates the model of Section 1 directly (also mirrored by the original MATLAB scripts in `validation/matlab/`) and compares it with the analytical results at 13 SINR thresholds, with and without RISs:

* with 4 million samples, the largest deviation is 2.1 standard errors over 26 points, and the ergodic rates agree within 1.1 standard errors;
* the simulation truncates interferers at 3 km, so the analytical side integrates over the same region. Integrating to infinity lowers coverage by about 0.2 % (relative), which the comparison detects;
* `tests/regression_test.cpp` checks the integrators on integrals with known values, the derivative functions against finite differences, the point-process identities $\mathcal{L} = e^{\lambda \cdot \texttt{eval\_exponent}}$, and the published example values.

## 7. The same engine in finance

With $X = \ln(S_T/K)$ under the risk-neutral measure, $\mathcal{L}_X(s) = K^{s}\,\varphi_{\ln S_T}(is)$ and

$$
\text{digital call} = e^{-rT}\,\mathbb{P}(X > 0), \qquad
\text{call} = e^{-rT} K \bigl(\mathbb{E}[e^{X}\mathbf{1}\{X>0\}] - \mathbb{P}(X > 0)\bigr).
$$

The second expectation is the tilted inversion with $a = -1$, the change to the stock numeraire, mirroring the $a = +1$ tilt that Rayleigh fading induces above. `examples/finance/pricing.hpp` implements Black-Scholes and Heston (in the "little Heston trap" form of Albrecher et al.). Black-Scholes prices match the closed form to $10^{-13}$; the Heston benchmark of Fang and Oosterlee (5.785155434) is reproduced to $4\times 10^{-10}$; and a Monte Carlo cross-check agrees at every strike.

## References

* J. Gil-Pelaez, "Note on the inversion theorem," *Biometrika*, 1951.
* A. C. Genz, A. A. Malik, "An adaptive algorithm for numerical integration over an N-dimensional rectangular region," *J. Comput. Appl. Math.*, 1980.
* J. Berntsen, T. O. Espelid, A. Genz, "An adaptive algorithm for the approximate calculation of multiple integrals," *ACM Trans. Math. Softw.*, 1991.
* T. Hahn, "Cuba: a library for multidimensional numerical integration," *Comput. Phys. Commun.*, 2005.
* S. L. Heston, "A closed-form solution for options with stochastic volatility with applications to bond and currency options," *Rev. Financ. Stud.*, 1993.
* H. Albrecher, P. Mayer, W. Schoutens, J. Tistaert, "The little Heston trap," *Wilmott Magazine*, 2007.
* F. Fang, C. W. Oosterlee, "A novel pricing method for European options based on Fourier-cosine series expansions," *SIAM J. Sci. Comput.*, 2008.
* G. Sun, F. Baccelli, K. Feng, L. Uzeda Garcia, S. Paris, "A stochastic geometry framework for performance analysis of RIS-assisted OFDM cellular networks," *IEEE Trans. Wireless Commun.*, 2025.
