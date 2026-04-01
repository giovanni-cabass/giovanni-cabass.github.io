---
layout: post
title: "An exact formula for nonlinear inference and its perturbative consequences"
date: 2026-03-31
---

This post works through Appendix A of [arXiv:2307.04706](https://arxiv.org/abs/2307.04706), joint work with Marko Simonović and Matias Zaldarriaga. The question is: given a nonlinear forward model, what does the likelihood look like once the latent field has been marginalized out in the small-noise limit? The short answer is that the integral can be done exactly by a change of variables, and that two perturbative consequences — unbiasedness of the score and equality of the Fisher information with the linear estimator — both follow from a single integration by parts. The cosmological application is what prompted it, but the argument is purely probabilistic.

---

## The setup

Let $x \in \mathbb{R}^N$ be a latent variable drawn from a Gaussian prior,

$$
\mathcal{P}_0(x) \propto \exp\!\left(-\tfrac{1}{2}\, (C^{-1})_{\mu\nu} x^\mu x^\nu\right)\,\,,
$$

where $C$ is the $N \times N$ prior covariance. Think of $N$ as very large — in the cosmological application, $N \sim 10^6$ is the number of Fourier modes resolved by a galaxy survey.

We have a nonlinear forward model $F(x, \bm{\theta})$ relating $x$ to the observed variables $y$, where $\bm{\theta}$ are parameters of interest. The true latent $\hat{x}$ and true parameters $\hat{\bm{\theta}}$ generate the observed data via

$$
\hat{y} = F(\hat{x},\, \hat{\bm{\theta}})\,\,.
$$

We do not observe $x$ directly so we need to marginalize over them. The quantity of interest is then the posterior

$$
{\cal P}(\bm{\theta}\mid\hat{y}) \propto \int{\cal D}x\,{\cal L}(\hat{y}\mid x,\bm{\theta})\,\mathcal{P}_0(x)\,\,,
$$

where $\mathcal{D}x = \prod_{\mu=1}^N {\rm d}x^\mu$ is the flat measure on $\mathbb{R}^N$, ${\cal L}$ is the likelihood and $\hat{y}$ is the observed data.

The key assumption now is that the likelihood has a well-defined zero-noise (signal-dominated) limit in which it becomes a Dirac delta:

$$
{\cal L}(y\mid x,\bm{\theta}) \;\xrightarrow{\text{zero noise}}\; \delta^{(N)}_{\rm D}\!\bigl(y - F(x, \bm{\theta})\bigr)\,\,.
$$

This is the idealization in which the data perfectly constrain $x$ given $\bm{\theta}$. Marginalizing over $x$ gives the posterior for $\bm{\theta}$:

$$
\mathcal{P}(\hat{y} \mid \bm{\theta}) \propto \int \mathcal{D}x\; \delta^{(N)}_{\rm D}\!\bigl(\hat{y} - F(x, \bm{\theta})\bigr)\, \mathcal{P}_0(x)\,\,.
$$

We assume $F$ is smooth and that the equation $F(x, \bm{\theta}) = \hat{y}$ has a discrete set of solutions. Among these we take the branch connected to the linear approximation. This is inspired by cosmological perturbation theory where, for small $x$, we can expand

$$
F^\mu(x,\bm{\theta}) = M^\mu_{\hp{\mu}\nu}(\bm{\theta})\,x^\nu + {\cal O}(x^2)\,\,.
$$

On this branch the Jacobian matrix $\mathcal{D}^\mu_{\hp{\mu}\nu} = \partial F^\mu / \partial x^\nu$ is invertible; local invertibility is all the argument requires.

The question: **what are the first and second derivatives of ${-\ln \mathcal{P}(\hat{y} \mid \bm{\theta})}$ with respect to $\bm{\theta}$, and what do they look like on average?**

**Note on $N \to \infty$.** The argument below carries through formally when $N \to \infty$, with the ordinary integral replaced by the path integral of statistical field theory. This is exactly the setting of the original paper. Care is required with regularization and consequently renormalization, which we do not discuss here.

---

## Index notation

We use a compact index notation in which the algebra looks like ordinary matrix calculus and, as will become clear, carries a natural differential-geometric interpretation.

- **Greek indices** $(\mu, \nu, \rho, \ldots)$ run from $1$ to $N$ and label components of $x$ (or $y$). Derivatives with respect to $x^\mu$ are written $\partial_\mu$.
- **Latin indices** $(i, j, \ldots)$ label the parameters $\theta_i$. The derivative $\partial_i$ is always evaluated at the fiducial value $\hat{\bm{\theta}}$.

Repeated Greek indices are summed (Einstein convention). Quantities decorated with a hat are evaluated at the true point $(\hat{x}, \hat{\bm{\theta}})$ — the same $\hat{x}$ and $\hat{\bm{\theta}}$ that generated $\hat{y} = F(\hat{x}, \hat{\bm{\theta}})$. Quantities with a tilde are evaluated at $(\hat{x}, \bm{\theta})$ — same true $\hat{x}$, but with $\bm{\theta}$ left general.

The forward model and its higher derivatives at the truth:

$$
\hat F^\mu = F^\mu(\hat{x}, \hat{\bm{\theta}})\,\,, \qquad \hat{\mathcal{D}}^\mu_{\hp{\mu}\nu} = \partial_\nu \hat F^\mu\,\,, \qquad \hat{\mathcal{D}}^\mu_{\hp{\mu}\nu_1\cdots\nu_n} = \partial_{\nu_1}\cdots\partial_{\nu_n} \hat F^\mu\,\,.
$$

The matrix $\hat{\mathcal{I}}^\mu_{\hp{\mu}\nu}$ is the inverse of the Jacobian:

$$
\hat{\mathcal{I}}^\mu_{\hp{\mu}\rho}\, \hat{\mathcal{D}}^\rho_{\hp{\rho}\nu} = \delta^\mu_{\hp{\mu}\nu}\,\,.
$$

The *mismatch* between the forward model at $\bm{\theta}$ and at the true point:

$$
X^\mu = \hat F^\mu - \tilde F^\mu\,\,, \qquad \tilde F^\mu = F^\mu(\hat{x}, \bm{\theta})\,\,.
$$

The tilded Jacobian and its inverse are defined analogously at $(\hat{x}, \bm{\theta})$:

$$
\tilde{\mathcal{D}}^\mu_{\hp{\mu}\nu} = \partial_\nu \tilde F^\mu\,\,, \qquad \tilde{\mathcal{D}}^\mu_{\hp{\mu}\nu_1\cdots\nu_n} = \partial_{\nu_1}\cdots\partial_{\nu_n} \tilde F^\mu\,\,, \qquad \tilde{\mathcal{I}}^\mu_{\hp{\mu}\rho}\, \tilde{\mathcal{D}}^\rho_{\hp{\rho}\nu} = \delta^\mu_{\hp{\mu}\nu}\,\,.
$$

Note $X^\mu = 0$ when $\bm{\theta} = \hat{\bm{\theta}}$; this is the key small quantity for a perturbative expansion around the true point.

---

## Evaluating the integral

We shift $x^\mu = \hat{x}^\mu + \Delta^\mu$ around the true latent $\hat{x}$ and expand the argument of the delta function in $\Delta$:

$$
F^\mu(\hat{x}+\Delta,\, \bm{\theta}) - \hat{F}^\mu = \underbrace{\tilde{\mathcal{D}}^\mu_{\hp{\mu}\nu}\Delta^\nu + \frac{1}{2}\tilde{\mathcal{D}}^\mu_{\hp{\mu}\nu\rho}\Delta^\nu\Delta^\rho + \cdots}_{\equiv\; G(\Delta)^\mu} - X^\mu\,\,.
$$

The delta function enforces $G(\Delta)^\mu = X^\mu$. Perform the change of variables $\Delta^\mu \to Y^\mu = G(\Delta)^\mu$, a smooth invertible map in $\mathbb{R}^N$ (on the chosen branch). The Jacobian of $G$ contributes $\lvert\det \partial G/\partial \Delta\rvert$; the delta function sets $Y^\mu = X^\mu$. The result is exact:

$$
-\ln \mathcal{P}(\hat{y} \mid \bm{\theta}) = -\ln \mathcal{P}\!\bigl(\hat{x} + G^{-1}(X)\bigr) - \ln \left|\det \frac{\partial G^{-1}(X)}{\partial X}\right|\,\,.
$$

Two terms: the **prior** evaluated at the inferred latent variable $\hat{x} + G^{-1}(X)$, and a **log-Jacobian** from the change of variables. No approximation has been made; this is exact for any smooth $F$.

A practical consequence: all derivatives of $-\ln \mathcal{P}(\hat{y} \mid \bm{\theta})$ with respect to $\bm{\theta}$ follow from derivatives of $F$ with respect to $x$ and $\bm{\theta}$ — no additional integrals are needed. In an ML setting, this means the score and Fisher information are accessible by automatic differentiation through $F$.

The explicit series for $G^{-1}(X)$ follows from the inverse function theorem:

$$
G^{-1}(X)^\mu = \tilde{\mathcal{I}}^\mu_{\hp{\mu}\nu} X^\nu - \frac{1}{2} \tilde{\mathcal{I}}^\nu_{\hp{\nu}\beta}\, \tilde{\mathcal{I}}^\mu_{\hp{\mu}\rho}\, \tilde{\mathcal{I}}^\sigma_{\hp{\sigma}\alpha}\, \tilde{\mathcal{D}}^\rho_{\hp{\rho}\sigma\nu}\, X^\alpha X^\beta + \cdots
$$

(using $\partial_\nu \tilde{\mathcal{I}}^\mu_{\hp{\mu}\alpha} = -\tilde{\mathcal{I}}^\mu_{\hp{\mu}\rho}\, \tilde{\mathcal{D}}^\rho_{\hp{\rho}\sigma\nu}\, \tilde{\mathcal{I}}^\sigma_{\hp{\sigma}\alpha}$). All coefficients are evaluated at $(\hat{x}, \bm{\theta})$ (hence all the matrices carry still a dependence on $\bm{\theta}$). This is a series in the mismatch $X^\mu$, which is small when $\bm{\theta}$ is close to $\hat{\bm{\theta}}$, but is exact in $\hat{x}$ — the expansion is perturbative in the parameters, not in the latent field.

**Beyond the zero-noise limit.** The formula above assumes the likelihood collapses exactly to a Dirac delta. One can go beyond this for any likelihood of the form $\mathcal{L}(y \mid x,\bm{\theta}) = p(y - F(x,\bm{\theta}))$, where $p$ is the noise distribution. By Fourier inversion that

$$p(\epsilon) = e^{W({\rm i}\partial_\epsilon)}\,\delta^{(N)}_{\rm D}(\epsilon)\,\,,$$

where (being slightly careless about overall normalization factors)

$$W(J) = \ln Z(J)\,\,, \qquad Z(J) = \int {\cal D}\epsilon\; e^{\mathrm{i} J_\mu \epsilon^\mu}\, p(\epsilon)$$

is the cumulant generating function of the noise, and the substitution $J_\mu \to {\rm i}\partial_{\hat{y}^\mu}$ is understood. For a Gaussian with covariance $N^{\mu\nu}$: $W(J) = -\frac{1}{2}N^{\mu\nu}J_\mu J_\nu$, giving $W({\rm i}\partial_{\hat{y}^\mu}) = \frac{1}{2}N^{\mu\nu}\partial_\mu\partial_\nu$, which reproduces the familiar result. For non-Gaussian noise, higher cumulants contribute additional derivative operators. Since $\partial_{\hat{y}^\mu}$ acts on $\hat{y}$ and not on $x$, it commutes with the integral over $x$, and the full posterior is

$$\mathcal{P}_\epsilon(\hat{y} \mid \bm{\theta}) = e^{W(\mathrm{i}\partial_{\hat{y}^\mu})}\,\mathcal{P}(\hat{y} \mid \bm{\theta})\,\,,$$

where $\mathcal{P}$ is the zero-noise posterior given by the exact formula above. The expansion is controlled whenever $W$ is expandable in powers of $J$ with a small parameter (the noise-to-signal ratio). In particular, one can work through the unbiasedness and the Fisher matrix calculations of the following sections at each order in this expansion.

---

## Unbiasedness via integration by parts

We want to show $\langle {-\partial_i \ln \mathcal{P}(\hat{y} \mid \hat{\bm{\theta}})} \rangle = 0$, where the expectation is over $\hat{x} \sim \hat{\mathcal{P}}_0(\hat{x})$ (we use the hat to denote the prior evaluated at the fiducial $x$, i.e. $\hat{x}$).

Differentiating the series for $G^{-1}(X)$ with respect to $\theta_i$ and evaluating at $\hat{\bm{\theta}}$ (where $X^\mu = 0$), one can rewrite the result somewhat suggestively as

$$
-\partial_i \ln \mathcal{P}(\hat{y} \mid \hat{\bm{\theta}}) = \underbrace{(\partial_\nu \ln \hat{\mathcal{P}})}_{\hp{\Gamma^\rho_{\rho\nu}\,}=\,\Gamma^\rho_{\rho\nu}} \hat{\mathcal{I}}^\nu_{\hp{\nu}\mu}\, \partial_i \hat F^\mu + \partial_\nu\big\{\underbrace{\hat{\mathcal{I}}^\nu_{\hp{\nu}\mu}\, \partial_i \hat F^\mu}_{\hp{V^\nu_i\,}=\,V^\nu_i}\big\} = \nabla_\nu V^\nu_i\,\,.
$$

Here $\partial_\nu \ln \hat{\mathcal{P}}$ plays the role of the Christoffel symbol $\Gamma^\rho_{\rho\nu}$ for the metric with volume element ${\cal D}\hat{x}\,\hat{\mathcal{P}}(\hat{x})$, and $V^\nu_i = \hat{\mathcal{I}}^\nu_{\hp{\nu}\mu}\,\partial_i\hat{F}^\mu$ is a vector field. The expression is therefore the covariant divergence $\nabla_\nu V^\nu_i$. Averaging over $\hat{x} \sim \hat{\mathcal{P}}_0(x)$:

$$
\left\langle {-\partial_i \ln \mathcal{P}(\hat{y} \mid \hat{\bm{\theta}})} \right\rangle = \int \mathcal{D}\hat{x}\; \hat{\mathcal{P}}_0(\hat{x})\,\nabla_\nu V^\nu_i = \int \mathcal{D}\hat{x}\; \partial_\nu\!\left[\hat{\mathcal{P}}_0(\hat{x})\, V^\nu_i\right] = 0\,\,.
$$

This is Stokes' theorem in the space weighted by the prior, or equivalently integration by parts in $\mathbb{R}^N$: the boundary term vanishes because the Gaussian prior $\hat{\mathcal{P}}_0(\hat{x})$ decays exponentially as $\lvert\hat{x}\rvert \to \infty$. More generally, the same argument goes through for any prior that decays sufficiently fast at infinity — Gaussianity is not essential, only the vanishing of the boundary term.

An interesting avenue is to push the geometric structure further. The weighted integrals above live naturally on a space whose metric, for these purposes, can be taken to be conformally flat with conformal factor given by a suitable power of the prior:

$$g_{\mu\nu}(x) = \mathcal{P}_0^{2/N}(x)\,\delta_{\mu\nu}\,\,.$$

Whether the full Riemannian geometry of this space — geodesics, curvature, natural frames — offers a richer language for the structure of the posterior is an open question worth exploring.

---

## Second derivatives and Fisher information

The second derivative $-\partial_i\partial_j \ln \mathcal{P}$ has four terms: two from differentiating the prior piece, two from differentiating the log-Jacobian. Writing $Z^\mu_{\hp{\mu}\nu} = \partial G^{-1}(X)^\mu / \partial X^\nu$, one finds

$$
{-\partial_i\partial_j \ln \mathcal{P}} = {-(\partial_\mu\partial_\nu \ln \hat{\mathcal{P}}_0)}\,\partial_i G^{-1}(X)^\mu\, \partial_j G^{-1}(X)^\nu - (\partial_\mu \ln \hat{\mathcal{P}}_0)\,\partial_i\partial_j G^{-1}(X)^\mu + J^{11}_{ij} + J^{02}_{ij}\,\,.
$$

The key structural fact is that the "second-derivative" terms — involving $\partial_\mu \ln \hat{\cal P}_0$ and

$$J^{02}={-(Z^{-1})^\nu_{\hp{\nu}\mu}\, \partial_i\partial_j Z^\mu_{\hp{\mu}\nu}}\,\,,$$

all now evaluated at $\hat{\bm\theta}$ — combine into a total divergence:

$$
{-(\partial_\mu \ln \hat{\mathcal{P}}_0)}\,\partial_i\partial_j G^{-1}(X)^\mu - (Z^{-1})^\nu_{\hp{\nu}\mu}\, \partial_i\partial_j Z^\mu_{\hp{\mu}\nu} = \nabla_\mu\!\left\{\partial_i\partial_j G^{-1}(X)^\mu\right\}\,\,.
$$

By the same integration by parts against the Gaussian prior, this vanishes on average. What remains is the **Fisher information**:

$$
F_{ij} = \langle{-\partial_i\partial_j \ln \mathcal{P}}\rangle = \Big\langle({-\partial_\mu\partial_\nu\ln\hat{\cal P}_0}) \; \hat{\cal I}^\mu_{\hp{\mu}\rho}\hat{\cal I}^\nu_{\hp{\nu}\sigma}\partial_i\hat{F}^\rho\partial_j\hat{F}^\sigma + J^{11}_{ij}\Big\rangle\,\,,
$$

where we recognize the usual term involving first derivatives of the forward model squared, i.e.

$$
\partial_i G^{-1}[X]^\mu\,\partial_j G^{-1}[X]^\nu = \hat{\cal I}^\mu_{\hp{\mu}\rho}\hat{\cal I}^\nu_{\hp{\nu}\sigma}\partial_i\hat{F}^\rho\partial_j\hat{F}^\sigma
$$

at $\bm{\theta}=\hat{\bm{\theta}}$, and the Jacobian term

$$J_{ij}^{11}=(Z^{-1})^\nu_{\hp{\nu}\rho}(Z^{-1})^\sigma_{\hp{\sigma}\mu}\partial_iZ^\mu_{\hp{\mu}\nu}\,\partial_jZ^\rho_{\hp{\rho}\sigma}$$

that takes a complicated form in terms of first derivatives of $\hat{F}^\mu$, $\hat{\mathcal{D}}^\mu_{\hphantom{\mu}\nu}$, etc. Crucially though the Fisher information is controlled -- as expected -- entirely by *first* derivatives of the forward model. Second derivatives of $\hat{F}$ only appear in total-divergence terms that vanish on average.

---

## Reproducing the linear estimator

To make the Fisher formula explicit, we specialize -- inspired by cosmology -- to a forward model that is a function of $x$ rescaled component-by-component:

$$
F^\mu(x, \bm{\theta}) = F^\mu\big(\tau(\bm{\theta})\cdot x\big)\,\,,
$$

where $\tau(\bm{\theta})\cdot x$ is

$$
\tau(\bm{\theta})^\mu\,x^\mu \quad \text{(no sum on } \mu\text{)}\,\,.
$$

In the cosmological application $\tau^\mu$ is a $\vert{\bm k}\vert$-dependent transfer function $\mathcal{M}(k,\bm{\theta})/\mathcal{M}(k,\hat{\bm{\theta})}$, but the algebra is identical in any dimension. Define $\gamma^\mu_i = x^\mu\partial_i \ln \tau^\mu\big\vert_{\hat{\bm{\theta}}}$ (no sum on $\mu$). The key identity is

$$
\partial_i F^\mu = \gamma^\nu_i\, \mathcal{D}^\mu_{\hp{\mu}\nu}\,\,,
$$

i.e. differentiating $F$ with respect to a parameter can be "exchanged" with a derivative w.r.t. $x^\mu$.

Substituting these formulas in terms of the fiducial $\hat{x}^\mu$ into the Fisher formula and using $\langle \hat{x}^\mu \hat{x}^\nu \rangle = C^{\mu\nu}$, all terms combine -- at all orders in $\hat{x}$ -- to give

$$
F_{ii} = 2\sum_\mu \left(\frac{\partial_i \tau^\mu}{\tau^\mu}\right)^2\,\,,
$$

exactly the **optimal linear estimator** result. The conclusion: inference with a nonlinear forward model achieves the same Fisher information as the optimal linear estimator in the zero-noise limit. Crucially, though, we have to recall that this result hinges on the *implicit assumption* that the inverse of the forward model is uniquely defined *everywhere* (since we are integrating over all $\hat{x}^\mu$, once we start computing averages of the log-posterior).

---

## Summary

Three results, each following from the same two moves — a change of variables and integration by parts:

1. **Exact posterior.** In the zero-noise limit, the integral over $x$ can be performed exactly by a change of variables, giving ${-\ln \mathcal{P}(\hat{y} \mid \bm{\theta})}$ as a prior term plus a log-Jacobian. Beyond the zero-noise limit, the result extends order by order in the noise via the cumulant generating function of the likelihood. All $\bm{\theta}$-derivatives of the log-posterior follow from derivatives of $F$ alone — no additional integrals — making score and Fisher information accessible via autodiff through $F$.

2. **Unbiasedness.** The gradient of the log-posterior is a covariant divergence in the space weighted by the prior, so its expectation vanishes by integration by parts. This holds to each order in the noise expansion.

3. **Fisher = linear.** The second-derivative contributions to the Fisher information are total divergences and vanish on average. The remaining terms match the optimal linear estimator exactly for any nonlinear forward model of a rescaled field — provided the forward model is globally invertible.
