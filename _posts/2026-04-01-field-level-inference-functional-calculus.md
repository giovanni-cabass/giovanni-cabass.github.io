---
layout: post
title: "Field-level inference and functional calculus: an exact inversion"
date: 2026-04-01
---

This post is an exposition of Appendix A of [arXiv:2307.04706](https://arxiv.org/abs/2307.04706), joint work with Marko Simonović and Matias Zaldarriaga. The appendix contains a self-contained functional-calculus argument that I find elegant enough to deserve its own writeup.

---

## The setup

In *field-level inference* one writes down the full posterior over cosmological parameters $\vec\theta$ given the observed galaxy field $\hat\delta_g$. In the idealized limit of zero noise (cosmic-variance-limited data), the conditional likelihood is a Dirac delta functional,

$$
\mathcal{P}[\hat\delta_g \mid \delta, \vec\theta] = \delta^{(\infty)}_D\!\bigl(\hat\delta_g - F[\delta, \vec\theta]\bigr),
$$

where $F[\delta,\vec\theta]$ is the *forward model*: it takes the primordial initial conditions $\delta(\mathbf{k})$ and outputs the predicted galaxy field. Marginalizing over initial conditions gives the posterior

$$
\mathcal{P}[\hat\delta_g \mid \vec\theta] = \int \mathcal{D}\delta\; \delta^{(\infty)}_D\!\bigl(\hat\delta_g - F[\delta,\vec\theta]\bigr)\, \mathcal{P}[\delta].
$$

The prior $\mathcal{P}[\delta]$ is Gaussian with power spectrum $P(k)$.

The question we want to answer is: **what is the posterior over $\vec\theta$, to all orders in perturbation theory?** The main text of the paper answers this perturbatively, mode by mode. Appendix A gives an exact, non-perturbative answer using functional calculus.

---

## Index notation

The key move is to introduce a compact notation that makes the functional algebra look like finite-dimensional matrix algebra.

- **Greek indices** $(\mu, \nu, \rho, \ldots)$ label momentum modes. The initial conditions $\delta(\mathbf{k})$ become coordinates $\delta^\mu$, the fiducial initial conditions $\hat\delta(\mathbf{k})$ become $\hat\delta^\mu$. Functional derivatives with respect to $\hat\delta(\mathbf{k})$ are written $\partial_\mu$.
- **Latin indices** $(i, j, \ldots)$ label the parameters $\vec\theta$. The derivative $\partial_i$ is always evaluated at the fiducial $\hat{\vec\theta}$.

We define the forward model and its derivatives evaluated at the fiducial point:

$$
\hat F^\mu = F[\hat\delta, \hat{\vec\theta}](\mathbf{k}), \qquad \hat{\mathcal{D}}^\mu_{\ \nu} = \partial_\nu \hat F^\mu, \qquad \hat{\mathcal{D}}^\mu_{\ \nu_1\nu_2\cdots\nu_n} = \partial_{\nu_1}\cdots\partial_{\nu_n} \hat F^\mu.
$$

The matrix $\hat{\mathcal{I}}^\mu_{\ \nu}$ is the functional inverse of $\hat{\mathcal{D}}^\mu_{\ \nu}$:

$$
\hat{\mathcal{I}}^\mu_{\ \rho}\, \hat{\mathcal{D}}^\rho_{\ \nu} = \delta^\mu_{\ \nu}.
$$

In finite dimensions this would just be an ordinary matrix inverse. Here it is the inverse of the Jacobian of the forward model — the object that tells you how a small change in initial conditions propagates to a change in the galaxy field.

For the prior we write

$$
\partial_\mu \ln \hat{\mathcal{P}} = -P_{\mu\nu}\hat\delta^\nu, \qquad \partial_\mu\partial_\nu \ln \hat{\mathcal{P}} = -P_{\mu\nu},
$$

where $P_{\mu\nu} = P^{-1}(k)\,\delta_D(\mathbf{k}+\mathbf{k}')$ is the inverse power spectrum as a matrix.

Finally, we define the "mismatch" between the forward model evaluated at $\vec\theta$ and at the fiducial $\hat{\vec\theta}$:

$$
X^\mu = \hat F^\mu - \tilde F^\mu, \qquad \tilde F^\mu = F[\hat\delta, \vec\theta](\mathbf{k}).
$$

Note that $X^\mu$ vanishes when $\vec\theta = \hat{\vec\theta}$. This will be the key small quantity for the perturbative expansion in $\vec\theta$.

---

## Evaluating the path integral

We shift $\delta^\mu = \hat\delta^\mu + \Delta^\mu$, so the integrand peaks at $\Delta^\mu = 0$, and expand the argument of the delta functional:

$$
F^\mu - \hat F^\mu = \underbrace{\tilde{\mathcal{D}}^\mu_{\ \nu}\Delta^\nu + \frac{1}{2!}\tilde{\mathcal{D}}^\mu_{\ \nu\rho}\Delta^\nu\Delta^\rho + \cdots}_{\equiv\; G[\Delta]^\mu} - X^\mu.
$$

The delta functional forces $G[\Delta]^\mu = X^\mu$. We change variables from $\Delta^\mu$ to $Y^\mu = G[\Delta]^\mu$ — a functional change of variables, with a Jacobian — and the delta functional sets $Y^\mu = X^\mu$. The result is exact:

$$
-\ln \mathcal{P}[\hat\delta_g \mid \vec\theta] = -\ln \mathcal{P}\!\left[\hat\delta + G^{-1}[X]\right] - \ln \left|\frac{\partial G^{-1}[X]}{\partial X}\right|.
$$

Two terms appear: the **prior** evaluated at the inferred initial conditions $\hat\delta + G^{-1}[X]$, and a **Jacobian** term from the change of variables. We only consider the solution connected to linear theory, which is the physical branch for a forward model built from a filtered field.

---

## Unbiasedness via Stokes' theorem

We want to show that the posterior is unbiased: $\langle -\partial_i \ln \mathcal{P} \rangle = 0$. We expand $G^{-1}[X]$ as a power series in $X^\mu$:

$$
G^{-1}[X]^\mu = \tilde{\mathcal{I}}^\mu_{\ \nu} X^\nu - \frac{1}{2} \tilde{\mathcal{I}}^\nu_{\ \beta}\, \tilde{\mathcal{I}}^\mu_{\ \rho}\, \tilde{\mathcal{I}}^\sigma_{\ \alpha}\, \tilde{\mathcal{D}}^\rho_{\ \sigma\nu}\, X^\alpha X^\beta + \cdots,
$$

which uses the identity $\partial_\nu \tilde{\mathcal{I}}^\mu_{\ \alpha} = -\tilde{\mathcal{I}}^\mu_{\ \rho}\, \tilde{\mathcal{D}}^\rho_{\ \sigma\nu}\, \tilde{\mathcal{I}}^\sigma_{\ \alpha}$ (differentiation of a matrix inverse). Since $X^\mu$ vanishes at $\vec\theta = \hat{\vec\theta}$, differentiating with respect to $\theta_i$ and evaluating at the fiducial gives

$$
-\partial_i \ln \mathcal{P}[\hat\delta_g \mid \hat{\vec\theta}] = (\partial_\nu \ln \hat{\mathcal{P}})\, \hat{\mathcal{I}}^\nu_{\ \mu}\, \partial_i \hat F^\mu + \partial_\nu\!\left\{\hat{\mathcal{I}}^\nu_{\ \mu}\, \partial_i \hat F^\mu\right\}.
$$

Now comes the elegant step. Define $V^\nu_i \equiv \hat{\mathcal{I}}^\nu_{\ \mu}\, \partial_i \hat F^\mu$. Then the expression above is

$$
-\partial_i \ln \mathcal{P} = \underbrace{(\partial_\nu \ln \hat{\mathcal{P}})}_{\Gamma^\rho_{\ \rho\nu}} V^\nu_i + \partial_\nu V^\nu_i = \nabla_\nu V^\nu_i,
$$

where $\nabla_\nu$ is a *covariant* divergence. Averaging over $\hat\delta$ is equivalent to integrating over $\hat\delta$ with measure $\hat{\mathcal{P}}[\hat\delta]\, \mathcal{D}\hat\delta$. The term $\partial_\nu \ln \hat{\mathcal{P}}$ plays the role of the Christoffel symbol $\Gamma^\rho_{\ \rho\nu}$ for the metric with determinant $\hat{\mathcal{P}}$. By **Stokes' theorem in functional space**, and using the fact that the Gaussian measure decays exponentially on the "boundary" at $|\hat\delta| \to \infty$,

$$
\left\langle -\partial_i \ln \mathcal{P}[\hat\delta_g \mid \hat{\vec\theta}] \right\rangle = \int \mathcal{D}\hat\delta\; \hat{\mathcal{P}}\, \nabla_\nu V^\nu_i = 0.
$$

This is the non-perturbative proof of unbiasedness.

---

## Second derivatives and Fisher information

The second derivative $-\partial_i\partial_j \ln \mathcal{P}$ has four terms: two from differentiating the prior piece, two from differentiating the Jacobian. Writing $Z^\mu_{\ \nu} = \partial G^{-1}[X]^\mu / \partial X^\nu$, one finds

$$
-\partial_i\partial_j \ln \mathcal{P} = -(\partial_\mu\partial_\nu \ln \hat{\mathcal{P}})\,\partial_i G^{-1}[X]^\mu\, \partial_j G^{-1}[X]^\nu - (\partial_\mu \ln \hat{\mathcal{P}})\,\partial_i\partial_j G^{-1}[X]^\mu + J^{11}_{ij} + J^{02}_{ij},
$$

where $J^{11}$ and $J^{02}$ come from the Jacobian. The explicit expressions for each term are lengthy, but there is a key structural simplification: the two "second-derivative" terms combine into a total covariant divergence,

$$
-(\partial_\mu \ln \hat{\mathcal{P}})\,\partial_i\partial_j G^{-1}[X]^\mu - (Z^{-1})^\nu_{\ \mu}\, \partial_i\partial_j Z^\mu_{\ \nu} = \nabla_\mu\!\left\{\partial_i\partial_j G^{-1}[X]^\mu\right\}.
$$

By the same Stokes' theorem argument, this term vanishes on average. The upshot is that **on average, the Fisher information is controlled entirely by first derivatives of the forward model** — not second derivatives. This is the functional-calculus analogue of the standard Fisher matrix formula $F_{ij} = \langle \partial_i \ln \mathcal{L}\; \partial_j \ln \mathcal{L} \rangle$.

---

## Reproducing linear theory

As a non-trivial check, consider a parameter $\theta$ that enters only through the linear matter power spectrum, via a rescaling

$$
\delta(\mathbf{k}) \to \tau(k, \vec\theta)\, \delta(\mathbf{k}), \qquad \tau(k, \vec\theta) = \frac{\mathcal{M}(k, \vec\theta)}{\mathcal{M}(k, \hat{\vec\theta})}.
$$

Define the field $\gamma(\mathbf{k}) = (\tau'/\tau)\, \delta(\mathbf{k})$ and the matrix $\Gamma(\mathbf{k},\mathbf{k}') = (\tau'/\tau)\, \delta_D(\mathbf{k}-\mathbf{k}')$. The key identity is

$$
\partial_{i_0} F^\mu = \gamma^\nu\, \mathcal{D}^\mu_{\ \nu},
$$

which says: differentiating the forward model with respect to a power-spectrum parameter is the same as acting on the initial conditions with a diagonal rescaling $\gamma$ and then pushing forward through the Jacobian $\mathcal{D}$. This follows from the fact that $\tau$ and $\delta$ always enter the forward model in the combination $\tau\delta$.

Substituting into the Fisher information and using $\langle \delta^\mu \delta^\nu \rangle = P^{\mu\nu}$, all terms combine to give

$$
\left\langle -\partial^2_{i_0} \ln \mathcal{P} \right\rangle = 2 \int \frac{\mathrm{d}^3 k}{(2\pi)^3}\, V \left(\frac{\tau'(k, \hat{\vec\theta})}{\tau(k, \hat{\vec\theta})}\right)^2,
$$

exactly the linear-theory Fisher information. This is the non-perturbative statement that field-level inference with a perturbative forward model does not lose information compared to the optimal linear estimator, **regardless of the order at which the forward model is truncated**.

---

## Summary

The argument in Appendix A achieves three things:

1. **Exact posterior.** The path integral over initial conditions can be done at all orders by the change of variables $G[\Delta] = X$, yielding a closed-form expression for $-\ln \mathcal{P}$.

2. **Unbiasedness.** The gradient of the log-posterior is a covariant divergence in functional space, so its expectation vanishes by Stokes' theorem.

3. **Fisher = linear theory.** Second derivatives of the posterior simplify dramatically — the "second-derivative" terms are total divergences — and the remaining "first-derivative squared" terms reproduce the linear-theory Fisher information exactly.

The functional-calculus notation makes all three steps look like their finite-dimensional analogues. The cost is keeping track of the Jacobians; the reward is results that hold to all orders.
