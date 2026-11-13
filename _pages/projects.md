---
layout: page
title: ML Projects
permalink: /projects/
---

<ul class="project-list">

  <li class="project-card">
    <div class="project-title">Metric-preserving autoencoders</div>
    <div class="project-desc">Autoencoders trained to preserve distances, not just reconstruct inputs — so two things far apart in reality don't land oddly close together once compressed. The distance-preservation is enforced two different ways: locally, by penalizing the decoder's Jacobian directly at each encoded point (demonstrated on a toy cylinder manifold); and globally, by matching pairwise distances between encoded points to their original high-dimensional distances, weighted toward nearby pairs. Applied to real gravitational-wave template banks, where this lets template-bank coverage be reasoned about directly in the compressed coordinates instead of the original, expensive-to-compare waveform space.</div>
    <div class="project-tags">
      <span>PyTorch</span>
      <span>Riemannian geometry</span>
      <span>Autoencoders</span>
      <span>Gravitational waves</span>
    </div>
    <div class="project-links">
      <a href="https://github.com/giovanni-cabass/autoencoder_GWs">GitHub</a>
      <a href="https://arxiv.org/abs/2511.10466">Paper</a>
    </div>
  </li>

  <li class="project-card">
    <div class="project-title">Geodesics on autoencoder manifolds</div>
    <div class="project-desc">Geodesic distances on the latent manifold induced by a trained decoder — one decoder reused three ways: as a distillation teacher for a separate autoencoder, as the map whose Jacobian defines the geodesic metric, and as the target of a metric-flattening experiment. Applied to a 3-parameter cosmological model compressed into a 16-dimensional representation of the matter power spectrum: geodesics computed by minimizing a discretized path-energy functional, a Ricci-scalar curvature diagnostic to check whether the manifold can even be flattened globally, and a neural-ODE test of whether a learned coordinate change can flatten it back to Euclidean. Companion to the project above.</div>
    <div class="project-tags">
      <span>PyTorch</span>
      <span>Differential geometry</span>
      <span>Autoencoders</span>
      <span>Neural ODEs</span>
    </div>
    <div class="project-links">
      <a href="https://github.com/giovanni-cabass/autoencoder_manifold_geodesics">GitHub</a>
      <a href="https://arxiv.org/abs/2511.10466">Paper</a>
    </div>
  </li>

  <li class="project-card">
    <div class="project-title">TFT for crypto forecasting</div>
    <div class="project-desc">Relatively advanced stage, investigating how to fold physics ideas (e.g. separate-universe picture in cosmology) into addressing long-short timescale couplings. Built on cryptocurrency markets in part for the much higher data cadence available there versus traditional markets — which meant real data engineering on two separate fronts: a training dataset built by unifying crypto exchange, macro, and on-chain data onto one coherent ~8-year timeline, and, separately, a resilient, crash-safe collection pipeline — running unattended for months with zero data loss — for the options data used to develop and test the trading strategy itself. The multi-GPU training infrastructure I built for this project is public — see the link below.</div>
    <div class="project-tags">
      <span>PyTorch</span>
      <span>Temporal Fusion Transformer</span>
      <span>Time series</span>
    </div>
    <div class="project-links">
      <a href="https://github.com/giovanni-cabass/ddp-training-harness">Training infrastructure (on GitHub)</a>
    </div>
  </li>

</ul>
