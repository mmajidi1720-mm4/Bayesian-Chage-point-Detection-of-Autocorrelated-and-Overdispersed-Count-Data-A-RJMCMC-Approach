# Bayesian Change-Point Detection in Autocorrelated and Overdispersed Count Data
> **Supplementary Code Repository**  
> Companion code for the manuscript:  
> *"Bayesian Chage-point Detection of Autocorrelated and Overdispersed Count Data: A Reversible Jump MCMC Approach with Application to COVID-19 Dynamics"*.

---

## 📌 Overview

This repository provides a self-contained, standalone implementation of the **Negative Binomial Autoregressive (NB-AR)** Bayesian change-point model using **Trans-dimensional Reversible Jump Markov Chain Monte Carlo (RJMCMC)**. 

The codebase reproduces the Monte Carlo simulation study (Section 3 of the paper), comparing the proposed **NB-AR** model against the standard **Poisson-IID** change-point baseline under:
1. **Autocorrelated counts** ($\gamma = 0.6$)
2. **Non-autocorrelated counts** ($\gamma = 0.0$)

The implementation is written strictly in **pure base-R (version $\ge$ 4.0.0)** without external compiled dependencies, ensuring cross-platform stability and long-term reproducibility.

---

## 🔬 Methodological Framework

### 1. Data-Generating Process & Likelihood
Let $Y = \{y_1, \dots, y_n\}$ denote a time series of counts partitioned by $k$ change-points $\boldsymbol{\tau} = (\tau_1, \dots, \tau_k)^T$ into $K = k + 1$ distinct regimes, with boundary conventions $\tau_0 = 0$ and $\tau_K = n$. Conditional on past observations $\mathcal{F}_{t-1}$, observations in segment $s \in \{1, \dots, K\}$ satisfy:

$$y_t \mid \mathcal{F}_{t-1}, \boldsymbol{\theta}_s \sim \text{NB}(\mu_{t,s}, \phi_s), \quad \tau_{s-1} < t \le \tau_s$$

where $\phi_s > 0$ denotes the segment-specific dispersion parameter, and the conditional mean is modeled via a log-linear autoregressive specification:

$$\log(\mu_{t,s}) = \beta_{0,s} + \gamma_s \log(y_{t-1} + 1)$$

* The autoregressive parameter $\gamma_s$ captures local serial dependence (momentum).
* Stationarity within each segment is enforced via the truncated prior $\gamma_s \sim \mathcal{N}(0, \sigma^2_\gamma)\,\mathbb{I}(-1 < \gamma_s < 1)$.

### 2. Prior Formulations
* **Number of segments:** Truncated Poisson prior $K \sim \text{Poisson}(\lambda_K = 3)$ for $1 \le K \le K_{\max}$.
* **Change-point locations:** Discrete uniform over valid ordered partitions: $p(\boldsymbol{\tau} \mid K) \propto 1$, subject to minimum segment length constraint $n_{\min}$ ($\tau_s - \tau_{s-1} \ge n_{\min}$). The exact combinatorial normalizer is evaluated via stars-and-bars:
  
  $$\log C(K, n, n_{\min}) = \log \binom{n - K \cdot n_{\min} + K - 1}{K - 1}$$

* **Segment parameters:** Weakly informative priors:
  $$\beta_{0,s} \sim \mathcal{N}(0, 10^2), \quad \gamma_s \sim \mathcal{N}(0, 1^2)\,\mathbb{I}(-1 < \gamma_s < 1), \quad \log(\phi_s) \sim \mathcal{N}(0, \sigma^2_\phi)$$

### 3. Trans-Dimensional RJMCMC Mechanics
The sampler alternates between **Between-Model moves** and **Within-Model moves**:

* **Birth Proposal ($K \to K+1$):**
  A segment $s$ is selected uniformly at random and split at candidate point $\tau^* \sim \mathcal{U}\{\tau_{s-1} + n_{\min}, \dots, \tau_s - n_{\min}\}$. Auxiliary variables $\boldsymbol{u} = (u_\beta, u_\gamma, u_\phi)^T \sim \mathcal{N}(\mathbf{0}, \sigma_u^2 \mathbf{I})$ generate split parameters via the deterministic bijection:
  
  $$\boldsymbol{\theta}_{s,1}' = \boldsymbol{\theta}_s - \boldsymbol{u}, \quad \boldsymbol{\theta}_{s,2}' = \boldsymbol{\theta}_s + \boldsymbol{u}$$

  The transformation decomposes into $p+2 = 3$ uncoupled linear systems with Jacobian determinant:
  
  $$|\mathcal{J}| = 2^1 \times 2 \times 2 = 2^3 = 8$$

* **Metropolis-Hastings-Green Acceptance Probability:**
  $$\alpha_{\text{birth}} = \min\left(1, \, \frac{\mathcal{L}(\boldsymbol{\Theta}' \mid \mathbf{y})}{\mathcal{L}(\boldsymbol{\Theta} \mid \mathbf{y})} \times \frac{p(K+1)\,p(\boldsymbol{\tau}' \mid K+1)\,p(\boldsymbol{\Theta}')}{p(K)\,p(\boldsymbol{\tau} \mid K)\,p(\boldsymbol{\Theta})} \times \frac{d_{K+1}}{b_K} \frac{N_s}{q(\boldsymbol{u})} \times |\mathcal{J}|\right)$$

* **Death Proposal ($K \to K-1$):**
  An active change-point is removed at random. Detailed balance is guaranteed via the exact inverse ratio:
  
  $$A_{\text{death}} = A_{\text{birth}}^{-1}$$

* **Within-Model Proposals:**
  * **Parameter Updates:** Random-walk Metropolis updates for $(\beta_{0,s}, \gamma_s, \log \phi_s)$.
  * **Relocation:** Perturbation of change-point location $\tau_j$ within local window $[\tau_j - h, \tau_j + h]$ respecting $n_{\min}$ boundaries.

---

## 📂 Repository Contents
```text
.
├── simulation_nbar.R       # Complete simulation script and RJMCMC sampling engine
├── LICENSE                 # Open-source MIT License
└── README.md               # Methodological documentation and reproduction guide

---

## 🚀 Getting Started

### Prerequisites
* **R (version $\ge$ 4.0.0)**.
* No additional packages or compiled toolchains (e.g., Rcpp, Stan, JAGS) are required. All computations run in native base-R.

### Execution

#### 1. Quick Verification Run (Default: $N_{\text{sim}} = 20$)
To verify the algorithmic pipeline and inspect acceptance diagnostics:
bash
Rscript simulation_nbar.R
*Approximate runtime: ~10–15 minutes on a modern 8-core CPU.*

#### 2. Full Monte Carlo Study ($N_{\text{sim}} = 100$)
To reproduce the exact published numbers from Table 1 of the manuscript:
1. Open `simulation_nbar.R`.
2. Set `N_SIM <- 100` in **Section 1 (Line 23)**.
3. Run the script:
bash
Rscript simulation_nbar.R
*Approximate runtime: ~1–2 hours depending on CPU performance.*

---

## 📊 Benchmark Results

Reproducing **Table 1 (Section 3.3)** from the manuscript:

| Scenario | Model | Correct Detection Rate (CDR % for $K=2$) | Expected Segments $E[\hat{K}]$ | MAE of $\hat{\tau}$ |
| :--- | :--- | :---: | :---: | :---: |
| **1. Autocorrelated** ($\gamma = 0.6$) | **Poisson-IID** | **0%** | **8.00** | **NA** |
| | **NB-AR (Proposed)** | **75%** | **2.40** | **0.3** |
| **2. Non-Autocorrelated** ($\gamma = 0.0$) | **Poisson-IID** | 80% | 2.35 | 0.0 |
| | **NB-AR (Proposed)** | **100%** | **2.00** | **0.0** |

### Key Findings
* **Unmodeled Serial Correlation:** When counts exhibit serial dependence ($\gamma = 0.6$), the baseline Poisson model fails completely ($\text{CDR} = 0\%$), saturating the state space at $E[\hat{K}] = 8.00$ ($K_{\max}$). Standard models mistake short-term momentum for structural breaks.
* **Proposed Robustness:** The proposed **NB-AR** model achieves a **75% CDR** with high localization accuracy ($\text{MAE} = 0.3$) under strong autocorrelation, and **100% CDR** when observations are conditionally independent.

---

## ⚙️ Key Simulation Parameters

All simulation settings are controlled in **Section 1** of `simulation_nbar.R`:

| Parameter | Default Value | Description |
| :--- | :---: | :--- |
| `N_SIM` | `20` (set to `100` for paper) | Monte Carlo replications per scenario |
| `T_LEN` | `400` | Length of time series |
| `TRUE_TAU` | `200` | True change-point location |
| `TRUE_K` | `2` | True number of segments ($K_{\text{true}} = 2 \implies 1$ change-point) |
| `N_MCMC` | `15,000` | Total MCMC iterations per chain |
| `N_BURN` | `5,000` | Burn-in iterations discarded |
| `THIN` | `5` | Thinning interval for posterior draws |
| `K_MAX` | `8` | Maximum allowable number of segments |
| `N_MIN` | `10` | Minimum allowed segment length |
| `SIGMA_U` | `0.10` | Proposal SD for dimension-matching auxiliary variables |
| `TUNING_SCALE` | `0.05` | Random Walk scale for within-model parameter updates |

---

