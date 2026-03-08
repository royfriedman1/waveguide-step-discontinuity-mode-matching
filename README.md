# Waveguide Discontinuity Analysis — Mode-Matching Method

**Course:** Microwave Engineering (4th year, Electrical Engineering, Tel Aviv University)

## Overview

This project analyzes electromagnetic wave propagation through a **rectangular waveguide step discontinuity** using the **Mode-Matching Method (MMM)**. The discontinuity connects two waveguide sections with different heights and dielectric fillings.

### Problem Setup

| Parameter | Region 1 (Input) | Region 2 (Output) |
|-----------|-----------------|-------------------|
| Height    | H₁ = 7.9 mm     | H₂ = 10.6 mm      |
| Relative permittivity | εᵣ = 2 | εᵣ = 1 |

Frequency range: ~17.3 – 25 GHz (above TE₁ cutoff in both regions)

## Method

The Mode-Matching Method expands the fields in each region as a superposition of waveguide modes. Boundary conditions at the discontinuity yield a linear system for the scattering coefficients. Key steps:

1. **Compute transverse wavenumbers** for TE modes in each region
2. **Build the overlap integral matrix** I(m,n) — inner products of mode profiles across the aperture
3. **Assemble and solve** the 2N×2N linear system for reflection/transmission amplitudes
4. **Extract S-parameters**: S₁₁ (reflection), S₂₁ (transmission)
5. **Verify power conservation**: |S₁₁|² + |S₂₁|² ≈ 1

N = 100 modes are used for convergence.

## Files

| File | Description |
|------|-------------|
| `WG_ex.m` | MATLAB simulation — mode-matching solver, S-parameter computation, power check |
| `wg_ex.cst` | CST Microwave Studio simulation file (full-wave reference) |
| `תרגיל מחשב.pdf` | Full written report (Hebrew) with derivations, results, and analysis |
| `Computer Exercise.docx` | Written report (English/mixed) |

## Results

The simulation produces:
- **S-parameter magnitude plots** (|S₁₁|, |S₂₁| in dB) vs. frequency
- **Power conservation check** — deviation from unity quantifies numerical error
- Comparison with CST full-wave simulation validates the mode-matching results

## How to Run

Open `WG_ex.m` in MATLAB and run. No additional toolboxes are required.

```matlab
% Key outputs:
% - Figure 1: |S11| and |S21| in dB vs frequency
% - Figure 2: Power conservation (|S11|² + |S21|²)
% - Figure 3: Deviation from power conservation
```

## Technologies

- MATLAB (numerical simulation)
- CST Microwave Studio (full-wave EM solver)
