# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GEMB (Glacier Energy and Mass Balance Model) is a 1-dimensional physical model for simulating surface energy balance and vertical firn evolution of glaciers and ice sheets. It resolves temperature, density, water content, and grain properties over time. Written in MATLAB; Octave compatibility is untested.

**Citation**: Gardner et al., 2023 — Geosci. Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023

## Running Tests

Tests use MATLAB's built-in unit testing framework. CI runs them via GitHub Actions (`.github/workflows/main.yml`) against the `/test` folder recursively.

To run all tests locally in MATLAB:
```matlab
results = runtests('test', 'RecurseSubfolders', true);
```

To run a single test file:
```matlab
results = runtests('test/test_calculate_albedo.m');
```

## Basic Model Workflow

The model requires four steps in sequence:

```matlab
addpath("src")

% 1. Define climate forcing
ClimateForcing = model_initialize_forcing(...);
% or use synthetic data:
ClimateForcing = simulate_climate_forcing("test_1", 3);  % 3-hourly

% 2. Configure model parameters
ModelParam = model_initialize_parameters();

% 3. Initialize firn column
Profile = model_initialize_profile(ModelParam, ClimateForcing);

% 4. Run the model
OutData = gemb(Profile, ClimateForcing, ModelParam);
```

See `examples/` for complete working examples (synthetic data and ERA5 reanalysis).

## Architecture

### Call Hierarchy

```
gemb()                          ← main driver; loops over time
  └── gemb_core()               ← single time-step physics engine
        ├── calculate_accumulation()
        ├── calculate_albedo()
        ├── calculate_shortwave_radiation()
        ├── calculate_temperature()   ← heat equation with phase change
        │     ├── thermal_conductivity()
        │     └── turbulent_heat_flux()
        ├── calculate_melt()          ← percolation and refreezing
        ├── calculate_density()       ← Herron-Langway / Arthern / Ligtenberg
        ├── calculate_grain_size()
        └── manage_layers()           ← Lagrangian grid merging/splitting
```

### Key Data Structures

- **`ClimateForcing`** — MATLAB timetable with time-series surface forcing (temperature, pressure, precipitation, shortwave/longwave radiation, wind speed, vapor pressure)
- **`ModelParam`** — struct with ~38 configuration fields (densification method, albedo scheme, grid geometry, output controls); initialized by `model_initialize_parameters()`
- **`Profile`** — struct representing the initial firn column state (temperature, density, layer thicknesses, grain properties, water content)
- **`OutData`** — struct of time-series output arrays for the full model run

### Directory Layout

| Path | Contents |
|------|----------|
| `src/` | All model source code |
| `src/fit_simulated_climate_to_data/` | Functions for fitting synthetic climate to observations |
| `examples/` | Example scripts (synthetic and ERA5 forcing) |
| `test/` | Unit tests (one file per `src/` function) |
| `docs/` | Markdown documentation for every public function |

### Key Source Files

- `src/gemb.m` — main driver; manages time-stepping and output collection
- `src/gemb_core.m` — single time-step integration; orchestrates all physical process functions
- `src/manage_layers.m` — dynamic vertical grid management (merge/split layers to maintain thickness constraints)
- `src/model_initialize_parameters.m` — defines all tunable model parameters with defaults
- `src/model_initialize_profile.m` — builds the initial column state
- `src/model_initialize_forcing.m` — constructs a timetable from raw climate data arrays

### Spinup Workflow

For multi-millennial spinups, use `gemb_profile.m` to extract the final column state from one run as the initial profile for the next, and `gemb_interp.m` to regularize output onto a consistent vertical grid.

## Documentation

Every public function has a corresponding Markdown doc in `docs/`. The variable reference (`docs/GEMB_variables.md`) lists all 86+ model variables with units and descriptions. `docs/GEMB_overview.md` contains the function dependency diagram.
