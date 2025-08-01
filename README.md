# priestley-taylor

[![CI](https://github.com/gregory-halverson-jpl/priestley-taylor/actions/workflows/ci.yml/badge.svg)](https://github.com/gregory-halverson-jpl/priestley-taylor/actions/workflows/ci.yml)

The `priestley-taylor` Python package provides robust, peer-reviewed utilities for physically-based evapotranspiration modeling, focusing on the Priestley-Taylor and Penman-Monteith equations. It includes scientific constants and functions for calculating the slope of the saturation vapor pressure curve, the psychrometric constant, and the epsilon parameter, all essential for hydrology, meteorology, and agricultural science applications.

[Gregory H. Halverson](https://github.com/gregory-halverson-jpl) (they/them)<br>
[gregory.h.halverson@jpl.nasa.gov](mailto:gregory.h.halverson@jpl.nasa.gov)<br>
NASA Jet Propulsion Laboratory 329G

## Installation

This package is available on PyPI as `priestley-taylor`.

```
pip install priestley-taylor
```

## Usage

Import this package as `priestley_taylor`:

```python
import priestley_taylor
```

### 1. `GAMMA_KPA` and `GAMMA_PA`
- **Description:** Psychrometric constant (γ) in kPa/°C and Pa/°C, respectively. γ quantifies the relationship between vapor pressure and temperature, accounting for atmospheric pressure and the specific heat of air. It is a key parameter in energy balance approaches to evapotranspiration.
- **Reference:** Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). FAO Irrigation and Drainage Paper 56, Table 2.2.

### 2. `delta_kPa_from_Ta_C(Ta_C)`
- **Description:** Calculates the slope of the saturation vapor pressure curve (Δ) at a given air temperature (°C), in kPa/°C. Δ quantifies the sensitivity of saturation vapor pressure to temperature changes, a key parameter in Priestley-Taylor and Penman-Monteith equations.
- **Parameters:** `Ta_C` (numpy array or Raster): Air temperature in Celsius.
- **Returns:** Slope of saturation vapor pressure curve (kPa/°C).
- **References:**
	- Allen et al. (1998), Eq. 2.18
	- Monteith, J. L. (1965). Evaporation and environment. In The State and Movement of Water in Living Organisms (pp. 205–234). Academic Press.

### 3. `delta_Pa_from_Ta_C(Ta_C)`
- **Description:** Converts Δ from kPa/°C to Pa/°C for use in models requiring Pascals.
- **Parameters:** `Ta_C` (numpy array or Raster): Air temperature in Celsius.
- **Returns:** Slope of saturation vapor pressure curve (Pa/°C).
- **Reference:** Allen et al. (1998).

### 4. `calculate_epsilon(delta, gamma)`
- **Description:** Computes the dimensionless ratio ε (epsilon), defined as ε = Δ / (Δ + γ), representing the relative importance of energy supply versus atmospheric demand in controlling evapotranspiration.
- **Parameters:**
	- `delta` (numpy array or Raster): Slope of saturation vapor pressure curve (Pa/°C or kPa/°C)
	- `gamma` (numpy array or Raster): Psychrometric constant (same units as delta)
- **Returns:** Epsilon (dimensionless).
- **References:**
	- Priestley, C. H. B., & Taylor, R. J. (1972). Monthly Weather Review, 100(2), 81–92.
	- Allen et al. (1998), Eq. 6.2

### 5. `epsilon_from_Ta_C(Ta_C, gamma_Pa=GAMMA_PA)`
- **Description:** Calculates ε from air temperature (°C) and the psychrometric constant (Pa/°C). This combines the calculation of Δ (in Pa/°C) and ε, providing a direct way to obtain the key parameter for the Priestley-Taylor equation.
- **Parameters:**
	- `Ta_C` (numpy array or Raster): Air temperature in Celsius.
	- `gamma_Pa` (float, numpy array, or Raster): Psychrometric constant in Pa/°C (default: GAMMA_PA).
- **Returns:** Epsilon (dimensionless).
- **References:**
	- Priestley, C. H. B., & Taylor, R. J. (1972). Monthly Weather Review, 100(2), 81–92.
	- Allen et al. (1998)

## References

- Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration – Guidelines for computing crop water requirements – FAO Irrigation and Drainage Paper 56. FAO, Rome. https://www.fao.org/3/x0490e/x0490e00.htm
- Monteith, J. L. (1965). Evaporation and environment. In The State and Movement of Water in Living Organisms (pp. 205–234). Academic Press.
- Priestley, C. H. B., & Taylor, R. J. (1972). On the assessment of surface heat flux and evaporation using large-scale parameters. Monthly Weather Review, 100(2), 81–92. https://doi.org/10.1175/1520-0493(1972)100<0081:OTAOSH>2.3.CO;2
