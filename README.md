# Priestley-Taylor Potential Latent Heat Flux Python Package

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

## Scientific Methodology

### Theoretical Foundation

The Priestley-Taylor equation provides a simplified yet physically-based approach to estimating potential evapotranspiration (PET) under energy-limited conditions. This method eliminates the need for wind speed and vapor pressure deficit measurements required by the more complex Penman-Monteith equation, making it particularly valuable for regional-scale applications and data-sparse environments.

### Core Equation

The Priestley-Taylor equation for potential latent heat flux is:

**LE = α × ε × (Rn - G)**

Where:
- **LE** = Latent heat flux [W/m²] - energy used for evapotranspiration
- **α** = Priestley-Taylor coefficient [dimensionless] - typically 1.26 for well-watered vegetation
- **ε** = Energy partitioning coefficient [dimensionless] - fraction of available energy converted to latent heat under equilibrium
- **Rn** = Net radiation [W/m²] - total radiative energy available at the surface
- **G** = Soil heat flux [W/m²] - energy conducted into or out of the soil

### Step-by-Step Methodology

#### Step 1: Calculate Saturation Vapor Pressure Slope (Δ)

The slope of the saturation vapor pressure curve quantifies how rapidly the atmosphere's water-holding capacity increases with temperature:

**Δ = 4098 × [0.6108 × exp(17.27 × Ta / (237.3 + Ta))] / (Ta + 237.3)²**

Where:
- **Δ** = Slope of saturation vapor pressure curve [kPa/°C or Pa/°C]
- **Ta** = Air temperature [°C]
- **4098** = Derived coefficient from differentiation [dimensionless]
- **0.6108** = Magnus-Tetens coefficient [kPa]
- **17.27** = Magnus-Tetens temperature coefficient [dimensionless]
- **237.3** = Magnus-Tetens temperature offset [°C]

*Scientific rationale:* This relationship is derived from the Clausius-Clapeyron equation combined with the Magnus-Tetens approximation, providing accurate vapor pressure calculations for temperatures from -40°C to +50°C (Murray, 1967).

#### Step 2: Apply Psychrometric Constant (γ)

The psychrometric constant relates vapor pressure changes to temperature in the context of evaporation:

**γ = (cp × P) / (ε × λ) = 0.0662 [kPa/°C] at standard conditions**

Where:
- **γ** = Psychrometric constant [kPa/°C or Pa/°C]
- **cp** = Specific heat of moist air ≈ 1.013 [kJ/(kg·°C)]
- **P** = Atmospheric pressure [kPa] - 101.3 kPa at sea level
- **ε** = Ratio of molecular weights (water vapor/dry air) = 0.622 [dimensionless]
- **λ** = Latent heat of vaporization ≈ 2.45 [MJ/kg] at 20°C

*Scientific rationale:* This parameter represents the fundamental thermodynamic relationship between sensible and latent heat in the atmosphere, derived from the kinetic theory of gases (Monteith & Unsworth, 2013).

#### Step 3: Compute Energy Partitioning Coefficient (ε)

The epsilon ratio determines how available energy is partitioned between sensible and latent heat:

**ε = Δ / (Δ + γ)**

Where:
- **ε** = Dimensionless ratio [0 < ε < 1]
- **Δ** = Vapor pressure slope [same units as γ]
- **γ** = Psychrometric constant [same units as Δ]

*Scientific rationale:* This ratio represents the theoretical maximum fraction of available energy that can be converted to evapotranspiration under equilibrium conditions, based on atmospheric thermodynamics (Priestley & Taylor, 1972).

#### Step 4: Apply Priestley-Taylor Coefficient (α)

The empirical coefficient accounts for non-equilibrium conditions and surface characteristics:

**Typical values:**
- **α = 1.26** for well-watered vegetation (Priestley & Taylor, 1972)
- **α = 1.05-1.15** for open water surfaces
- **α = 0.8-1.2** for water-stressed vegetation
- **α = 1.3-1.5** for advective conditions (desert oases)

*Scientific rationale:* This coefficient compensates for the aerodynamic and surface resistance terms that are explicitly calculated in the Penman-Monteith equation, allowing the simplified energy-balance approach to maintain accuracy (Stewart & Rouse, 1977).

#### Step 5: Energy Balance Components

**Net Radiation (Rn):**
Rn = (SWin - SWout) + (LWin - LWout) = SWin(1-α) + εσ(Ta⁴ - Ts⁴)

Where:
- **SWin** = Incoming shortwave radiation [W/m²]
- **α** = Surface albedo [dimensionless, 0-1]
- **ε** = Surface emissivity [dimensionless, 0.85-0.99]
- **σ** = Stefan-Boltzmann constant = 5.67×10⁻⁸ [W/(m²·K⁴)]
- **Ta, Ts** = Air and surface temperatures [K]

**Soil Heat Flux (G):**
Typically estimated as G ≈ 0.1 × Rn for daily averages, or calculated using empirical relationships with vegetation indices and surface temperature.

### Physical Interpretation

The Priestley-Taylor approach assumes that under well-watered conditions:

1. **Energy limitation dominates** over atmospheric demand
2. **Surface resistance is minimal** (stomata fully open)
3. **Horizontal advection effects** are captured by the α coefficient
4. **Equilibrium conditions** approximate real evapotranspiration

This makes the method particularly suitable for:
- Dense vegetation canopies with adequate water supply
- Regional-scale applications where local advection averages out
- Climate change impact assessments
- Irrigation scheduling under non-limiting water conditions

### Advantages and Limitations

**Advantages:**
- Requires fewer meteorological inputs than Penman-Monteith
- More stable in data-sparse regions
- Less sensitive to measurement uncertainties
- Computationally efficient for large-scale applications

**Limitations:**
- Cannot account for stomatal control under water stress
- May overestimate ET in arid environments
- Less accurate under strong advective conditions
- Assumes sufficient water availability

## Usage

Import this package as `priestley_taylor`:

```python
import priestley_taylor
```

### 1. `GAMMA_KPA` and `GAMMA_PA`
- **Description:** Psychrometric constant (γ) representing the fundamental thermodynamic relationship between vapor pressure and temperature in atmospheric evaporation processes.
- **Values:** 
  - `GAMMA_KPA = 0.0662` [kPa/°C] - Standard conditions (101.3 kPa, 20°C)
  - `GAMMA_PA = 66.2` [Pa/°C] - SI unit equivalent
- **Physical meaning:** Ratio of specific heat of moist air to latent heat of vaporization, multiplied by the molecular weight ratio of water vapor to dry air
- **Formula:** γ = (cp × P) / (ε × λ)
  - cp = specific heat of moist air ≈ 1.013 kJ/(kg·°C)
  - P = atmospheric pressure [kPa]
  - ε = molecular weight ratio = 0.622 [dimensionless]
  - λ = latent heat of vaporization ≈ 2.45 MJ/kg
- **Elevation dependency:** γ(P) = 0.0662 × (P/101.3) for pressure corrections
- **References:** 
  - Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). FAO Irrigation and Drainage Paper 56, Table 2.2
  - Monteith, J. L., & Unsworth, M. H. (2013). Principles of Environmental Physics, 4th Ed.

### 2. `delta_kPa_from_Ta_C(Ta_C)`
- **Description:** Calculates the slope of the saturation vapor pressure curve (Δ) using the Magnus-Tetens approximation, representing how rapidly atmospheric water-holding capacity increases with temperature.
- **Formula:** Δ = 4098 × [0.6108 × exp(17.27 × Ta / (237.3 + Ta))] / (Ta + 237.3)²
- **Parameters:** 
  - `Ta_C` (numpy array or Raster): Air temperature [°C], valid range: -40 to +50°C
- **Returns:** Slope of saturation vapor pressure curve [kPa/°C], always positive
- **Physical significance:** Higher temperatures exponentially increase evaporative demand
- **Typical values:**
  - 0°C: Δ ≈ 0.189 kPa/°C
  - 20°C: Δ ≈ 1.45 kPa/°C  
  - 40°C: Δ ≈ 3.99 kPa/°C
- **Applications:** Central to all combination evapotranspiration equations (Penman, Priestley-Taylor, Penman-Monteith)
- **References:**
  - Allen et al. (1998), Equation 2.18 - FAO-56 standard methodology
  - Magnus, G. (1844). Annalen der Physik, 137(12), 225-247 - Original vapor pressure formulation
  - Murray, F. W. (1967). Journal of Applied Meteorology, 6(1), 203-204 - Accuracy validation

### 3. `delta_Pa_from_Ta_C(Ta_C)`
- **Description:** Calculates the slope of saturation vapor pressure curve in SI base units (Pascals) for consistency with meteorological modeling systems.
- **Unit conversion:** Δ[Pa/°C] = Δ[kPa/°C] × 1000 (exact conversion)
- **Parameters:** 
  - `Ta_C` (numpy array or Raster): Air temperature [°C]
- **Returns:** Slope of saturation vapor pressure curve [Pa/°C]
- **Typical values:**
  - 0°C: Δ ≈ 189 Pa/°C
  - 20°C: Δ ≈ 1450 Pa/°C
  - 40°C: Δ ≈ 3990 Pa/°C
- **Applications:** Preferred for models using Pascal pressure units throughout calculations
- **References:**
  - International System of Units (SI), BIPM (2019)
  - Allen et al. (1998) - Methodological foundation

### 4. `calculate_epsilon(delta, gamma)`
- **Description:** Computes the dimensionless energy partitioning coefficient (ε) that determines the theoretical fraction of available energy converted to latent heat under equilibrium conditions.
- **Formula:** ε = Δ / (Δ + γ)
- **Parameters:**
  - `delta` (numpy array or Raster): Vapor pressure slope [Pa/°C or kPa/°C]
  - `gamma` (numpy array or Raster): Psychrometric constant [same units as delta]
- **Returns:** Epsilon [dimensionless], range: 0 < ε < 1
- **Physical interpretation:** 
  - ε → 0: Cold conditions, minimal evaporative efficiency
  - ε → 1: Hot conditions, maximum energy conversion to ET
  - Typical range: 0.4-0.99 for terrestrial conditions
- **Temperature dependence:** Increases exponentially with temperature due to Δ growth
- **Applications:** 
  - Priestley-Taylor: LE = α × ε × (Rn - G)
  - Equilibrium evaporation: LEeq = ε × (Rn - G)
  - Decoupling coefficient calculations
- **References:**
  - Priestley, C. H. B., & Taylor, R. J. (1972). Monthly Weather Review, 100(2), 81-92
  - Penman, H. L. (1948). Proceedings of the Royal Society, 193(1032), 120-145
  - Jarvis, P. G., & McNaughton, K. G. (1986). Advances in Ecological Research, 15, 1-49

### 5. `epsilon_from_Ta_C(Ta_C, delta_Pa=None, gamma_Pa=GAMMA_PA)`
- **Description:** Direct calculation of energy partitioning coefficient from air temperature, combining vapor pressure slope calculation with epsilon computation for computational efficiency.
- **Process flow:** Ta[°C] → Δ[Pa/°C] → ε[dimensionless]
- **Parameters:**
  - `Ta_C` (numpy array or Raster): Air temperature [°C], valid range: -40 to +50°C
  - `delta_Pa` (optional): Pre-computed vapor pressure slope [Pa/°C] for efficiency
  - `gamma_Pa` (float/array): Psychrometric constant [Pa/°C], default: 66.2 Pa/°C
- **Returns:** Epsilon [dimensionless], typical range: 0.4-0.99
- **Example values:**
  - 0°C: ε ≈ 0.741
  - 20°C: ε ≈ 0.956
  - 40°C: ε ≈ 0.984
- **Temperature sensitivity:** Non-linear relationship, highest sensitivity at 0-10°C
- **Applications:** Direct input to Priestley-Taylor equation for potential ET estimation
- **References:**
  - Priestley, C. H. B., & Taylor, R. J. (1972). Monthly Weather Review, 100(2), 81-92
  - Shuttleworth, W. J. (1993). Handbook of Hydrology, McGraw-Hill, Ch. 4

### 6. `priestley_taylor(...)`
- **Description:** Main function implementing the complete Priestley-Taylor methodology for potential evapotranspiration estimation with automatic data retrieval and energy balance calculations.
- **Core equation:** LE = α × ε × (Rn - G)
- **Key parameters:**
  - `Ta_C`: Air temperature [°C] - drives vapor pressure relationships
  - `Rn_Wm2`: Net radiation [W/m²] - primary energy source
  - `G_Wm2`: Soil heat flux [W/m²] - energy not available for ET
  - `PT_alpha`: Priestley-Taylor coefficient [dimensionless], default: 1.26
  - `ST_C`, `albedo`, `emissivity`: Surface properties for radiation calculations
  - `NDVI`: Vegetation index for soil heat flux estimation
- **Automatic calculations:**
  - Net radiation using Verma et al. method when surface properties provided
  - Soil heat flux using SEBAL algorithm when vegetation data available
  - Meteorological data retrieval from GEOS-5 FP when geometry/time specified
- **Returns:** Dictionary with `LE_potential_Wm2`, `epsilon`, `G_Wm2`, and intermediate variables
- **Applications:**
  - Regional evapotranspiration mapping
  - Irrigation water requirements
  - Climate impact assessment
  - Hydrological modeling
- **References:**
  - Priestley & Taylor (1972) - Core methodology
  - Verma et al. (1976) - Net radiation algorithm
  - Bastiaanssen et al. (1998) - SEBAL soil heat flux method

## Example Usage

```python
import numpy as np
from datetime import datetime
import priestley_taylor as pt

# Basic calculation with provided inputs
Rn = np.array([[400, 500], [450, 550]])  # Net radiation [W/m²]
G = np.array([[40, 50], [45, 55]])       # Soil heat flux [W/m²]  
Ta = np.array([[20, 25], [22, 27]])      # Air temperature [°C]

# Calculate potential evapotranspiration
result = pt.priestley_taylor(Rn_Wm2=Rn, G_Wm2=G, Ta_C=Ta)
LE_potential = result["LE_potential_Wm2"]
print(f"Potential ET: {LE_potential} W/m²")

# Advanced usage with automatic data retrieval
from rasters import RasterGeometry

# Define spatial domain and time
geometry = RasterGeometry.from_bounds(...)  # Your study area
time_utc = datetime(2024, 6, 15, 12, 0)    # Noon on June 15

# Calculate with satellite/reanalysis inputs
result = pt.priestley_taylor(
    geometry=geometry,
    time_UTC=time_utc,
    ST_C=surface_temperature,     # From thermal remote sensing
    albedo=surface_albedo,        # From optical remote sensing  
    emissivity=surface_emissivity, # From spectral indices
    NDVI=vegetation_index         # From red/NIR bands
)
```

## References

### Foundational Papers
- **Priestley, C. H. B., & Taylor, R. J. (1972).** On the assessment of surface heat flux and evaporation using large-scale parameters. *Monthly Weather Review*, 100(2), 81-92. https://doi.org/10.1175/1520-0493(1972)100<0081:OTAOSH>2.3.CO;2
  - *Original derivation and validation of the Priestley-Taylor equation*

- **Penman, H. L. (1948).** Natural evaporation from open water, bare soil and grass. *Proceedings of the Royal Society of London Series A*, 193(1032), 120-145. https://doi.org/10.1098/rspa.1948.0037
  - *Theoretical foundation for combination evapotranspiration equations*

### Methodological Standards
- **Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).** Crop evapotranspiration: Guidelines for computing crop water requirements. *FAO Irrigation and Drainage Paper 56*. FAO, Rome. https://www.fao.org/3/x0490e/x0490e00.htm
  - *International standard for evapotranspiration calculations and reference methodology*

- **Monteith, J. L., & Unsworth, M. H. (2013).** Principles of Environmental Physics: Plants, Animals, and the Atmosphere. 4th Edition. Academic Press. ISBN: 978-0-12-386910-4
  - *Comprehensive treatment of environmental physics and energy balance principles*

### Vapor Pressure Relationships
- **Magnus, G. (1844).** Versuche über die Spannkräfte des Wasserdampfs. *Annalen der Physik*, 137(12), 225-247. https://doi.org/10.1002/andp.18441371202
  - *Original formulation of vapor pressure-temperature relationships*

- **Murray, F. W. (1967).** On the computation of saturation vapor pressure. *Journal of Applied Meteorology*, 6(1), 203-204. https://doi.org/10.1175/1520-0450(1967)006<0203:OTCOSV>2.0.CO;2
  - *Accuracy assessment of Magnus-Tetens approximation*

### Validation and Applications
- **Stewart, R. B., & Rouse, W. R. (1977).** Substantiating the Priestley and Taylor parameter α = 1.26 for potential evaporation in high latitudes. *Journal of Applied Meteorology*, 16(6), 649-650. https://doi.org/10.1175/1520-0450(1977)016<0649:STPTAP>2.0.CO;2
  - *Validation of Priestley-Taylor coefficient across different climates*

- **Flint, A. L., & Childs, S. W. (1991).** Use of the Priestley-Taylor evaporation equation for soil water limited conditions in a small forest clearcut. *Agricultural and Forest Meteorology*, 56(3-4), 247-260. https://doi.org/10.1016/0168-1923(91)90094-7
  - *Application under water-limited conditions and coefficient calibration*

- **Fisher, J. B., Tu, K. P., & Baldocchi, D. D. (2008).** Global estimates of the land-atmosphere water flux based on monthly AVHRR and ISLSCP-II data, validated at 16 FLUXNET sites. *Remote Sensing of Environment*, 112(3), 901-919. https://doi.org/10.1016/j.rse.2007.06.025
  - *Global-scale validation using eddy covariance measurements*

### Supporting Methods
- **Verma, S. B., Rosenberg, N. J., & Blad, B. L. (1976).** Turbulent exchange coefficients for sensible heat and water vapor under advective conditions. *Journal of Applied Meteorology*, 15(4), 330-338. https://doi.org/10.1175/1520-0450(1976)015<0330:TECFSH>2.0.CO;2
  - *Net radiation calculation methodology*

- **Bastiaanssen, W. G. M., Menenti, M., Feddes, R. A., & Holtslag, A. A. M. (1998).** A remote sensing surface energy balance algorithm for land (SEBAL). 1. Formulation. *Journal of Hydrology*, 212-213, 198-212. https://doi.org/10.1016/S0022-1694(98)00253-4
  - *SEBAL algorithm for soil heat flux estimation*

### Theoretical Background
- **Jarvis, P. G., & McNaughton, K. G. (1986).** Stomatal control of transpiration: scaling up from leaf to region. *Advances in Ecological Research*, 15, 1-49. https://doi.org/10.1016/S0065-2504(08)60119-1
  - *Decoupling coefficients and canopy-atmosphere interactions*

- **Shuttleworth, W. J. (1993).** Evaporation. In *Handbook of Hydrology* (pp. 4.1-4.53). McGraw-Hill. ISBN: 0-07-039732-5
  - *Comprehensive review of evapotranspiration methods and applications*
