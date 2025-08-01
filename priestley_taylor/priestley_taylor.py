from typing import Union
import numpy as np

from rasters import Raster

GAMMA_KPA = 0.0662  # kPa/C
"""
Psychrometric constant gamma in kiloPascal per degree Celsius (kPa/°C).
This value is for ventilated (Asmann type) psychrometers with an air movement of ~5 m/s.
It is a key parameter in physically-based evapotranspiration models, linking the energy and aerodynamic terms.
Reference: Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration – Guidelines for computing crop water requirements – FAO Irrigation and drainage paper 56. FAO, Rome. Table 2.2.
"""

GAMMA_PA = GAMMA_KPA * 1000
"""
Psychrometric constant gamma in Pascal per degree Celsius (Pa/°C).
This is a direct unit conversion from GAMMA_KPA (1 kPa = 1000 Pa).
Reference: Allen et al. (1998), FAO 56.
"""

def delta_kPa_from_Ta_C(Ta_C: Union[Raster, np.ndarray]) -> Union[Raster, np.ndarray]:
    """
    Calculate the slope of the saturation vapor pressure curve (Δ, delta) at a given air temperature (°C),
    returning the result in kPa/°C. This is a key parameter in the Penman-Monteith and Priestley-Taylor equations,
    quantifying the sensitivity of saturation vapor pressure to temperature changes.

    Δ = 4098 × [0.6108 × exp(17.27 × Ta / (237.7 + Ta))] / (Ta + 237.3)²

    Args:
        Ta_C: Air temperature in degrees Celsius (Raster or np.ndarray)
    Returns:
        Slope of saturation vapor pressure curve (kPa/°C)

    References:
        - Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). FAO Irrigation and Drainage Paper 56, Eq. 2.18.
        - Monteith, J. L. (1965). Evaporation and environment. In The State and Movement of Water in Living Organisms (pp. 205–234). Academic Press.
    """
    return 4098 * (0.6108 * np.exp(17.27 * Ta_C / (237.7 + Ta_C))) / (Ta_C + 237.3) ** 2

def delta_Pa_from_Ta_C(Ta_C: Union[Raster, np.ndarray]) -> Union[Raster, np.ndarray]:
    """
    Convert the slope of the saturation vapor pressure curve (Δ) from kPa/°C to Pa/°C.
    This is a unit conversion used in some formulations of evapotranspiration models.

    Args:
        Ta_C: Air temperature in degrees Celsius (Raster or np.ndarray)
    Returns:
        Slope of saturation vapor pressure curve (Pa/°C)

    Reference:
        - Allen et al. (1998), FAO 56.
    """
    return delta_kPa_from_Ta_C(Ta_C) * 1000

def calculate_epsilon(delta: Union[Raster, np.ndarray], gamma: Union[Raster, np.ndarray]) -> Union[Raster, np.ndarray]:
    """
    Compute the dimensionless ratio epsilon (ε), defined as ε = Δ / (Δ + γ),
    where Δ is the slope of the saturation vapor pressure curve and γ is the psychrometric constant.
    This ratio is fundamental in the Priestley-Taylor and Penman-Monteith equations, representing
    the relative importance of energy supply versus atmospheric demand in controlling evapotranspiration.

    Args:
        delta: Slope of saturation vapor pressure curve (Pa/°C or kPa/°C)
        gamma: Psychrometric constant (same units as delta)
    Returns:
        Epsilon (dimensionless)

    References:
        - Priestley, C. H. B., & Taylor, R. J. (1972). On the assessment of surface heat flux and evaporation using large-scale parameters. Monthly Weather Review, 100(2), 81–92.
        - Allen et al. (1998), FAO 56, Eq. 6.2
    """
    return delta / (delta + gamma)

def epsilon_from_Ta_C(
    Ta_C: Union[Raster, np.ndarray],
    gamma_Pa: Union[Raster, np.ndarray, float] = GAMMA_PA
) -> Union[Raster, np.ndarray]:
    """
    Calculate epsilon (ε) from air temperature (°C) and the psychrometric constant (Pa/°C).
    This function computes Δ in Pa/°C from temperature, then calculates ε = Δ / (Δ + γ).
    Epsilon is a key parameter in the Priestley-Taylor equation for potential evapotranspiration,
    determining the partitioning of available energy into latent heat flux.

    Args:
        Ta_C: Air temperature in degrees Celsius (Raster or np.ndarray)
        gamma_Pa: Psychrometric constant in Pa/°C (default: GAMMA_PA)
    Returns:
        Epsilon (dimensionless)

    References:
        - Priestley, C. H. B., & Taylor, R. J. (1972). Monthly Weather Review, 100(2), 81–92.
        - Allen et al. (1998), FAO 56
    """
    delta_Pa = delta_Pa_from_Ta_C(Ta_C)
    epsilon = calculate_epsilon(delta=delta_Pa, gamma=gamma_Pa)
    return epsilon
