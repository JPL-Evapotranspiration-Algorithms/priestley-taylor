from typing import Union
from datetime import datetime
import check_distribution
import numpy as np
from verma_net_radiation import verma_net_radiation
from SEBAL_soil_heat_flux import calculate_SEBAL_soil_heat_flux

from rasters import Raster, RasterGeometry
from GEOS5FP import GEOS5FP

RESAMPLING_METHOD = "cubic"

# Priestley-Taylor coefficient alpha (dimensionless)
# Typical value for unstressed vegetation
PT_ALPHA = 1.26

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

def calculate_epsilon(
        delta: Union[Raster, np.ndarray], 
        gamma: Union[Raster, np.ndarray]) -> Union[Raster, np.ndarray]:
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
    delta_Pa: Union[Raster, np.ndarray] = None,
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
    if delta_Pa is None:
        delta_Pa = delta_Pa_from_Ta_C(Ta_C)
    
    epsilon = calculate_epsilon(
        delta=delta_Pa, 
        gamma=gamma_Pa
    )

    return epsilon

def priestley_taylor(
    ST_C: Union[Raster, np.ndarray] = None,
    emissivity: Union[Raster, np.ndarray] = None,
    NDVI: Union[Raster, np.ndarray] = None,
    albedo: Union[Raster, np.ndarray] = None,
    Rn_Wm2: Union[Raster, np.ndarray] = None,
    G_Wm2: Union[Raster, np.ndarray] = None,
    SWin_Wm2: Union[Raster, np.ndarray] = None,
    Ta_C: Union[Raster, np.ndarray] = None,
    RH: Union[Raster, np.ndarray] = None,
    geometry: RasterGeometry = None,
    time_UTC: datetime = None,
    GEOS5FP_connection: GEOS5FP = None,
    resampling: str = RESAMPLING_METHOD,
    PT_alpha: float = PT_ALPHA,
    delta_Pa: Union[Raster, np.ndarray] = None,
    gamma_Pa: Union[Raster, np.ndarray, float] = GAMMA_PA
) -> Union[Raster, np.ndarray]:
    """
    Calculate the potential latent heat flux (LE) using the Priestley-Taylor equation.
    This method estimates the maximum possible evapotranspiration under optimal conditions,
    assuming sufficient water availability and minimal aerodynamic resistance.

    LE = α × ε × (Rn - G)

    where:
        LE = Latent heat flux (W/m²)
        α = Priestley-Taylor coefficient (dimensionless, typically 1.26 for unstressed vegetation)
        ε = Δ / (Δ + γ) (dimensionless)
        Rn = Net radiation at the surface (W/m²)
        G = Soil heat flux (W/m²)

    Args:
        Rn: Net radiation at the surface in W/m² (Raster or np.ndarray)
        G: Soil heat flux in W/m² (Raster or np.ndarray)
        Ta_C: Air temperature in degrees Celsius (Raster or np.ndarray)
        alpha: Priestley-Taylor coefficient (default: PT_ALPHA)
        gamma_Pa: Psychrometric constant in Pa/°C (default: GAMMA_PA)
    Returns:
        Potential latent heat flux LE in W/m²

    References:
        - Priestley, C. H. B., & Taylor, R. J. (1972). Monthly Weather Review, 100(2), 81–92.
        - Allen et al. (1998), FAO 56
    """
    results = {}

    # Create GEOS5FP connection if not provided
    if GEOS5FP_connection is None:
        GEOS5FP_connection = GEOS5FP()

    # Retrieve air temperature if not provided, using GEOS5FP and geometry/time
    if Ta_C is None and geometry is not None and time_UTC is not None:
        Ta_C = GEOS5FP_connection.Ta_C(
            time_UTC=time_UTC,
            geometry=geometry,
            resampling=resampling
        )

    if Ta_C is None:
        raise ValueError("air temperature (Ta_C) not given")

    # Compute net radiation if not provided, using albedo, ST_C, and emissivity
    if Rn_Wm2 is None and albedo is not None and ST_C is not None and emissivity is not None:
        # Retrieve incoming shortwave if not provided
        if SWin_Wm2 is None and geometry is not None and time_UTC is not None:
            SWin_Wm2 = GEOS5FP_connection.SWin(
                time_UTC=time_UTC,
                geometry=geometry,
                resampling=resampling
            )

        # Calculate net radiation using Verma et al. method
        Rn_results = verma_net_radiation(
            SWin_Wm2=SWin_Wm2,
            albedo=albedo,
            ST_C=ST_C,
            emissivity=emissivity,
            Ta_C=Ta_C,
            RH=RH,
            geometry=geometry,
            time_UTC=time_UTC,
            resampling=resampling,
            GEOS5FP_connection=GEOS5FP_connection
        )

        Rn_Wm2 = Rn_results["Rn_Wm2"]

    if Rn_Wm2 is None:
        raise ValueError("net radiation (Rn_Wm2) not given")

        # Compute soil heat flux if not provided, using SEBAL method
    if G_Wm2 is None and Rn_Wm2 is not None and ST_C is not None and NDVI is not None and albedo is not None:
        G_Wm2 = calculate_SEBAL_soil_heat_flux(
            Rn=Rn_Wm2,
            ST_C=ST_C,
            NDVI=NDVI,
            albedo=albedo
        )

    if G_Wm2 is None:
        raise ValueError("soil heat flux (G) not given")
    
    check_distribution(G_Wm2, "G_Wm2")
    results["G_Wm2"] = G_Wm2    

    epsilon = epsilon_from_Ta_C(
        Ta_C=Ta_C, 
        delta_Pa=delta_Pa,
        gamma_Pa=gamma_Pa
    )

    check_distribution(epsilon, "epsilon")
    results["epsilon"] = epsilon
    
    LE_potential_Wm2 = PT_alpha * epsilon * (Rn_Wm2 - G_Wm2)

    check_distribution(LE_potential_Wm2, "LE_potential_Wm2")
    results["LE_potential_Wm2"] = LE_potential_Wm2
    
    return results
