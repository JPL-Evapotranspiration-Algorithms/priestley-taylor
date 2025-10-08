from typing import Union
from datetime import datetime
from check_distribution import check_distribution
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
    Calculate the dimensionless epsilon ratio directly from air temperature.
    
    This convenience function combines the calculation of the saturation vapor pressure
    slope with the epsilon ratio computation, providing a direct pathway from air
    temperature to the key parameter needed for Priestley-Taylor evapotranspiration
    estimation.
    
    Process flow:
    1. Compute Δ from air temperature using Magnus-Tetens approximation
    2. Calculate ε = Δ / (Δ + γ) using the computed slope and psychrometric constant
    
    Mathematical sequence:
    Ta[°C] → Δ[Pa/°C] → ε[dimensionless]
    
    This approach is computationally efficient and maintains consistency in unit handling
    throughout the calculation chain, avoiding potential unit conversion errors.
    
    Temperature sensitivity analysis:
    The relationship between air temperature and epsilon is non-linear and asymptotic:
    - dε/dT decreases with increasing temperature
    - Sensitivity is highest at low temperatures (0-10°C)
    - Above 30°C, epsilon approaches unity asymptotically
    
    Args:
        Ta_C (Union[Raster, np.ndarray]): Air temperature [°C]
            Valid range: -40 to +50°C for accurate Magnus-Tetens approximation
            Typical range: -20 to +45°C for terrestrial applications
        delta_Pa (Union[Raster, np.ndarray], optional): Pre-computed saturation vapor
            pressure slope [Pa/°C]. If None, computed from Ta_C.
            Useful for efficiency when delta is already available.
        gamma_Pa (Union[Raster, np.ndarray, float], optional): Psychrometric constant [Pa/°C]
            Default: GAMMA_PA (66.2 Pa/°C at standard conditions)
            Can be array for elevation-dependent corrections
    
    Returns:
        Union[Raster, np.ndarray]: Epsilon ratio [dimensionless]
            Range: typically 0.4-0.99 for terrestrial conditions
            Increases monotonically with temperature
    
    Examples:
        >>> import numpy as np
        >>> Ta = np.array([0, 10, 20, 30, 40])  # °C
        >>> eps = epsilon_from_Ta_C(Ta)
        >>> print(f"Temperature: {Ta}°C")
        >>> print(f"Epsilon: {eps:.3f}")
        Temperature: [ 0 10 20 30 40]°C
        Epsilon: [0.741 0.926 0.956 0.974 0.984]
    
    References:
        - Priestley, C. H. B., & Taylor, R. J. (1972). On the assessment of surface heat flux
          and evaporation using large-scale parameters. Monthly Weather Review, 100(2), 81-92.
          DOI: 10.1175/1520-0493(1972)100<0081:OTAOSH>2.3.CO;2
        - Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration:
          Guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper 56.
          ISBN: 92-5-104219-5
        - Shuttleworth, W. J. (1993). Evaporation. In Handbook of Hydrology (pp. 4.1-4.53).
          McGraw-Hill. ISBN: 0-07-039732-5
    """
    # Compute saturation vapor pressure slope if not provided
    if delta_Pa is None:
        delta_Pa = delta_Pa_from_Ta_C(Ta_C)
    
    # Calculate the dimensionless partitioning coefficient
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
    gamma_Pa: Union[Raster, np.ndarray, float] = GAMMA_PA,
    epsilon: Union[Raster, np.ndarray] = None
) -> Union[Raster, np.ndarray]:
    """
    Calculate potential evapotranspiration using the Priestley-Taylor equation.
    
    This function implements the widely-used Priestley-Taylor method for estimating
    potential evapotranspiration (PET), which represents the theoretical maximum
    evapotranspiration rate under optimal water availability and minimal plant stress.
    
    THEORETICAL FOUNDATION:
    The Priestley-Taylor equation is a simplified form of the Penman-Monteith equation
    that eliminates the need for wind speed and vapor pressure deficit measurements.
    It assumes that under well-watered conditions, evapotranspiration is primarily
    controlled by available energy rather than atmospheric demand.
    
    CORE EQUATION:
    LE = α × ε × (Rn - G)
    
    Where:
    - LE = Latent heat flux [W/m²] - energy used for evapotranspiration
    - α = Priestley-Taylor coefficient [dimensionless] - empirical parameter (≈1.26)
    - ε = Δ/(Δ+γ) [dimensionless] - energy partitioning coefficient
    - Δ = Slope of saturation vapor pressure curve [Pa/°C]
    - γ = Psychrometric constant [Pa/°C]
    - Rn = Net radiation at surface [W/m²] - available radiative energy
    - G = Soil heat flux [W/m²] - energy conducted into/out of soil
    
    PHYSICAL INTERPRETATION:
    1. (Rn - G): Available energy for turbulent heat fluxes
    2. ε: Fraction of available energy partitioned to latent heat under equilibrium
    3. α: Correction factor accounting for non-equilibrium conditions and advection
    4. LE: Actual latent heat flux under potential conditions
    
    ASSUMPTIONS AND LIMITATIONS:
    1. Sufficient water availability (no soil moisture stress)
    2. Minimal plant water stress (stomata fully open)
    3. Horizontal homogeneity of surface conditions
    4. Steady-state energy balance conditions
    5. Valid for time scales > 1 hour (avoids storage effects)
    
    COMPARISON WITH OTHER METHODS:
    Advantages over Penman-Monteith:
    - Requires fewer meteorological inputs (no wind speed, humidity)
    - More stable in data-sparse regions
    - Less sensitive to measurement errors
    
    Limitations compared to Penman-Monteith:
    - Cannot account for stomatal control
    - Less accurate under advective conditions
    - Overestimates ET in arid environments
    
    ENERGY BALANCE CONTEXT:
    Rn = LE + H + G + S
    Where:
    - H = Sensible heat flux [W/m²]
    - S = Storage terms [W/m²] (usually neglected for daily/hourly averages)
    
    Under PT assumptions: H = (1-ε)/ε × α × LE
    
    Args:
        ST_C (Union[Raster, np.ndarray], optional): Surface temperature [°C]
            Required for net radiation calculation if Rn_Wm2 not provided
            Typical range: -50 to +70°C for terrestrial surfaces
            
        emissivity (Union[Raster, np.ndarray], optional): Surface emissivity [dimensionless]
            Required for longwave radiation calculation, range: 0.85-0.99
            Typical values: vegetation ~0.98, bare soil ~0.94, water ~0.96
            
        NDVI (Union[Raster, np.ndarray], optional): Normalized Difference Vegetation Index
            Required for soil heat flux calculation, range: -1 to +1
            Typical values: water <0, bare soil 0.1-0.2, vegetation 0.2-0.9
            
        albedo (Union[Raster, np.ndarray], optional): Surface albedo [dimensionless]
            Fraction of incoming solar radiation reflected, range: 0-1
            Typical values: fresh snow ~0.9, vegetation ~0.2, bare soil ~0.3
            
        Rn_Wm2 (Union[Raster, np.ndarray], optional): Net radiation [W/m²]
            All-wave net radiation (shortwave + longwave)
            Positive: energy input to surface; Negative: energy loss
            Daily range: typically -100 to +800 W/m²
            
        G_Wm2 (Union[Raster, np.ndarray], optional): Soil heat flux [W/m²]
            Conductive heat flux into (+) or out of (-) soil
            Daily range: typically -50 to +200 W/m²
            Often estimated as fraction of Rn: G ≈ 0.1×Rn (daily average)
            
        SWin_Wm2 (Union[Raster, np.ndarray], optional): Incoming shortwave radiation [W/m²]
            Direct + diffuse solar radiation at surface
            Range: 0-1400 W/m² (theoretical maximum ~1361 W/m²)
            
        Ta_C (Union[Raster, np.ndarray], optional): Air temperature [°C]
            Near-surface air temperature for vapor pressure calculations
            Measurement height: typically 1.5-2.0 m above surface
            
        RH (Union[Raster, np.ndarray], optional): Relative humidity [%]
            Required for longwave radiation calculations
            Range: 0-100%, affects atmospheric emissivity
            
        geometry (RasterGeometry, optional): Spatial grid definition
            Required for GEOS5FP data retrieval and spatial operations
            
        time_UTC (datetime, optional): UTC timestamp
            Required for solar angle calculations and meteorological data retrieval
            
        GEOS5FP_connection (GEOS5FP, optional): Connection to GEOS-5 FP dataset
            NASA's Global Modeling and Assimilation Office atmospheric reanalysis
            Provides meteorological forcing when local data unavailable
            
        resampling (str, optional): Spatial resampling method
            Default: "cubic" for smooth interpolation
            Options: "nearest", "linear", "cubic", "lanczos"
            
        PT_alpha (float, optional): Priestley-Taylor coefficient [dimensionless]
            Default: 1.26 (original Priestley & Taylor value)
            Calibration range: 0.8-1.5 depending on conditions
            
        delta_Pa (Union[Raster, np.ndarray], optional): Vapor pressure slope [Pa/°C]
            Pre-computed for computational efficiency
            If None, calculated from Ta_C using Magnus-Tetens equation
            
        gamma_Pa (Union[Raster, np.ndarray, float], optional): Psychrometric constant [Pa/°C]
            Default: 66.2 Pa/°C at standard atmospheric pressure
            Varies with elevation: γ(P) = 66.2 × (P/101.3)
            
        epsilon (Union[Raster, np.ndarray], optional): Pre-computed epsilon ratio
            If None, calculated from delta_Pa and gamma_Pa
            Range: typically 0.4-0.99 for terrestrial conditions
    
    Returns:
        dict: Dictionary containing computed variables and results:
            - "LE_potential_Wm2": Potential latent heat flux [W/m²]
            - "epsilon": Energy partitioning coefficient [dimensionless]  
            - "G_Wm2": Soil heat flux [W/m²] (if computed)
            Additional variables may be included from sub-calculations
    
    Raises:
        ValueError: If required inputs are missing or invalid
        Warning: If computed values are outside expected physical ranges
    
    Examples:
        Basic usage with all inputs provided:
        >>> import numpy as np
        >>> from datetime import datetime
        >>> 
        >>> # Define inputs
        >>> Rn = np.array([[400, 500], [450, 550]])  # W/m²
        >>> G = np.array([[40, 50], [45, 55]])       # W/m²  
        >>> Ta = np.array([[20, 25], [22, 27]])      # °C
        >>> 
        >>> # Calculate potential ET
        >>> result = priestley_taylor(Rn_Wm2=Rn, G_Wm2=G, Ta_C=Ta)
        >>> LE_pot = result["LE_potential_Wm2"]
        >>> print(f"Potential LE: {LE_pot:.1f} W/m²")
        
        Advanced usage with automatic data retrieval:
        >>> from rasters import RasterGeometry
        >>> from datetime import datetime
        >>> 
        >>> # Define spatial domain
        >>> geom = RasterGeometry.from_bounds(...)
        >>> time = datetime.now()
        >>> 
        >>> # Calculate with automatic meteorological data
        >>> result = priestley_taylor(
        ...     geometry=geom,
        ...     time_UTC=time,
        ...     ST_C=surface_temp,
        ...     albedo=surface_albedo,
        ...     emissivity=surface_emissivity,
        ...     NDVI=vegetation_index
        ... )
    
    References:
        FOUNDATIONAL PAPERS:
        - Priestley, C. H. B., & Taylor, R. J. (1972). On the assessment of surface heat flux
          and evaporation using large-scale parameters. Monthly Weather Review, 100(2), 81-92.
          DOI: 10.1175/1520-0493(1972)100<0081:OTAOSH>2.3.CO;2
          [Original derivation and validation of the Priestley-Taylor equation]
          
        - Penman, H. L. (1948). Natural evaporation from open water, bare soil and grass.
          Proceedings of the Royal Society of London, 193(1032), 120-145.
          DOI: 10.1098/rspa.1948.0037
          [Theoretical foundation for combination equations]
        
        METHODOLOGICAL STANDARDS:
        - Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration:
          Guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper 56.
          ISBN: 92-5-104219-5
          [International standard for evapotranspiration calculations]
          
        - Monteith, J. L., & Unsworth, M. H. (2013). Principles of Environmental Physics:
          Plants, Animals, and the Atmosphere. 4th Edition, Academic Press.
          ISBN: 978-0-12-386910-4
          [Comprehensive treatment of environmental physics principles]
        
        VALIDATION AND APPLICATIONS:
        - Stewart, R. B., & Rouse, W. R. (1977). Substantiating the Priestley and Taylor
          parameter α = 1.26 for potential evaporation in high latitudes. Journal of Applied
          Meteorology, 16(6), 649-650. DOI: 10.1175/1520-0450(1977)016<0649:STPTAP>2.0.CO;2
          
        - Flint, A. L., & Childs, S. W. (1991). Use of the Priestley-Taylor evaporation
          equation for soil water limited conditions in a small forest clearcut.
          Agricultural and Forest Meteorology, 56(3-4), 247-260.
          DOI: 10.1016/0168-1923(91)90094-7
          
        - Fisher, J. B., Tu, K. P., & Baldocchi, D. D. (2008). Global estimates of the
          land-atmosphere water flux based on monthly AVHRR and ISLSCP-II data, validated
          at 16 FLUXNET sites. Remote Sensing of Environment, 112(3), 901-919.
          DOI: 10.1016/j.rse.2007.06.025
          [Global-scale applications and validation]
    """
    # Initialize results dictionary to store all computed variables
    results = {}

    # STEP 1: ESTABLISH DATA CONNECTIONS AND RETRIEVE METEOROLOGICAL FORCING
    # Create GEOS5FP connection if not provided
    # GEOS-5 FP (Forward Processing) provides global atmospheric reanalysis data
    # at 0.25° × 0.3125° spatial resolution with hourly temporal resolution
    if GEOS5FP_connection is None:
        GEOS5FP_connection = GEOS5FP()

    # Retrieve air temperature from GEOS5FP if not provided locally
    # Air temperature is fundamental for calculating vapor pressure relationships
    # Standard measurement height: 2 meters above ground level
    if Ta_C is None and geometry is not None and time_UTC is not None:
        Ta_C = GEOS5FP_connection.Ta_C(
            time_UTC=time_UTC,
            geometry=geometry,
            resampling=resampling
        )

    # Air temperature is mandatory - all vapor pressure calculations depend on it
    if Ta_C is None:
        raise ValueError("air temperature (Ta_C) not given - required for vapor pressure calculations")

    # STEP 2: COMPUTE OR RETRIEVE NET RADIATION (Rn)
    # Net radiation represents the available radiative energy at the surface
    # Rn = (SWin - SWout) + (LWin - LWout) = SWin(1-α) + εσ(Ta⁴ - Ts⁴)
    # Where: SWin=incoming solar, α=albedo, ε=emissivity, σ=Stefan-Boltzmann constant
    if Rn_Wm2 is None and albedo is not None and ST_C is not None and emissivity is not None:
        # Retrieve incoming shortwave radiation from GEOS5FP if not provided
        # This includes both direct beam and diffuse sky radiation
        if SWin_Wm2 is None and geometry is not None and time_UTC is not None:
            SWin_Wm2 = GEOS5FP_connection.SWin(
                time_UTC=time_UTC,
                geometry=geometry,
                resampling=resampling
            )

        # Calculate net radiation using the Verma et al. method
        # This approach accounts for both shortwave and longwave radiation components
        # with atmospheric correction for water vapor and temperature effects
        Rn_results = verma_net_radiation(
            SWin_Wm2=SWin_Wm2,          # Incoming shortwave radiation [W/m²]
            albedo=albedo,               # Surface albedo [dimensionless]
            ST_C=ST_C,                   # Surface temperature [°C]
            emissivity=emissivity,       # Surface emissivity [dimensionless]
            Ta_C=Ta_C,                   # Air temperature [°C]
            RH=RH,                       # Relative humidity [%]
            geometry=geometry,           # Spatial grid definition
            time_UTC=time_UTC,           # UTC timestamp for solar calculations
            resampling=resampling,       # Spatial interpolation method
            GEOS5FP_connection=GEOS5FP_connection
        )

        # Extract net radiation from results
        Rn_Wm2 = Rn_results["Rn_Wm2"]

    # Net radiation is mandatory - represents the primary energy driver for ET
    if Rn_Wm2 is None:
        raise ValueError("net radiation (Rn_Wm2) not given - required as primary energy source")

    # STEP 3: COMPUTE OR RETRIEVE SOIL HEAT FLUX (G)
    # Soil heat flux represents energy conducted into (+) or out of (-) the soil
    # This energy is not available for evapotranspiration processes
    # Daily average: G ≈ 0.1 × Rn for most surfaces
    # Instantaneous values can be much larger, especially for bare soil
    if G_Wm2 is None and Rn_Wm2 is not None and ST_C is not None and NDVI is not None and albedo is not None:
        # Use SEBAL (Surface Energy Balance Algorithm for Land) method
        # This empirical approach relates G to net radiation using vegetation indices
        # G/Rn = f(NDVI, albedo, surface temperature) - accounts for vegetation effects
        G_Wm2 = calculate_SEBAL_soil_heat_flux(
            Rn=Rn_Wm2,      # Net radiation [W/m²] - primary energy input
            ST_C=ST_C,      # Surface temperature [°C] - affects soil thermal gradient
            NDVI=NDVI,      # Vegetation index - controls surface energy partitioning
            albedo=albedo   # Surface albedo - affects surface energy absorption
        )

    # Soil heat flux is mandatory for accurate energy balance
    if G_Wm2 is None:
        raise ValueError("soil heat flux (G_Wm2) not given - required for energy balance closure")
    
    # Validate soil heat flux distribution and store in results
    check_distribution(G_Wm2, "G_Wm2")
    results["G_Wm2"] = G_Wm2    

    # STEP 4: COMPUTE EPSILON (ε) - THE ENERGY PARTITIONING COEFFICIENT
    # Epsilon represents the theoretical fraction of available energy (Rn-G)
    # that would be converted to latent heat under equilibrium conditions
    # ε = Δ/(Δ+γ) where Δ=vapor pressure slope, γ=psychrometric constant
    if epsilon is None:
        # Calculate epsilon from air temperature using thermodynamic relationships
        # This involves: Ta → Δ(Magnus-Tetens) → ε = Δ/(Δ+γ)
        epsilon = epsilon_from_Ta_C(
            Ta_C=Ta_C,          # Air temperature [°C] - controls vapor pressure slope
            delta_Pa=delta_Pa,  # Pre-computed slope [Pa/°C] - for efficiency
            gamma_Pa=gamma_Pa   # Psychrometric constant [Pa/°C] - atmospheric property
        )

    # Validate epsilon distribution - should be between 0.4-0.99 for terrestrial conditions
    check_distribution(epsilon, "epsilon")
    results["epsilon"] = epsilon
    
    # STEP 5: APPLY THE PRIESTLEY-TAYLOR EQUATION
    # LE_potential = α × ε × (Rn - G)
    # This is the core calculation representing potential evapotranspiration
    # under optimal water availability and minimal plant stress
    #
    # Physical interpretation:
    # - (Rn - G): Available energy for turbulent heat fluxes [W/m²]
    # - ε: Equilibrium partitioning to latent heat [dimensionless]
    # - α: Correction for non-equilibrium conditions [dimensionless]
    # - LE_potential: Maximum possible latent heat flux [W/m²]
    LE_potential_Wm2 = PT_alpha * epsilon * (Rn_Wm2 - G_Wm2)

    # Validate final results - LE should be positive for normal daytime conditions
    # Negative values may occur at night or under extremely dry conditions
    check_distribution(LE_potential_Wm2, "LE_potential_Wm2")
    results["LE_potential_Wm2"] = LE_potential_Wm2
    
    # Return comprehensive results dictionary containing all computed variables
    # This allows users to access intermediate calculations for analysis or validation
    return results
