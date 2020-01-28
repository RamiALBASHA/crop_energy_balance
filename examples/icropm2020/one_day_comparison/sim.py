from math import radians
from json import load

from crop_irradiance.uniform_crops import inputs, params, shoot

with open('inputs.json', mode='r') as f:
    json_inputs = load(f)

hourly_direct_par = [0, 0, 0, 0, 0, 0, 0, 18.428820397314112, 86.779379638514726, 213.23527694182002,
                     332.53550250736947, 379.40012173450828, 359.77089500962774, 277.64770928541054,
                     133.08397549276251, 39.374833543598889, 1.1288678878776721, 0, 0, 0, 0, 0, 0, 0]
hourly_diffuse_par = [0, 0, 0, 0, 0, 0, 0, 18.428820397314112, 61.4976272200411, 74.918330146929719,
                      74.699625542626109,
                      85.077684125431659, 80.552408749094454, 65.365757925282779, 72.846107552867508,
                      37.432721347065836, 1.1288678878776721, 0, 0, 0, 0, 0, 0, 0]
hourly_solar_angle = [-22.15, - 22.91, -20.89, -16.32, -9.67, -0.67, 8.2, 18.34, 28.98, 39.81, 50.46, 60.27,
                      67.68, 69.63, 64.79, 56.01, 45.7, 34.91, 24.13, 13.67, 3.93, -5.28, -12.9, -18.68]

leaf_layers = {3: 1.0,
               2: 1.0,
               1: 1.0,
               0: 1.0}


def get_simulated_irradiance(leaves_category: str, hour: int):
    if leaves_category == 'sunlit-shaded':
        sim_inputs = inputs.SunlitShadedInputs(leaf_layers=leaf_layers,
                                               incident_direct_irradiance=hourly_direct_par[hour],
                                               incident_diffuse_irradiance=hourly_diffuse_par[hour],
                                               solar_inclination=radians(hourly_solar_angle[hour]))

        sim_params = params.SunlitShadedParams(leaf_reflectance=0.07, leaf_transmittance=0.00,
                                               leaves_to_sun_average_projection=0.5, sky_sectors_number=3,
                                               sky_type='soc',
                                               canopy_reflectance_to_diffuse_irradiance=0.057)
        sim_params.update(sim_inputs)
    else:
        sim_inputs = inputs.LumpedInputs(leaf_layers=leaf_layers,
                                         incident_irradiance=hourly_direct_par[hour] + hourly_diffuse_par[hour])
        sim_params = params.LumpedParams(extinction_coefficient=0.5)

    canopy = shoot.Shoot(leaves_category=leaves_category, inputs=sim_inputs, params=sim_params)

    absorbed_irradiance = {index: layer.absorbed_irradiance for index, layer in canopy.items()}
    soil_irradiance = (hourly_direct_par[hour] + hourly_diffuse_par[hour]) - (
        sum([sum(v.absorbed_irradiance.values()) for v in canopy.values()]))
    absorbed_irradiance.update({-1: {'lumped': soil_irradiance}})

    return absorbed_irradiance


def set_energy_balance_inputs(canopy_category: str,
                              leaves_category: str,
                              hour: int,
                              absorbed_irradiance: list) -> dict:
    simu_inputs = json_inputs.copy()
    if canopy_category == 'bigleaf':
        canopy_layers = {0: sum(leaf_layers.values())}
    else:
        canopy_layers = leaf_layers.copy()

    if leaves_category == 'lumped':
        incident_photosynthetically_active_radiation = {
            'lumped': hourly_direct_par[hour] + hourly_diffuse_par[hour]}
        solar_inclination = None
    else:
        incident_photosynthetically_active_radiation = {
            'direct': hourly_direct_par[hour],
            'diffuse': hourly_diffuse_par[hour]}
        solar_inclination = hourly_solar_angle[hour]

    simu_inputs.update(
        {"solar_inclination": solar_inclination,
         "leaf_layers": canopy_layers,
         "incident_photosynthetically_active_radiation": incident_photosynthetically_active_radiation,
         "absorbed_photosynthetically_active_radiation": absorbed_irradiance})

    return simu_inputs


abs_lumped_ls = []
abs_sunshade_ls = []
for hour in range(24):
    abs_lumped = get_simulated_irradiance(leaves_category='lumped', hour=hour)
    abs_sunshade = get_simulated_irradiance(leaves_category='sunlit-shaded', hour=hour)
    abs_lumped_ls.append(abs_lumped)
    abs_sunshade_ls.append(abs_sunshade)

    set_energy_balance_inputs(canopy_category='bigleaf',
                              leaves_category='lumped',
                              hour=0,
                              absorbed_irradiance=abs_lumped[0])
