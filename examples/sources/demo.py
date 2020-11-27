import matplotlib.pyplot as plt
import pandas as pd

from crop_energy_balance import utils


def get_weather_data() -> pd.DataFrame:
    raw_data = pd.read_csv('weather.csv', decimal='.', sep=';', skiprows=6).set_index('time')
    raw_data.loc[:, 'wind_speed'] = raw_data.apply(lambda x: x['wind_speed'] * 3600.0, axis=1)
    raw_data.loc[:, 'incident_direct_irradiance'] = raw_data['incident_global_irradiance'].apply(
        lambda x: utils.convert_global_irradiance_into_photosynthetically_active_radiation(x * 0.80))
    raw_data.loc[:, 'incident_diffuse_irradiance'] = raw_data['incident_global_irradiance'].apply(
        lambda x: utils.convert_global_irradiance_into_photosynthetically_active_radiation(x * 0.20))
    raw_data.loc[:, 'vapor_pressure_deficit'] = raw_data.apply(
        lambda x: utils.calc_vapor_pressure_deficit(x['air_temperature'], x['air_temperature'], x['relative_humidity']),
        axis=1)
    raw_data.loc[:, 'vapor_pressure'] = raw_data.apply(
        lambda x: x['vapor_pressure_deficit'] * x['relative_humidity'] / 100., axis=1)

    raw_data.drop(['incident_global_irradiance', 'relative_humidity'], axis=1, inplace=True)

    return raw_data


def plot_weather():
    actual_weather = get_weather_data()

    fig, ((ax_irradiance, ax_temperature), (ax_wind_speed, ax_vpd)) = plt.subplots(ncols=2, nrows=2, sharex='all')

    day_hours = range(24)
    ax_irradiance.plot(
        day_hours, actual_weather.loc[:, ['incident_direct_irradiance', 'incident_diffuse_irradiance']].sum(axis=1))
    ax_irradiance.set_ylabel(r'$\mathregular{irradiance\/[W_{PAR} \cdot m^{-2}_{ground}]}$')

    ax_temperature.plot(day_hours, actual_weather.loc[:, 'air_temperature'])
    ax_temperature.set_ylabel(r'$\mathregular{temperature\/[^\circ C]}$')

    ax_wind_speed.plot(day_hours, actual_weather.loc[:, 'wind_speed'] / 3600.)
    ax_wind_speed.set_ylabel(r'$\mathregular{wind\/speed\/[m \cdot s^{-1}]}$')

    ax_vpd.plot(day_hours, actual_weather.loc[:, 'vapor_pressure_deficit'])
    ax_vpd.set_ylabel('VPD [kPa]')

    fig.tight_layout()
    fig.savefig('weather.png')
    plt.close()


if __name__ == '__main__':
    plot_weather()
