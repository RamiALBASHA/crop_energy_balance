from ephem import Sun, Observer
from datetime import datetime


latitude = 43.16
longitude = 3.87


def get_solar_inclination(date: datetime) -> float:
    sun = Sun()
    obs = Observer()
    obs.lat = latitude
    obs.lon = longitude
    obs.date = date
    sun.compute(obs)
    return sun.

hourly_direct_par = [0, 0, 0, 0, 0, 0, 0, 18.428820397314112, 86.779379638514726, 213.23527694182002,
                     332.53550250736947, 379.40012173450828, 359.77089500962774, 277.64770928541054,
                     133.08397549276251, 39.374833543598889, 1.1288678878776721, 0, 0, 0, 0, 0, 0, 0]
hourly_diffuse_par = [0, 0, 0, 0, 0, 0, 0, 18.428820397314112, 61.4976272200411, 74.918330146929719,
                      74.699625542626109,
                      85.077684125431659, 80.552408749094454, 65.365757925282779, 72.846107552867508,
                      37.432721347065836, 1.1288678878776721, 0, 0, 0, 0, 0, 0, 0]
hourly_solar_angle = [-15.86]

sun = Sun()
sun.compute(datetime(2012,7,20, 12))


