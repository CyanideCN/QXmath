# cython: language_level=3
# distutils: sources = QXmathlib/QXmath.c
# distutils: include_dirs = QXmathlib/

from _qxmath cimport *

# Start of function exposed to Python
def lcl_temperature(double p, double t, double td):
    return Tc(p, t, td)

def saturated_vapor_pressure(double td):
    return E_WATER(td)

def saturated_vapor_pressure_ice(double td):
    return E_ICE(td)

def equivalent_potential_temperature(double p, double t, double td):
    return Qse(p, t, td)

def potential_temperature(double p, double t):
    return Qp(p, t)

def total_air_energy_surface(double p, double t, double td, double z):
    return Ttdm(p, t, td, z)

def total_air_energy(double p, double t, double td, double z, double v):
    return Ttgk(p, t, td, z, v)

def dewpoint_from_e(double E):
    return Etotd(E)

def specific_humidity(double p, double td):
    return qgk(p, td)

def relative_humidity(double p, double td):
    return rgk(p, td)

def condensation_function(double p, double td):
    return Fc(p, td)

def moist_lapse(double p, double t):
    return Rm(p, t)

def water_vapor_flux(double p, double td, double v):
    return sqtl(p, td, v)

def wind_to_uv(float fd, float ff):
    cdef float u = 0, v = 0
    fsfj(fd, ff, &u, &v)
    return u, v

def uv_to_wind(float u, float v):
    cdef float fd = 0, ff = 0
    fshc(u, v, &fd, &ff)
    return fd, ff

def showalter_index(double t850, double td850, double t500):
    return showalter(t850, td850, t500)

def richardson_number(double pdn, double tdn, double fddn, double ffdn, double pup, double tup, double fdup, double ffup):
    return richardson(pdn, tdn, fddn, ffdn, pup, tup, fdup, ffup)