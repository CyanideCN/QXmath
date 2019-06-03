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

def showalter_index(double t850, double td850, double t500):
    return showalter(t850, td850, t500)