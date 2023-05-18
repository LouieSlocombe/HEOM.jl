# Taken from CODATA 2018 https://physics.nist.gov/cuu/Constants/Table/allascii.txt
# SI units
si_k_b = 1.380649e-23 # Boltzmann constant
si_h = 6.62607015e-34 # Planks constant
si_h_bar = 1.0545718176461565e-34 # Reduced planks constant
si_me = 9.1093837015e-31 # Mass of the electrons
si_mp = 1.67262192369e-27 # Mass of a proton
si_R = 8.314462618 # Gas constant
si_Na = 6.02214076e+23 # Avogadro's number
si_c = 299792458.0 # Speed of light in m/s
si_angstrom = 1e-10 # Angstrom
si_ns = 1e-9 # Nano seconds
si_ps = 1e-12 # Pico seconds
si_fs = 1e-15 # femto seconds
si_c_cm = si_c * 100.0 # Speed of light in cm/s


# Atomic units
au_energy = 4.3597447222071e-18
au_length = 5.29177210903e-11
au_mass = 9.1093837015e-31
au_time = 2.4188843265857e-17
au_me_mp = 1836.1526734400013
au_temp = 315775.02480406675
au_m_dalt = 1822.888486209
au_m_H = 1.007825031898 * au_m_dalt
au_m_D = 2.014101777844 * au_m_dalt

# Conversions
ev_k_b = 8.617333262e-05
ev_2_hart = 27.211386245988
hart_2_jmol = 2625499.6394798253
cm_2_j = 1.9864458571489285e-23
cm_2_hart = 4.556335252912009e-06
cm_2_aut_omega = cm_2_hart #* 2.0 * pi # doesnt seem to be required
hart_2_inv_cm = 1.0 / cm_2_hart #219474.63
cm_2_hz = si_c_cm # multiply it by the speed of light
bohr_2_ang = au_length / si_angstrom

kcalmol_2_kjmol = 4.184
kcalmol_2_hart = 0.001593601437640628
aut_2_ns = au_time / si_ns
aut_2_ps = au_time / si_ps
aut_2_fs = au_time / si_fs


# More constants
h_bar = 1.0 # In hartrees
h = h_bar * 2.0 * pi
k_b = si_k_b / au_energy # In hartrees

# Other constants
# https://en.wikipedia.org/wiki/Euler%27s_constant
C_e = 0.57721566490153286060651209008240243104215933593992

# Variables
mass = au_me_mp
temperature = 298.15
beta = 1.0 / (temperature * k_b)
t_relax = 0.1 * 1.0e-12 / au_time
gamma = 1.0 / t_relax

# Precalculate terms
term_kin = (1.0im * h_bar) / (2.0 * mass)
term_pot = -1.0im / h_bar
term_dissi = -gamma
term_deco = (-2.0 * gamma * mass) / (beta * h_bar^2)
