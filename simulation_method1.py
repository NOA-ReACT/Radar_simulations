import numpy as np
import matplotlib.pyplot as plt
from pytmatrix.tmatrix import Scatterer
from pytmatrix import radar, orientation
from simulations.functions import parametrization_axis_ratio, velocity_diameters_rel, fast_plotting
import os
from scipy import interpolate
import random

dir_out = r'simulations/plots'

config_file = 'simulations/test_config_file.txt'

# In[Read the inputs from the configuration file]
config_vars = {} # This is a dictionary to store config variables

with open(config_file, 'r') as f:
    for line in f:
        line = line.strip()

        if not line or line.startswith('#'): # Ignore empty lines or lines starting with '#'
            continue

        key, value = line.split('=', 1) # Split the line into key and value at the '=' symbol
        key = key.strip()  # Remove whitespace around the key
        value = value.split('#', 1)[0].strip()  # Remove comments after the value

        try: # convert the value to an int or float
            config_vars[key] = float(value) if '.' in value else int(value)
        except ValueError:
            config_vars[key] = value  # Keep as a string if conversion fails

# Read the T-MATRIX inputs
frequency = config_vars['frequency']
wl = 299792458/ (frequency*1e6)
ri = complex(config_vars['ri_Re_part'], config_vars['ri_Im_part'])
k_square = config_vars['k_square']
phi_0, phi = config_vars['phi_0'], config_vars['phi']
d_start, d_end = config_vars['d_start'], config_vars['d_end']
theta_0, theta = config_vars['theta_0'], config_vars['theta']
std = config_vars['std']

# Read the RADAR inputs
M, int_time = config_vars['M'], config_vars['int_time'] # number of bins, integrasion time in s
PRF_kHz = config_vars['PRF_kHz'] # in Herz
SNR_V_dB = config_vars['SNR_V_dB'] # in linear
K = 2*int(int_time*PRF_kHz*1e3/M)
SNR_V_lin = 10**(SNR_V_dB/10)
v_nyquist = PRF_kHz*wl/4/np.sin(np.pi*45/180) ## PF in kHz, wl in mm --> v_nyquist in m/s

# Read the GAMMA DISTRIBUTION inputs
mu, D_0, N_0 = config_vars['mu'], config_vars['D_0'], config_vars['N_0']
lamda = (3.67+mu)/D_0
sigma_v = config_vars['sigma_v'] # in m/s
sigma_f = 2*sigma_v/wl


# In[Calculate the spectrum parameters]

diameters = np.arange(d_start, d_end, (d_end-d_start)/1000, dtype=float)
axis_ratio = parametrization_axis_ratio(diameters)
velocities = np.array([velocity_diameters_rel(D) for D in diameters/1000])

N = np.array([N_0 * np.exp(-lamda * D) * D ** (mu) for D in diameters])  # raindrops size distribution

sigma_bsc_vert, sigma_bsc_horiz, rho_hv, delta_hv, zdr  = ([] for l in range(5))

for d, ax in zip(diameters, axis_ratio):
    scatterer = Scatterer(radius=d / 2, wavelength=wl, m=ri, axis_ratio=ax, thet0=theta_0, thet=theta,
                          phi0=phi_0, phi=phi)

    # scatterer.orient = orientation.orient_averaged_fixed ###uncomment only of canting of drops is desirable
    # scatterer.or_pdf = orientation.gaussian_pdf(std=std)

    rho_hv.append(radar.rho_hv(scatterer))
    delta_hv.append(radar.delta_hv(scatterer))
    zdr.append(radar.Zdr(scatterer))
    sigma_bsc_vert.append(radar.radar_xsect(scatterer, h_pol=False))

    scatterer_h = Scatterer(radius=d / 2, wavelength=wl, m=ri, axis_ratio=ax, thet0=theta_0, thet=theta,
              phi0=phi_0, phi=phi, h_pol=True)

    sigma_bsc_horiz.append(radar.radar_xsect(scatterer_h, h_pol=True))

sigma_bsc_horiz, sigma_bsc_vert = np.array(sigma_bsc_horiz), np.array(sigma_bsc_vert)
rho_hv, delta_hv, zdr = np.array(rho_hv), np.array(delta_hv), np.array(zdr)


additional_points = 50   # points to extend the spectrum around the spectral characteristics
extended_vel = np.concatenate((np.linspace(-v_nyquist, 0, additional_points), velocities, np.linspace(velocities.max(), v_nyquist, additional_points)))
extended_sigma_v = np.concatenate((np.zeros(additional_points), sigma_bsc_vert, np.zeros(additional_points)))
extended_sigma_h = np.concatenate((np.zeros(additional_points), sigma_bsc_horiz, np.zeros(additional_points)))
extended_dDdV = np.concatenate((np.zeros(additional_points), np.gradient(diameters,velocities), np.zeros(additional_points)))
extended_N = np.concatenate((np.zeros(additional_points), N, np.zeros(additional_points)))
extended_rho = np.concatenate((np.ones(additional_points), rho_hv, np.zeros(additional_points)))
extended_delta = np.concatenate((np.zeros(additional_points), delta_hv, np.zeros(additional_points)))
extended_zdr = np.concatenate((np.ones(additional_points), zdr, np.ones(additional_points)))

extended_Spectra_vert = wl**4/(np.pi**5*k_square) * extended_N * extended_sigma_v * extended_dDdV * wl/2  # spectral power in freq domain mm6m-3*s (frequency domain)
extended_Spectra_horiz = wl**4/(np.pi**5*k_square) * extended_N * extended_sigma_h * extended_dDdV * wl/2

frequencies = 2 * extended_vel/wl
frequencies_interpolated = np.linspace(-2*v_nyquist/wl, 2*v_nyquist/wl, M, dtype=float) #kHz
velocities_interpolated = frequencies_interpolated*wl/2

# Interpolations of spectra and polarimetric variables
f = interpolate.interp1d(frequencies, extended_Spectra_vert, fill_value="extrapolate")
interpolated_spectrum_V = f(frequencies_interpolated)
f = interpolate.interp1d(frequencies, extended_Spectra_horiz, fill_value="extrapolate")
interpolated_spectrum_H = f(frequencies_interpolated)
f = interpolate.interp1d(frequencies, extended_rho, fill_value="extrapolate")
interpolated_rho = f(frequencies_interpolated)
f = interpolate.interp1d(frequencies, extended_zdr, fill_value="extrapolate")
interpolated_zdr = f(frequencies_interpolated)
f = interpolate.interp1d(frequencies, extended_delta, fill_value="extrapolate")
interpolated_delta = f(frequencies_interpolated)

extended_Spectra_hv = wl**4/(np.pi**5*k_square) * extended_N * np.sqrt(extended_sigma_v*extended_sigma_h) * extended_dDdV * wl/2 * extended_rho * np.exp(1j*extended_delta)
f = interpolate.interp1d(frequencies, extended_Spectra_hv, fill_value="extrapolate") # (x,y)
interpolated_spectrum_HV = f(frequencies_interpolated)

# In[Turbulence Implementation]

v_turb=4 # m/s The turbulent velocities
frequencies_turbulence = np.arange(-v_turb*2/wl, v_turb*2/wl, step=np.diff(frequencies_interpolated)[0], dtype=float)

S_air = np.array(1/(np.sqrt(2*np.pi)*sigma_f)*np.exp(-frequencies_turbulence**2 / (2*sigma_f**2)))

Svv_turb = np.convolve(interpolated_spectrum_V, S_air, mode='same')/S_air.sum(axis=0)
Shh_turb = np.convolve(interpolated_spectrum_H, S_air, mode='same')/S_air.sum(axis=0)
Shv_turb = np.convolve(interpolated_spectrum_HV, S_air, mode='same')/S_air.sum(axis=0)

cross_variable = Shv_turb/np.sqrt(Svv_turb*Shh_turb)

rho_hv_turb = np.abs(cross_variable)
delta_hv_turb = np.angle(cross_variable)
zdr_turb = Shh_turb/Svv_turb

# In[Noise Implementation]

V_v, V_h, V_v_noise, V_h_noise = ([] for l in range(4))

noise_level_linear = (np.trapz(Svv_turb, x=frequencies_interpolated))/(2*frequencies.max()*SNR_V_lin) # spectral noise in the frequency domain
mean_value = noise_level_linear # spectral noise

for k in range(K):
    randoms_u = np.random.uniform(0,1,size=M)
    randoms_theta = -np.pi + 2*np.pi*np.random.uniform(0,1,size=M)

    #### Fourier Transform of V-channel
    V_fv = np.sqrt(-Svv_turb * np.log(randoms_u)) * np.exp(1j * randoms_theta)
    V_v.append(V_fv)

    #### Vertical
    randoms = []
    [randoms.append(random.random()) for i in range(M)]
    randoms_log = -np.log(randoms)
    white_noise_v = mean_value * randoms_log

    randoms_theta = -np.pi + 2 * np.pi * np.random.uniform(0, 1, size=M)
    randoms_u = np.random.uniform(0, 1, size=M)
    V_v_noise.append(np.sqrt(-white_noise_v * np.log(randoms_u)) * np.exp(1j * randoms_theta))

    #### Horizontal
    randoms = []
    [randoms.append(random.random()) for i in range(M)]
    randoms_log = -np.log(randoms)
    white_noise_h = mean_value * randoms_log
    randoms_theta = -np.pi + 2 * np.pi * np.random.uniform(0, 1, size=M)
    randoms_u = np.random.uniform(0, 1, size=M)
    V_h_noise.append(np.sqrt(-white_noise_h*np.log(randoms_u))*np.exp(1j*randoms_theta))

    # new u and theta for R(k,m)
    randoms_u = np.random.uniform(0, 1, size=M)
    randoms_theta = -np.pi + 2*np.pi*np.random.uniform(0,1,size=M)
    R = np.sqrt(-Svv_turb*np.log(randoms_u))*np.exp(1j*randoms_theta)

    V_fh = np.sqrt(zdr_turb)*(rho_hv_turb * V_fv + np.sqrt(1 - rho_hv_turb**2 ) * R) * np.exp(1j * delta_hv_turb)
    V_h.append(V_fh)

V_v = np.array(V_v)
V_h = np.array(V_h)
V_v_noise = np.array(V_v_noise)
V_h_noise = np.array(V_h_noise)

# complex numbers that represent the simulation of the noisy I and Qs in the frequency domain
Svv = V_v+V_v_noise
Shh = V_h+V_h_noise

# In[Calculate the spectral polarimetric variables]
Svv_conj = np.conjugate(Svv)
Shh_conj = np.conjugate(Shh)

a_V = Svv*Svv_conj
a_H = Shh*Shh_conj

Sx = Shh*Svv_conj

nominator = (Sx).mean(axis=0)
a = (np.abs(Shh)**2).mean(axis=0)
b = (np.abs(Svv)**2).mean(axis=0)
variable = nominator/np.sqrt(a*b)

rho_hv_2 = np.abs(variable)
delta_hv_2 = np.angle(variable)
zdr_ratio_2 = a_H.mean(axis=0)/a_V.mean(axis=0)

plt.plot(velocities_interpolated, delta_hv_2*180/np.pi)
plt.xlim(0,10)
plt.ylim(-5,10)
plt.show()

# In[Plot the polarimetric variables]

# Ideal Variables
title = None # "Ideal Variables"

fast_plotting(velocities_interpolated, interpolated_delta*180/np.pi, color='b', xlabel="Velocity (m/s)", ylabel="$δ_{hv}$ (degrees)",
              title=title, save_dir=dir_out, filename="delta_ideal.png", ylimits=(-5, 10))

fast_plotting(velocities_interpolated, 10*np.log10(interpolated_zdr), color='b', xlabel="Velocity (m/s)", ylabel="$Z_{DR}$ (dB)",
              title=title, save_dir=dir_out, filename="zdr_ideal.png", ylimits=(-1.5, 1.5))

fast_plotting(velocities_interpolated, interpolated_rho, color='b', xlabel="Velocity (m/s)", ylabel="$ρ_{hv}$",
              title=title, save_dir=dir_out, filename="rho_ideal.png", ylimits=(0.95, 1.02))

# Turbulent Variables
title = None # "Ideal Variables"

fast_plotting(velocities_interpolated, delta_hv_turb*180/np.pi, color='g', xlabel="Velocity (m/s)", ylabel="$δ_{hv}$ (degrees)",
              title=title, save_dir=dir_out, filename="delta_turbulent.png", ylimits=(-5, 10))

fast_plotting(velocities_interpolated, 10*np.log10(zdr_turb), color='g', xlabel="Velocity (m/s)", ylabel="$Z_{DR}$ (dB)",
              title=title, save_dir=dir_out, filename="zdr_turbulent.png", ylimits=(-1.5, 1.5))

fast_plotting(velocities_interpolated, rho_hv_turb, color='g', xlabel="Velocity (m/s)", ylabel="$ρ_{hv}$",
              title=title, save_dir=dir_out, filename="rho_turbulent.png", ylimits=(0.95, 1.02))

# Noise and Turbulence
title = None # "Noisy & Turbulent Variables"
fast_plotting(velocities_interpolated, delta_hv_2*180/np.pi, color='r', xlabel="Velocity (m/s)", ylabel="$δ_{hv}$ (degrees)",
              title=title, save_dir=dir_out, filename="delta_noise_turb.png", ylimits=(-5, 10))

fast_plotting(velocities_interpolated, 10*np.log10(zdr_ratio_2), color='r', xlabel="Velocity (m/s)", ylabel="$Z_{DR}$ (dB)",
              title=title, save_dir=dir_out, filename="zdr_noise_turb.png", ylimits=(-1.5, 1.5))

fast_plotting(velocities_interpolated, rho_hv_2, color='r', xlabel="Velocity (m/s)", ylabel="$ρ_{hv}$",
              title=title, save_dir=dir_out, filename="rho_noise_turb.png", ylimits=(0.95, 1.02))


