#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 David Lilien and Nicholas Rathmann <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Fit constants similar to Richards et al, 2021.

We use exponential in temperature for migration and linear for rotation recrystallization.

Fit is dimensionless, so then  Γ = tilde{Γ} * strainrate, ditto for lambda
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

year_in_seconds = 365.25 * 24.0 * 60.0 * 60.0

# Data from Richards et al., table 1.
# Independent variables for regression
measurement_temps_celcius = np.array([-30.0, -13.6, -10.2, -9.5, -30.3, -7.0, -5.5])
measurement_gamma_dots = np.array([1.26e-5, 7.5e-6, 1.03e-5, 1.66e-4, 1.33e-4, 1.14e-6, 1.26e-4])

measurement_temps = measurement_temps_celcius + 273.15
temp_covariances_percent = np.array([0.0, 0.0, 8.9, 3.0, 1.2, 0.0, 3.5])
temp_sigma = np.sqrt(np.abs(temp_covariances_percent / 100.0 * measurement_temps_celcius))

# dependent variables
measurement_lambda_tildes = np.array([0.173, 0.198, 0.126, 0.343, 0.153, 0.139, 0.178])
measurement_beta_tildes = np.array([0.62, 4.25, 5.92, 2.75, 0.763, 4.12, 5.51])

# measurement_lambdas = measurement_lambda_tildes * measurement_gamma_dots
# measurement_betas = measurement_beta_tildes * measurement_gamma_dots


# Define the functions that will be used for fitting. Just use curve_fit in case we go more complex in the future.
def arrhenius_relation(A, T=measurement_temps):
    return A[0] * np.exp(-A[1] / T)


def linear_relation(A, T=measurement_temps_celcius):
    return A[1] * T + A[0]


def mig(A):
    return measurement_beta_tildes - arrhenius_relation(A, measurement_temps)


migration_fits, migration_cov = leastsq(mig, np.array([1.0, 1.0]), maxfev=10000)
migration_fits_lin, migration_cov_lin = leastsq(lambda x: measurement_beta_tildes - linear_relation(x), np.array([1.0, 1.0]))
rotation_fits, rotation_cov = leastsq(lambda x: measurement_lambda_tildes - linear_relation(x), np.array([1.0, 1.0]))

print("Migration recrystallization: Γ_0 = {:e} x exp(-{:e} / T)".format(migration_fits[0], migration_fits[1]))
print("Migration recrystallization: Γ_0 = {:e} x T + {:e}".format(migration_fits_lin[1], migration_fits_lin[0]))
print("Rotation recrystallization: Λ_0 = {:e} x T + {:e}".format(rotation_fits[1], rotation_fits[0]))

print("Zeros at {:e}, {:e}".format(-migration_fits_lin[0] / migration_fits_lin[1], -rotation_fits[0] / rotation_fits[1]))
print("At 0C {:e}, {:e}, {:e}".format(arrhenius_relation(migration_fits, 273.15), linear_relation(migration_fits_lin, 0.0), linear_relation(rotation_fits, 0.0)))
print("At -10C {:e}, {:e}, {:e}".format(arrhenius_relation(migration_fits, 263.15), linear_relation(migration_fits_lin, -10.0), linear_relation(rotation_fits, -10.0)))
print("At -25C {:e}, {:e}, {:e}".format(arrhenius_relation(migration_fits, 248.15), linear_relation(migration_fits_lin, -25.0), linear_relation(rotation_fits, -25.0)))

strainrate = 0.1
print("Assuming a characteristic strain of {:e}/a".format(strainrate))
print("At 0C Γ_0 = {:e}, Λ_0 = {:e}".format(strainrate * arrhenius_relation(migration_fits, 273.15), strainrate * linear_relation(rotation_fits, 0.0)))
print("At -5C Γ_0 = {:e}, Λ_0 = {:e}".format(strainrate * arrhenius_relation(migration_fits, 268.15), strainrate * linear_relation(rotation_fits, -5.0)))
print("At -30C Γ_0 = {:e}, Λ_0 = {:e}".format(strainrate * arrhenius_relation(migration_fits, 243.15), strainrate * linear_relation(rotation_fits, -30.0)))

# Reproduce figure 4
dummy_temps_c = np.linspace(-31.0, 0.0)
dummy_temps = dummy_temps_c + 273.15
fig, ax = plt.subplots()
ax.plot(measurement_temps_celcius, measurement_beta_tildes, color="firebrick", marker="^", linestyle="none")
ax.plot(measurement_temps_celcius, measurement_lambda_tildes, color="gold", marker="s", linestyle="none")

ax.plot(dummy_temps_c, arrhenius_relation(migration_fits, dummy_temps), color="firebrick", linestyle="dashed")
ax.plot(dummy_temps_c, linear_relation(migration_fits_lin, dummy_temps_c), color="firebrick", linestyle="dashed")
ax.plot(dummy_temps_c, linear_relation(rotation_fits, dummy_temps_c), color="gold", linestyle="dotted")
plt.show()
