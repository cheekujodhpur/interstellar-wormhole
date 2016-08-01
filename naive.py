#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Number of simulation points
N = 1000

# Constants
rho = 0.5

# Functions
def nl(theta, phi):
    return -np.sin(theta)*np.cos(phi)

def nph(theta, phi):
    return -np.sin(theta)*np.sin(phi)

def nt(theta, phi):
    return np.cos(theta)

def r(l):
    return np.sqrt(rho*rho + l*l)

def dr_dl(l):
    return l/r(l)

def p_l(theta, phi):
    return nl(theta, phi)

def p_theta(r,theta,phi):
    return r*nt(theta, phi)

def p_phi(r, theta, phi, theta_c):
    return r * np.sin(theta_c) * nph(theta, phi)

def b(r, theta, phi, theta_c):
    return p_phi(r, theta, phi, theta_c)

def B_sq(r, theta, phi):
    M = nt(theta, phi)
    N = nph(theta, phi)
    return r*r*(M*M + N*N)

ct_b = 1
ct_Bsq = 1

# Time for simulation
T = 100
t = np.linspace(0, T, N)

def dy_dt(y, t):
    return np.array([
            y[3],
            y[4]/(r(y[0])**2),
            ct_b/((r(y[0])**2)*(np.sin(y[1])**2)),
            ct_Bsq*((dr_dl(y[0]))/(r(y[0])**3)),
            ((ct_b**2)/(r(y[0])**2))*(np.cos(y[1])/(np.sin(y[1])**3))
            ])


# Trial simulation
theta_cs = 0.0
phi_cs = 0.0

with open("transform.csv", "w") as outfile:
    while theta_cs<=np.pi:
        phi_cs = 0.0
        while phi_cs <= 2*np.pi:

            # Initial
            l = 1.0
            theta = np.pi/2
            phi = 0

            ct_b = b(r(l), theta_cs, phi_cs, theta)
            ct_Bsq = B_sq(r(l), theta_cs, phi_cs)

            y0 = np.array([
                    l,
                    theta,
                    phi,
                    p_l(theta_cs, phi_cs),
                    p_theta(r(l), theta_cs, phi_cs)
                    ])


#            print "Entering for", theta_cs, phi_cs
            result = odeint(dy_dt, y0, t)
#            print "Exit with", result[-1][1], result[-1][2]
        
            outfile.write(str(theta_cs) + "," +
                        str(phi_cs) + "," +
                        str(result[-1][0]) + "," +
                        str(result[-1][1]) + "," +
                        str(result[-1][1]) + "\n")
            
            phi_cs = phi_cs + np.pi/30

        theta_cs  = theta_cs + np.pi/30

