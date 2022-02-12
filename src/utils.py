import numpy as np
import pandas as pd


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)

def jdays2years(epochs):
    if (epochs.shape[0] > 0 and epochs[0] > 30000.0000):
        return 1900.0 + (epochs - 15020.31352) / 365.242198781
    else:
        return epochs
    
def kepler(M, e, eps=1e-4):
    En  = M
    Ens = En - (En - e * np.sin(En)- M) / (1 - e * np.cos(En))
    k = 1
    while (abs(Ens - En) > eps):
        En = Ens
        Ens = En - (En - e * np.sin(En) - M) / (1 - e * np.cos(En))
        k += 1
        if (k > 400000):
            print('Kepler Eq (%f, %f) not solved in less than 400000 iterations\n', M, e)
            E = np.nan
            return E
    E = Ens
    return E

def getOrbit(t, T, P, e, a, w, Omega, i):
    N = len(t)
    E = np.empty(N)
    # Iterate over epochs
    for j in range(N):
        # Mean anomaly
        M = 2 * np.pi * (t[j] - T) / P
        # Eccentric anomaly
        E[j] = kepler(M, e)
    # Thiele-Innes Constants
    A = a * ( np.cos(w) * np.cos(Omega) - np.sin(w) * np.sin(Omega) * np.cos(i))
    B = a * ( np.cos(w) * np.sin(Omega) + np.sin(w) * np.cos(Omega) * np.cos(i))
    F = a * (-np.sin(w) * np.cos(Omega) - np.cos(w) * np.sin(Omega) * np.cos(i))
    G = a * (-np.sin(w) * np.sin(Omega) + np.cos(w) * np.cos(Omega) * np.cos(i))
    # Auxiliary normalized coordinates
    x = np.cos(E) - e
    y = np.sin(E) * np.sqrt(1 - e ** 2)
    # Apparent Orbit
    X = A * x + F * y
    Y = B * x + G * y
    return (X, Y)

def getRV(t, T, P, e, a, w, i, V0, plx, q, K1=None, K2=None):
    N = len(t)
    E = np.empty(N)
    # Iterate over epochs
    for j in range(N):
        # Mean anomaly
        M = 2 * np.pi * (t[j] - T) / P
        # Eccentric anomaly
        E[j] = kepler(M, e)
    # True anomaly
    nu = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))
    # Amplitudes
    convCoeff = 149597870.660 / (365.25 * 86400.0)
    if (K1 is None):
        K1 = 2 * np.pi * (a / plx) * (q / (1 + q)) * convCoeff * np.sin(i) / (P * np.sqrt(1 - e ** 2))
    if (K2 is None):
        K2 = 2 * np.pi * (a / plx) * (1 / (1 + q)) * convCoeff * np.sin(i) / (P * np.sqrt(1 - e ** 2))
    # Radial velocities
    V1 = V0 + K1 * (np.cos(w + nu) + e * np.cos(w))
    V2 = V0 - K2 * (np.cos(w + nu) + e * np.cos(w))
    return (V1, V2)

def getRVSingle(t, T, P, e, a, w, i, V0, f_plx, K1=None):
    N = len(t)
    E = np.empty(N)
    # Iterate over epochs
    for j in range(N):
        # Mean anomaly
        M = 2 * np.pi * (t[j] - T) / P
        # Eccentric anomaly
        E[j] = kepler(M, e)
    # True anomaly
    nu = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))
    # Amplitude
    convCoeff = 149597870.660 / (365.25 * 86400.0)
    if (K1 is None):
        K1 = 2 * np.pi * a * f_plx * convCoeff * np.sin(i) / (P * np.sqrt(1 - e ** 2))
    # Radial velocity
    V1 = V0 + K1 * (np.cos(w + nu) + e * np.cos(w))
    return V1

def read_in3(filepath, polar_coords=True):
    # Astrometric observations (polar coordinates)
    epoch_as = []
    rho_as   = []
    theta_as = []
    sigma_as = []
    # Spectroscopic observations
    epoch_v1 = []
    obs_v1   = []
    sigma_v1 = []
    epoch_v2 = []
    obs_v2   = []
    sigma_v2 = []
    # Priors values (Gaussian)
    plx_obs = 0
    plx_err = 0
    m1_obs  = 0
    m1_err  = 0
    with open(filepath) as file:
        line = file.readline()
        while(line):
            tokens = line.split()
            if tokens:
                if (tokens[0] != 'C') and (tokens[0] != 'CC'):
                    if 'I1' in tokens:
                        rho_as.append(float(tokens[2]))
                        theta = float(tokens[1])
                        theta_as.append(theta * np.pi / 180.0)
                        epoch_as.append(float(tokens[0]))
                        sigma_as.append(abs(float(tokens[3])))
                    if 'Va' in tokens:
                        epoch_v1.append(float(tokens[0]))
                        obs_v1.append(float(tokens[1]))
                        sigma_v1.append(abs(float(tokens[2])))
                    if 'Vb' in tokens:
                        epoch_v2.append(float(tokens[0]))
                        obs_v2.append(float(tokens[1]))
                        sigma_v2.append(abs(float(tokens[2])))
                    if 'parallax:' in tokens:
                        plx_obs = float(tokens[1])
                        plx_err = abs(float(tokens[2]))
                    if 'm1:' in tokens:
                        m1_obs = float(tokens[1])
                        m1_err = abs(float(tokens[2]))      
            line = file.readline()
    # Polar to cartesian coordinates
    rho_as = np.array(rho_as)
    theta_as = np.array(theta_as)
    X_as, Y_as = pol2cart(rho_as, theta_as)
    # Julian days to Besselian years
    epoch_as = jdays2years(np.array(epoch_as))
    epoch_v1 = jdays2years(np.array(epoch_v1))
    epoch_v2 = jdays2years(np.array(epoch_v2))
    # Dataframes
    df_as = pd.DataFrame(np.vstack([epoch_as, X_as, Y_as, sigma_as, sigma_as]).T, columns=['epoch', 'X', 'Y', 'X_err', 'Y_err'])
    df_v1 = pd.DataFrame(np.vstack([epoch_v1, obs_v1, sigma_v1]).T, columns=['epoch', 'RV', 'err'])
    df_v2 = pd.DataFrame(np.vstack([epoch_v2, obs_v2, sigma_v2]).T, columns=['epoch', 'RV', 'err'])
    df_priors = pd.DataFrame(np.vstack([['plx', 'm1'], [plx_obs, m1_obs], [plx_err, m1_err]]).T, columns=['param', 'mean', 'std'])
    
    return (df_as, df_v1, df_v2, df_priors)
