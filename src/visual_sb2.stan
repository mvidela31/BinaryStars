functions {
    real kepler_eq(real M, real e);
    vector[] orbit(int N, vector t, real T, real P, real e, real a, real w, real Omega, real i) {
        // Variables declaration
        real M; real A; real B; real F; real G;
        vector[N] E; vector[N] x; vector[N] y; vector[N] pos[2];
        // Iterate over epochs
        for (j in 1:N) {
            // Mean anomaly
            M = 2 * pi() * (t[j] - T) / P;
            // Eccentric anomaly
            E[j] = kepler_eq(M, e);
        }
        // Auxiliary normalized coordinates
        x = cos(E) - e;
        y = sin(E) * sqrt(1 - e^2);
        // Thiele-Innes constants
        A = a * (cos(w) * cos(Omega) - sin(w) * sin(Omega) * cos(i));
        B = a * (cos(w) * sin(Omega) + sin(w) * cos(Omega) * cos(i));
        F = a * (-sin(w) * cos(Omega) - cos(w) * sin(Omega) * cos(i));
        G = a * (-sin(w) * sin(Omega) + cos(w) * cos(Omega) * cos(i));
        // Apparent orbit
        pos[1] = A * x + F * y;
        pos[2] = B * x + G * y;
        return pos;
    }
    vector radial_velocity(int N, vector t, real T, real P, real e, real a, real w, real i, real V0, real plx, real q, int primary) {
        // Variables declaration
        real M; real K; vector[N] E; vector[N] nu; vector[N] V;
        real convCoeff = 149597870.660 / (365.25 * 86400.0); // [AU] to [arcsec] 
        // Iterate over epochs
        for (j in 1:N) {
            // Mean anomaly
            M = 2 * pi() * (t[j] - T) / P;
            // Eccentric anomaly
            E[j] = kepler_eq(M, e);
        }
        // True anomaly
        nu = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
        // Primary object radial velocity
        if (primary) {
            // Amplitude
            K = 2 * pi() * (a / plx) * (q / (1+q)) * convCoeff * sin(i) / (P * sqrt(1 - e^2));
            // Radial velocity
            V = V0 + K * (cos(w + nu) + e * cos(w));
        }
        // Companion object radial velocity
        else {
            // Amplitude
            K = 2 * pi() * (a / plx) * (1 / (1+q)) * convCoeff * sin(i) / (P * sqrt(1 - e^2));
            // Radial velocity
            V = V0 - K * (cos(w + nu) + e * cos(w));
        }
        return V;
    }
}
data {
    int<lower=0>   N_as;       // number of samples
    int<lower=0>   N_v1;       // number of samples
    int<lower=0>   N_v2;       // number of samples
    vector[N_as]   t_as;       // epochs (years)
    vector[N_v1]   t_v1;       // epochs (years)
    vector[N_v2]   t_v2;       // epochs (years)
    vector[N_as]   x_obs;      // x coordinate observation
    vector[N_as]   y_obs;      // y coordinate observation
    vector[N_as]   x_err;      // x coordinate error
    vector[N_as]   y_err;      // y coordinate error
    vector[N_v1]   v1_obs;     // rv1 observations
    vector[N_v1]   v1_err;     // rv1 observation uncertainties
    vector[N_v2]   v2_obs;     // rv2 observations
    vector[N_v2]   v2_err;     // rv2 observation uncertainties
}
transformed data {
    real   min_t     = fmin(min(t_as), fmin(min(t_v1), min(t_v2)));   // first epoch
    vector[N_as] t0  = t_as - min_t;                                  // normalized as epochs
    vector[N_v1] t10 = t_v1 - min_t;                                  // normalized rv1 epochs
    vector[N_v2] t20 = t_v2 - min_t;                                  // normalized rv2 epochs  
}
parameters {
    real<lower=0, upper=1>        T0;      // normalized time of periastron passage (T - t0) / P
    real                          log_P;   // ln(period)
    real<lower=0, upper=1>        e;       // eccentricity
    real<lower=0>                 a;       // semimajor axis (seconds of arc)
    real<lower=0, upper=2*pi()>   w;       // argument of periapsis (rad)
    real<lower=0, upper=2*pi()>   Omega;   // longitude of the ascending node (rad)
    real<lower=0, upper=pi()>     i;       // inclination (rad)
    real                          V0;      // velocity of center of mass (km/s)
    real<lower=0>                 plx;     // parallax (seconds of arc)
    real<lower=0, upper=1>        q;       // mass ratio
}
transformed parameters {
    real P = exp(log_P);   // period
}
model {
    // Variables declaration
    real T; vector[N_as] pos[2]; vector[N_v1] V1; vector[N_v2] V2;
    // Period
    T = T0 * P;
    // Orbit
    pos = orbit(N_as, t0, T, P, e, a, w, Omega, i);
    x_obs ~ normal(pos[1], x_err);
    y_obs ~ normal(pos[2], y_err);
    // Primary object radial velocity
    V1 = radial_velocity(N_v1, t10, T, P, e, a, w, i, V0, plx, q, 1);
    v1_obs ~ normal(V1, v1_err);
    // Companion object radial velocity
    V2 = radial_velocity(N_v2, t20, T, P, e, a, w, i, V0, plx, q, 0);
    v2_obs ~ normal(V2, v2_err);
}
generated quantities {
    real w_deg     = w * 180 / pi();                      // argument of periapsis (degrees)
    real Omega_deg = Omega * 180 / pi();                  // longitude of the ascending node (degrees)
    real i_deg     = i * 180 / pi();                      // inclination (degrees)
    real T         = T0 * P + min_t;                      // unnormalized time of periastron passage
    real m1        = (a / plx)^3 * 1 / (P^2 * (1 + q));   // primary object mass (solar masses)
    real m2        = q * m1;                              // companion object mass (solar masses)
    real plx_mas   = plx * 1000;                          // parallax (miliseconds of arc)
    real f_plx     = q / (1 + q) * (1 / plx);             // f/plx = q/(1+q)/plx (parsecs)
}
