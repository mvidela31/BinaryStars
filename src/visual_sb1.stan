functions {
    real kepler_eq(real M, real e);
    vector[] orbit(int N, vector t, real T, real P, real e, real a, real w, real Omega, real i) {
        // Variables declaration
        real M; real A; real B; real F; real G;
        vector[N_as] E; vector[N_as] x; vector[N_as] y; vector[N_as] pos[2];
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
    vector radial_velocity(int N, vector t, real T, real P, real e, real a, real w, real i, real V0, real alpha) {
        // Variables declaration
        real M; real K; vector[N_as] E; vector[N_as] nu; vector[N_as] V;
        real convCoeff = 149597870.660 / (365.25 * 86400.0); 
        // Iterate over epochs
        for (j in 1:N) {
            // Mean anomaly
            M = 2 * pi() * (t[j] - T) / P;
            // Eccentric anomaly
            E[j] = kepler_eq(M, e);
        }
        // True anomaly
        nu = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
        // Amplitude
        K = 2 * pi() * a * alpha * convCoeff * sin(i) / (P * sqrt(1 - e^2));
        // Radial velocity
        V = V0 + K * (cos(w + nu) + e * cos(w));
        return V;
    }
}
data {
    int<lower=0>    N_as;      // number of as samples
    int<lower=0>    N_v1;      // number of rv1 samples
    vector[N_as]    t_as;      // epochs (years)
    vector[N_v1]    t_v1;      // epochs (years)
    vector[N_as]    x_obs;     // x coordinate observation
    vector[N_as]    y_obs;     // y coordinate observation
    vector[N_as]    x_err;     // x coordinate uncertainties
    vector[N_as]    y_err;     // y coordinate uncertainties
    vector[N_v1]    v1_obs;    // rv1 observations
    vector[N_v1]    v1_err;    // rv1 observation uncertainties
}
transformed data {
    real       min_t = fmin(min(t_as), min(t_v1));    // first epoch
    vector[N_as]  t0 = t_as - min_t;                  // normalized as epochs
    vector[N_v1] t10 = t_v1 - min_t;                  // normalized RV1 epochs
}
parameters {
    real<lower=0, upper=1>        T0;      // normalized time of periastron passage (T - t0) / P [yr]
    real                          log_P;   // ln(period [yr])
    real<lower=0, upper=1>        e;       // eccentricity
    real<lower=0>                 a;       // semimajor axis [seconds of arc]
    real<lower=0, upper=2*pi()>   w;       // argument of periapsis [rad]
    real<lower=0, upper=2*pi()>   Omega;   // longitude of the ascending node [rad]
    real<lower=0, upper=pi()>     i;       // inclination [rad]
    real                          V0;      // velocity of center of mass [km/s]
    real<lower=0, upper=200>      f_plx;   // f/plx = q/(1+q)/plx [parsecs]
}
transformed parameters {
    real P   = exp(log_P);    // period [yr]
}
model {
    // Variables declaration
    real T; vector[N_as] pos[2]; vector[N_v1] V1;
    // Period
    T = T0 * P;
    // Orbit
    pos = orbit(N_as, t0, T, P, e, a, w, Omega, i);
    x_obs ~ normal(pos[1], x_err);
    y_obs ~ normal(pos[2], y_err);
    // Primary object radial velocity
    V1 = radial_velocity(N_v1, t10, T, P, e, a, w, i, V0, f_plx);
    v1_obs ~ normal(V1, v1_err);
}
generated quantities {
    real w_deg     = w * 180 / pi();       // argument of periapsis [degrees]
    real Omega_deg = Omega * 180 / pi();   // longitude of the ascending node [degrees]
    real i_deg     = i * 180 / pi();       // inclination [degrees]
    real T         = T0 * P + min_t;       // unnormalized time of periastron passage [yr]
}
