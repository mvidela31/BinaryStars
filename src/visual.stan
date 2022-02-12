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
}
data {
    int<lower=0>   N_as;       // number of samples
    vector[N_as]   t_as;       // epochs (years)
    vector[N_as]   x_obs;      // x coordinate observation
    vector[N_as]   y_obs;      // y coordinate observation
    vector[N_as]   x_err;      // x coordinate error
    vector[N_as]   y_err;      // y coordinate error
}
transformed data {
    real   min_t     = min(t_as);      // first epoch
    vector[N_as] t0  = t_as - min_t;   // normalized as epochs
}
parameters {
    real<lower=0, upper=1>        T0;      // normalized time of periastron passage (T - t0) / P
    real                          log_P;   // ln(period)
    real<lower=0, upper=1>        e;       // eccentricity
    real<lower=0>                 a;       // semimajor axis (seconds of arc)
    real<lower=0, upper=2*pi()>   w;       // argument of periapsis (rad)
    real<lower=0, upper=pi()>     Omega;   // longitude of the ascending node (rad)
    real<lower=0, upper=pi()>     i;       // inclination (rad)
}
transformed parameters {
    real P = exp(log_P);   // period
}
model {
    // Variables declaration
    real T; vector[N_as] pos[2];
    // Period
    T = T0 * P;
    // Orbit
    pos = orbit(N_as, t0, T, P, e, a, w, Omega, i);
    x_obs ~ normal(pos[1], x_err);
    y_obs ~ normal(pos[2], y_err);
}
generated quantities {
    real w_deg     = w * 180 / pi();         // argument of periapsis (degrees)
    real Omega_deg = Omega * 180 / pi();     // longitude of the ascending node (degrees)
    real i_deg     = i * 180 / pi();         // inclination (degrees)
    real T         = T0 * P + min_t;         // unnormalized time of periastron passage
}
