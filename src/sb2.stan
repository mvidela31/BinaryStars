functions {
    real kepler_eq(real M, real e);
    vector radial_velocity(int N, vector t, real T, real P, real e, real asin_i, real w, real V0, real plx, real q, int primary) {
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
            K = 2 * pi() * (asin_i / plx) * (q / (1+q)) * convCoeff / (P * sqrt(1 - e^2));
            // Radial velocity
            V = V0 + K * (cos(w + nu) + e * cos(w));
        }
        // Companion object radial velocity
        else {
            // Amplitude
            K = 2 * pi() * (asin_i / plx) * (1 / (1+q)) * convCoeff / (P * sqrt(1 - e^2));
            // Radial velocity
            V = V0 - K * (cos(w + nu) + e * cos(w));
        }
        return V;
    }
}
data {
    int<lower=0>   N_v1;       // number of samples
    int<lower=0>   N_v2;       // number of samples
    vector[N_v1]   t_v1;       // epochs (years)
    vector[N_v2]   t_v2;       // epochs (years)
    vector[N_v1]   v1_obs;     // rv1 observations
    vector[N_v1]   v1_err;     // rv1 observation uncertainties
    vector[N_v2]   v2_obs;     // rv2 observations
    vector[N_v2]   v2_err;     // rv2 observation uncertainties
}
transformed data {
    real   min_t     = fmin(min(t_v1), min(t_v2));   // first epoch
    vector[N_v1] t10 = t_v1 - min_t;                 // normalized rv1 epochs
    vector[N_v2] t20 = t_v2 - min_t;                 // normalized rv2 epochs  
}
parameters {
    real<lower=0, upper=1>        T0;      // normalized time of periastron passage (T - t0) / P
    real                          log_P;   // ln(period)
    real<lower=0, upper=1>        e;       // eccentricity
    real                          asin_i;  // semimajor axis * sin(inclination) (seconds of arc)
    real<lower=0, upper=2*pi()>   w;       // argument of periapsis (rad)
    real                          V0;      // velocity of center of mass (km/s)
    real<lower=0>                 plx;     // parallax (seconds of arc)
    real<lower=0, upper=1>        q;       // mass ratio
}
transformed parameters {
    real P = exp(log_P);   // period
}
model {
    // Variables declaration
    real T; vector[N_v1] V1; vector[N_v2] V2;
    // Period
    T = T0 * P;
    // Primary object radial velocity
    V1 = radial_velocity(N_v1, t10, T, P, e, asin_i, w, V0, plx, q, 1);
    v1_obs ~ normal(V1, v1_err);
    // Companion object radial velocity
    V2 = radial_velocity(N_v2, t20, T, P, e, asin_i, w, V0, plx, q, 0);
    v2_obs ~ normal(V2, v2_err);
}
generated quantities {
    real w_deg     = w * 180 / pi();                      // argument of periapsis (degrees)
    real T         = T0 * P + min_t;                      // unnormalized time of periastron passage
    real plx_mas   = plx * 1000;                          // parallax (miliseconds of arc)
    real f_plx     = q / (1 + q) * (1 / plx);             // f/plx = q/(1+q)/plx (parsecs)
}
