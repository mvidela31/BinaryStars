#include <cmath>
#include <ostream>
#include <stan/math/rev/core.hpp>

namespace default_model_namespace {

	using namespace stan::math;

	inline double kepler_eq(double M, double e, std::ostream* pstream) {
		double En;
		double E = M;
		for (int n = 0; n < 200; ++n) {
			// Newton's method iteration
			En = E - (E - e * sin(E) - M) / (1.0 - e * cos(E));
			// Convergence check
			if (fabs((En - E) / En) <= 1.234e-10) {
				return En;
			}
			E = En;
		}
		// Return the best estimate despite non-convergence
		return En;
	}

	inline var kepler_eq(const var& M_var, const var& e_var, std::ostream* pstream) {
		// Extract values from inputs
	  	double M = M_var.val();
		double e = e_var.val();
		double E = kepler_eq(M, e, pstream);
		// Compute the partial derivatives
		double dE_dM = 1.0 / (1.0 - e * cos(E));
		double dE_de = sin(E) * dE_dM;
		// Return the new parameter with values of the function and gradients
		return make_callback_var(E, [M_var, e_var, dE_dM, dE_de](auto& vi) {
			M_var.adj() += vi.adj() * dE_dM;
			e_var.adj() += vi.adj() * dE_de;
		});
	}

	inline var kepler_eq(double M, const var& e_var, std::ostream* pstream) {
		// Extract values from inputs
	  	double e = e_var.val();
		double E = kepler_eq(M, e, pstream);
		// Compute the partial derivatives
		double dE_de = sin(E) / (1.0 - e * cos(E));
		// Return the new parameter with values of the function and gradients
		return make_callback_var(E, [e_var, dE_de](auto& vi) {
			e_var.adj() += vi.adj() * dE_de;
		});
	}

	inline var kepler_eq(const var& M_var, double e, std::ostream* pstream) {
		// Extract values from inputs
	  	double M = M_var.val();
		double E = kepler_eq(M, e, pstream);
		// Compute the partial derivatives
		double dE_dM = 1.0 / (1.0 - e * cos(E));
		// Return the new parameter with values of the function and gradients
		return make_callback_var(E, [M_var, dE_dM](auto& vi) {
			M_var.adj() += vi.adj() * dE_dM;
		});
	}
}
