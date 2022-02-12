#include <cmath>
#include <ostream>
#include <stan/math/rev/core.hpp>

namespace visual_sb1_priors_model_namespace {

	using namespace stan::math;

	inline double kepler_eq(double M, double e, std::ostream* pstream) {
		double E0 = M;
		double E = M;
		double g, gp;
		for (int i = 0; i < 200; ++i) {
			g = E0 - e * sin(E0) - M;
			gp = 1.0 - e * cos(E0);
			E = E0 - g / gp;
			// Convergence check
			if (fabs((E - E0) / E) <= 1.234e-10) {
				return E;
			}
			E0 = E;
		}
		// Return the best estimate despite non-convergence
		return E;
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