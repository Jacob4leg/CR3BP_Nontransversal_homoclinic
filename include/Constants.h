#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include "capd/capdlib.h"

namespace constants {
    // Relative mass of the smaller body (Jupiter) in the Sun-Jupiter system.
    inline constexpr long double MU = 0.00095388114032796904L;

    // V0: approximate point on the vertical Lyapunov orbit that has a
    // homoclinic connection (lies in Fix(R), symmetry components zeroed).
    // W0: approximate point of the unstable manifold of V0 whose image under
    // P_y P_x lies in Fix(R) (w~ in Sec. 5); near but not equal to V0.
    inline const capd::LDVector W0{0.9468923401720671061132517L,-4.072102120831082499146823e-24L,0.05316795353478707980175375L,
        -5.553112274845604899656097e-08L,-0.01115319054270743243833971L,-9.199674025000000205149813e-08L};
    inline const capd::LDVector V0{0.94689235111997293925295213L,0,0.053167950300199904636584778L,0,-0.011153205408683909493953946L,0};

    // (tau~, e~) in Sec. 5: approximate (tau, e) of W0, the center of Omega
    inline constexpr long double TAU_TILDE = 9.1996740268429362011e-08L;
    inline constexpr long double E_TILDE = 3.0268216642156945825L;

    // Cone condition parameters for the horizontal/vertical quadratic forms.
    inline constexpr double ALPHA_H = 0.001;
    inline constexpr double ALPHA_V = 0.001;
    inline constexpr double M_H = 1000.0 * 1000.0;
    inline constexpr double M_V = 1000.0 * 1000.0 + 1;

    // Omega = (TAU_TILDE,E_TILDE) + (-DELTA_TAU,DELTA_TAU) x (-DELTA_E,DELTA_E) (Sec. 5);
    // double literals on purpose: the computation was done with these values
    inline constexpr long double DELTA_TAU = 1e-10;
    inline constexpr long double DELTA_E = 1e-8;

    // records per wall of Omega in init_data.txt: 300 x 8000 pieces
    inline constexpr int TAU_DIV = 301;
    inline constexpr int E_DIV = 8001;

    // N = [0,TAU_MAX] x [-TAU_MAX sqrt(ALPHA_H), TAU_MAX sqrt(ALPHA_H)]^3
    inline constexpr double TAU_MAX = 1e-7;
}

#endif
