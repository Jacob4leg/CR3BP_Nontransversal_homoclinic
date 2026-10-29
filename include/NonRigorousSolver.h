#include<iostream>
#include<fstream>
#include <iomanip>
#include "capd/capdlib.h"
#include "Constants.h"
#include "cr3bp.h"
#include "MatrixTools.h"

#ifndef _NON_RIGOROUS_SOLVER_H_
#define _NON_RIGOROUS_SOLVER_H_


struct NonRigorousSolver
{
    typedef capd::LDVector LDVector;
    typedef capd::LDMatrix LDMatrix;

    LDVector w0 = constants::W0;
    LDVector v0 = constants::V0;

    long double tau_0; long double tau_radius;
    long double e0; long double e_radius;
    int tau_div; int e_div;
    long double tau_delta; long double E_delta;

    CR3BP<long double> vf;

    NonRigorousSolver(long double tau_radius, long double e_radius, int tau_div, int e_div) :
        tau_radius(tau_radius), e_radius(e_radius), tau_div(tau_div), e_div(e_div) {

        tau_delta = tau_radius * 2 / (tau_div - 1);
        E_delta = e_radius * 2 / (e_div - 1);
        e0 = constants::E_TILDE;
        tau_0 = constants::TAU_TILDE;
    }

    /**
     * Fixed point of the vertical Lyapunov orbit family at energy e0 + e.
     * @param [in] e - offset from the reference energy level e0
     * @return - approximate fixed point (x,0,z,0,dy,0)
     */
    LDVector get_fixed_point(long double e) {
        return vf.findVerticalLyapunovOrbit(v0,e + e0);
    }

    /**
     * Derivative of two consecutive y-section Poincare returns at fixed point v_e,
     * used to find the dominant (unstable) eigendirection.
     * @param [in] v_e - a fixed point of (P_y)^2
     * @return - the 6x6 derivative of (P_y)^2 at v_e
     */
    LDMatrix derivative(LDVector v_e) {
        LDMatrix D = LDMatrix::Identity(6), T(6,6);
        LDVector u = v_e;
        for(int i = 0; i < 2; i++) {
            u = vf.pm_y(u,T);
            D = vf.pm_y.computeDP(u,T) * D;
        }
        return D;
    }

    /**
     * A point displaced by tau along the dominant eigendirection from the fixed
     * point at energy e - a candidate near the homoclinic point.
     * @param [in] tau - displacement along the unstable eigendirection
     * @param [in] e - offset from the reference energy level e0
     * @return - approximate point on the unstable manifold
     */
    LDVector get_initial_point(long double tau, long double e) {
        LDVector v_e = get_fixed_point(e);
        LDMatrix D = derivative(v_e);
        
        LDVector rV(6),iV(6);
        LDMatrix rVec(6,6), iVec(6,6);
        capd::computeEigenvaluesAndEigenvectors(D,rV,iV,rVec,iVec);

        LDVector v_res = v_e + (tau_0 + tau) * rVec.column(0);

        return v_res;
    }

    /**
     * Isomorphism straightening the Poincare map's image, built from a
     * finite-difference estimate of d(dx,dz)/d(tau,E).
     * @return - the 6x6 isomorphism matrix J_inv (J is its rigorous inverse)
     */
    LDMatrix get_isomorphism() {
        LDVector u1 = get_initial_point(tau_radius,0);
        LDVector u2 = get_initial_point(0,e_radius);

        LDVector Pu1_x = vf.pm_x(u1); LDVector Pu1_y = vf.pm_y(Pu1_x);
        LDVector Pu2_x = vf.pm_x(u2); LDVector Pu2_y = vf.pm_y(Pu2_x);

        // finite-difference estimate of d(dx,dz)/d(tau,E) - the Jacobian of G
        LDMatrix res = LDMatrix::Identity(6);
        res[3][3] = Pu1_y[3]; res[3][5] = Pu2_y[3];
        res[5][3] = Pu1_y[5]; res[5][5] = Pu2_y[5];

        LDVector u = vf.pm_y(vf.pm_x(w0));

        // column 1 (y) absorbs the flow direction at the image of w0, to
        // account for the time-shift freedom along the section
        res.column(1) = vf.vf(u);

        return res;
    }

    /**
     * Real eigenbasis change-of-basis matrix L at energy e, used to set up
     * the cone conditions in the proof.
     * @param [in] e - offset from the reference energy level e0
     * @return - the 6x6 change-of-basis matrix L
     */
    LDMatrix get_change_of_basis(long double e) {
        LDVector v_e = get_fixed_point(e);
        LDMatrix T = LDMatrix::Identity(6);

        LDVector dE = vf.E.derivative(v_e)[0];
        LDVector z_grad = -LDVector{dE[0],dE[1],-1,dE[3],dE[4],dE[5]} / dE[2];

        T.row(2) = z_grad;

        // change basis by T: slot 2 (z) now represents motion along E; T^{-1} on the left
        // would change only row 2, which is dropped below
        LDMatrix D = derivative(v_e) * T;

        // L maps (tau,E,gamma1,gamma2,gamma3) -> (x,E,dx,dy,dz); drop the E slot and y
        // (section Pi_y = {y=0} has no independent y direction) -> 4x4 over (x,dx,dy,dz)
        LDMatrix free_block = matrix_erase_cord(matrix_erase_cord(D,2),1);
        LDMatrix basis = eigenbasis(free_block);            // real eigenbasis on free directions
        LDMatrix L = matrix_add_cord(matrix_add_cord(basis, 1), 1); // reinsert E, then y

        L[2][2] = 1.; // L acts as identity on E; y stays fully zeroed

        return L;
    }

    /**
     * Formats a vector for file output (scientific notation, full precision).
     * @param [in] v - the vector to format
     * @return - a space-separated string of v's components
     */
    std::string vectorToString(LDVector v) {
        std::ostringstream oss;
        oss << std::scientific
            << std::setprecision(std::numeric_limits<long double>::max_digits10);
        for(int i = 0; i < v.dimension(); i++) {
            oss << v[i];
            if(i != v.dimension() - 1) oss << " ";
        }
        return oss.str();
    }

    /**
     * Saves initial data to rigorous solver, one record per line in format
     * tau e x y z dx dy dz
     * where v_e = (x,y,z,dx,dy,dz) is an approximate fixed point on energy level e
     */
    void save_init_data_to_file() {
        // ends of the (tau,e) rectangle; every wall takes its corners from these,
        // so adjacent walls meet exactly
        const long double tau_l = -tau_radius + tau_0, tau_r = tau_radius + tau_0;
        const long double e_l = -e_radius + e0, e_r = e_radius + e0;
        auto tau_at = [&](int i) { return i == tau_div - 1 ? tau_r : tau_l + i * tau_delta; };
        auto e_at = [&](int i) { return i == e_div - 1 ? e_r : e_l + i * E_delta; };

        std::ofstream file("init_data.txt");
        file << std::scientific
             << std::setprecision(std::numeric_limits<long double>::max_digits10);

        long double fixed_point_error;
        long double total_fixed_point_error = 0.;
        long double energy_error;
        long double total_energy_error = 0.;

        capd::LDMaxNorm norm;

        LDVector v_e1 = get_fixed_point(-e_radius);
        fixed_point_error = norm(v_e1 - vf.pm_y(vf.pm_y(v_e1)));
        energy_error = abs(vf.E(v_e1)[0] + e_radius - e0);

        if(total_fixed_point_error < fixed_point_error) total_fixed_point_error = fixed_point_error;
        if(total_energy_error < energy_error) total_energy_error = energy_error;

        LDVector v_e2 = get_fixed_point(e_radius);
        fixed_point_error = norm(v_e2 - vf.pm_y(vf.pm_y(v_e2)));
        energy_error = abs(vf.E(v_e2)[0] - e_radius - e0);

        if(total_fixed_point_error < fixed_point_error) total_fixed_point_error = fixed_point_error;
        if(total_energy_error < energy_error) total_energy_error = energy_error;

        for(int i = 0; i < tau_div; i++)
            file << tau_at(i) << " " << e_l << " " << vectorToString(v_e1) << std::endl;

        file << std::endl;

        for(int i = 0; i < tau_div; i++)
            file << tau_at(i) << " " << e_r << " " << vectorToString(v_e2) << std::endl;

        file << std::endl;

        for(int i = 0; i < e_div; i++) {
            long double e_i = e_at(i);
            LDVector v_e = get_fixed_point(e_i - e0);
            fixed_point_error = norm(v_e - vf.pm_y(vf.pm_y(v_e)));
            energy_error = abs(vf.E(v_e)[0] - e_i);

            if(total_fixed_point_error < fixed_point_error) total_fixed_point_error = fixed_point_error;
            if(total_energy_error < energy_error) total_energy_error = energy_error;

            file << tau_l << " " << e_i << " " << vectorToString(v_e) << std::endl;
        }

        file << std::endl;

        for(int i = 0; i < e_div; i++) {
            long double e_i = e_at(i);
            LDVector v_e = get_fixed_point(e_i - e0);
            file << tau_r << " " << e_i << " " << vectorToString(v_e) << std::endl;
        }

        std::cout << "data saved to init_data.txt with error: " << std::endl;
        std::cout << "fixed_point_error: " << total_fixed_point_error << std::endl;
        std::cout << "energy_error: " << total_energy_error << std::endl;

        file.close();
    }
};

#endif
