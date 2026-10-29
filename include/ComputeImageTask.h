#ifndef _COMPUTE_IMAGE_TASK_H_
#define _COMPUTE_IMAGE_TASK_H_

#include<iostream>
#include<atomic>
#include<tuple>
#include "capd/capdlib.h"
#include "Constants.h"
#include "cr3bp.h"
#include "MatrixTools.h"
#include "NonRigorousSolver.h"



class ComputeImageTask : public capd::threading::Task {
    typedef capd::Interval Interval;
    typedef capd::LDVector LDVector; typedef capd::IVector IVector;
    typedef capd::LDMatrix LDMatrix; typedef capd::IMatrix IMatrix;

    Interval tau, E; // tau parametrizes the unstable manifold, E is the energy level; this box is one piece of the boundary of Omega or one interior strip
    Interval E0,E_remainder;
    IVector V,V0,V_remainder; // rigorous (x,y,z,dx,dy,dz) box for periodic-orbit fixed point, parameterized by E: V0 is the center, V_remainder the uncertainty around it
    // L: change of basis into the eigendirections; 
    // J: isomorphism into the basis where the final zero-finding argument applies; 
    // all three are fixed once in main() and shared
    // read-only across every task, so they're held by reference instead of copied per task.
    const IMatrix &L, &J, &L_inv;
    
    capd::IMap* E_map = nullptr; // per-thread pooled energy map/Poincare maps, fetched once in run() to avoid cross-thread solver state sharing
    capd::LDMap* LDE_map = nullptr;
    capd::IPoincareMap* pm_x = nullptr;
    capd::IPoincareMap* pm_y = nullptr;
    capd::IPoincareMap* pm_z = nullptr;

    
    // Implicit-function gradient of z(x,y,dx,dy,dz) on the level set E=e, packaged as
    // a change-of-basis matrix: row 2 (z) becomes -dE/d(other coords) / dE/dz, all
    // other rows stay identity. 
    IMatrix get_energy_change_of_basis(const IVector& dE) {
        IVector z_grad = -IVector{dE[0],dE[1],-1,dE[3],dE[4],dE[5]} / dE[2];
        IMatrix T = IMatrix::Identity(6);
        T.row(2) = z_grad;
        return T;
    }

    // Standard 1D Newton: finds approximate z such that E(u0)=e for a 
    // fixed (x,y,dx,dy,dz). Only used to get a good 
    // center z0 for the rigorous checks
    double find_approximate_z(LDVector u0, long double e) {

        long double z_delta;
        do{
            long double f = (*LDE_map)(u0)[0] - e;
            LDVector E_grad = LDE_map->derivative(u0)[0];
            z_delta = f / E_grad[2];
            u0[2] -= z_delta;
        }while (abs(z_delta) > 1e-15);

        return u0[2];
    }

    /**
     First re-centers W[2] at an approximate root z0 at the center (non-rigorous Newton).
     Then checks (C_z_derivative) from Sec. 3.3 on the box W: dE/dz < 0 on W, and
     E(k,W[2].left) > E, E(k,W[2].right) < E for every k in pi(W). Since E decreases in z,
     for every k in pi(W) and e in E there is a unique z(k,e) in W[2].
     Then encloses z(k,E0) for all k in pi(center) by one mean-value step from z0.
     @param [in] W - a box (x,y,z,dx,dy,dz)
     @param [in] center - the (interval) center of the set S, contained in W; its z is ignored
     @return - an enclosure of z(pi(center),E0), together with T over the re-centered W
     */
    std::pair<Interval,IMatrix> get_proper_bound_for_z(IVector W, IVector center) {
        LDVector c0(6);
        for(int i = 0; i < 6; i++) c0[i] = center[i].mid().leftBound();
        double z0 = find_approximate_z(c0,E0.leftBound());

        // the graph of z drifts with the displacement L*eta, while W[2] is built around v(e);
        // re-center W[2] at z0, with a radius wide enough to still cover the original W[2]
        Interval r = capd::max(z0 - W[2].left(), W[2].right() - z0);
        W[2] = z0 + r * Interval(-1,1);

        IVector dE = E_map->derivative(W)[0];
        IVector W_l = W, W_r = W;
        W_l[2] = W[2].left(); W_r[2] = W[2].right();
        if(!(dE[2] < 0 && (*E_map)(W_l)[0] > E && (*E_map)(W_r)[0] < E))
            throw std::runtime_error("Implicit function z(...) is not proven");

        // z(k,E0) and z0 lie in W[2], so by the mean value theorem
        // z(k,E0) = z0 - (E(k,z0) - E0) / dE/dz(k,theta) with (k,theta) in W
        center[2] = z0;
        Interval z = z0 - ((*E_map)(center)[0] - E0) / dE[2];

        return {z, get_energy_change_of_basis(dE)};
    }

    /**
     Encloses Phi_e(eta) = iota_e(pi^(v(e)) + L*eta) for all e in E, v(e) in V and eta in
     the given box (Sec. 4.4), as the set (CAPD tripleton)
         center + (T*L) * eta_remainder + remainder.
     - center: V0 + L*eta0 (eta0 = mid of eta) with its z replaced by an enclosure of
       z(center, E0) from get_proper_bound_for_z;
     - T = D iota: the identity except row 2 = gradient of z, enclosed over W; its column 2
       is dz/de, since slot 2 of the input is E (L acts as the identity on slot 2);
     - eta_remainder: eta - eta0, with slot 2 = E - E0;
     - remainder: T * (V - V0), the spread of v(e) from Step 1, with the z slot set to 0
       (z of the point is not taken from v(e), it follows from iota_e).
     W = V + L*eta, widened in z by z_margin, is the box on which z(k,e) is proven to exist
     and T is enclosed. Used by Step 2 (cone_conditions) and Step 3 (run).
     @param [in] eta - a box (tau,0,0,gamma1,gamma2,gamma3) in the eigenbasis; slot 1 is y,
                       slot 2 must be 0 (E enters through eta_remainder)
     @param [in] z_margin - extra half-width of the z-range of W
     @return - center, T, eta_remainder and remainder of the set
     */
    std::tuple<IVector,IMatrix,IVector,IVector> initial_set(IVector eta, double z_margin) {
        IVector W = V + L * eta;
        W[2] += z_margin * Interval(-1,1);

        IVector eta0 = eta, eta_remainder(6);
        split(eta0,eta_remainder);
        eta_remainder[2] = E_remainder; // E's own uncertainty, around the exact center E0

        IVector center = V0 + L * eta0;
        auto [z,T] = get_proper_bound_for_z(W, center);
        center[2] = z;

        IVector remainder = V_remainder;
        remainder[2] = 0; // z is recomputed in center[2]; slot 2 of T is E, not z
        return {center, T, eta_remainder, T * remainder};
    }

    /**
        The vertical Lyapunov orbit is double symmetric with respect to reflection z->-z and
        R-symmetry. Its point (x,0,z,0,dy,0) in Fix(R) is mapped by P_z into Fix(S o R):
        (x,0,z,0,dy,0) -> P_z(x,0,z,0,dy,0) = (x1,y1,z1,dx1,dy1,dz1) = (x1,0,0,0,dy1,dz1)
        Setting the Poincare section z=0 gives two equations, y1=0 and dx1=0; together with the
        energy constraint E(x,0,z,0,dy,0)=e (here z is also unknown) this gives three equations
        for the three unknowns (x,z,dy). This proves existence and uniqueness of the implicit
        function v(e)=(x,z,0,dy,0) within the box U, for every e in the given interval, via the
        interval Newton operator for parameterized equations.
        @param [in] U - a box (x,y,z,dx,dy,dz) enclosing the periodic orbits for e in the given interval E
        @return - true if existence and uniqueness of v(e) in U is proven for every e in E
     */
    bool prove_periodic_orbit(IVector U) {
        IVector u0 = U;
        IVector U_remainder(6);
        split(u0,U_remainder); // u0 = midpoint, U_remainder = box around it

        IMatrix D(6,6);
        capd::C1HORect2Set S(U);   // full box, to get the Jacobian over the whole box
        capd::C1HORect2Set S0(u0); // midpoint, to get the function value at a single point

        Interval e_val = (*E_map)(u0)[0] - E;     // how far the midpoint's energy is from E
        IVector e_grad = E_map->derivative(U)[0]; // energy gradient over the full box

        IVector Y = (*pm_z)(S,D);  // image of P_z over the full box (also fills D)
        IVector Y0 = (*pm_z)(S0);  // image of P_z at the midpoint

        IVector p0{Y0[1],Y0[3],e_val}; // (y1,dx1,E-e) evaluated at the midpoint — should be ~0

        D = pm_z->computeDP(Y,D); // Jacobian of P_z w.r.t. the initial point
        IMatrix M{{D[1][0],D[1][2],D[1][4]},
                  {D[3][0],D[3][2],D[3][4]},
                  {e_grad[0],e_grad[2],e_grad[4]}}; // Jacobian of (y1,dx1,E) w.r.t. (x,z,dy)

        IVector N = - capd::matrixAlgorithms::gauss(M,p0); // Newton step from the midpoint
        IVector delta{U_remainder[0],U_remainder[2],U_remainder[4]};
        // existence & uniqueness: the Newton step stays inside the original box
        return subset(N,delta);
    }


    /**
     Cone-conditions check (Capinski, Wasieczko-Zajac): verifies that the fixed alpha_h,
     alpha_v make the quadratic forms Q_h, Q_v positive definite for M = D(P_y)^2 restricted
     to the 4 free directions (tau and the 3 center-stable ones gamma), i.e. that
       M^T Q_h M - m_h Q_h > 0   and   M^T Q_v M - m_v Q_v > 0
     hold across a small tau-sweep around V, split into tau_div sub-boxes. This only
     confirms the constants alpha_h, alpha_v, m_h, m_v already chosen in Constants.h
     @param [in] tau_div - number of sub-boxes the small tau-sweep is split into
     @return - true if the cone conditions hold
     */
    bool cone_conditions(int tau_div) {
        // Total tau-width (in the eigenbasis) chosen, such that
        // tau~ + DELTA_TAU < sqrt(1 - alpha_v) * tau_max, cut into tau_div sub-boxes
        // so each one's image stays tight enough for the cone check to pass.
        double tau_max = constants::TAU_MAX;
        double tau_delta = tau_max / tau_div;

        Interval alpha_h = constants::ALPHA_H, m_h = constants::M_H;
        Interval alpha_v = constants::ALPHA_V, m_v = constants::M_V;

        
        // Coordinate 0 is the unstable direction tau, coordinates 1-3 the center-stable
        // ones gamma: Q_h = alpha_h tau^2 - |gamma|^2, Q_v = tau^2 - alpha_v |gamma|^2.
        IMatrix Q_h = -IMatrix::Identity(4);
        Q_h[0][0] = alpha_h;
        IMatrix Q_v = IMatrix::Identity(4);
        Q_v[1][1] = -alpha_v; Q_v[2][2] = -alpha_v; Q_v[3][3] = -alpha_v;

        // Q_h(v)>=0 means alpha_h*v0^2 >= v1^2+v2^2+v3^2, i.e. |v1,2,3| <= sqrt(alpha_h)*v0;
        // gamma is exactly that bound for v0 up to tau_max, keeping w (below) inside the cone.
        Interval gamma = tau_max * sqrt(alpha_h) * Interval(-1,1);
        Interval tau_i(0, tau_delta);

        // One sub-box per iteration: re-derive its implicit z and energy-consistent
        // basis T, propagate P_y's derivative seeded with that basis, then express it
        // in the eigenbasis (L_inv) and fold it into the running hull D_total.
        IMatrix D_total(5,5);
        for(int i = 0; i < tau_div; i++) {
            IVector eta{tau_i,0,0, gamma,gamma,gamma}; // this sub-box's offset in the eigenbasis
            auto [center, T, eta_remainder, remainder] = initial_set(eta, 2.0e-8);

            IMatrix T_total = T * L; // eigenbasis shift, then z corrected for the energy constraint

            // T_total is also the initial condition for the variational equation, since the
            // propagated derivative needs to start from the energy-consistent basis.
            capd::C1HORect2Set S(capd::C1Rect2Set::C0BaseSet(center, T_total, eta_remainder, remainder),
                                  capd::C1Rect2Set::C1BaseSet(T_total));
            
            IMatrix D(6,6);
            IVector y = (*pm_y)(S,D,2);
            if(!(y[2] > 0)) throw std::runtime_error("pi_z F(Phi_e(N)) > 0 is not verified");
            D = pm_y->computeDP(y,D);
            
            // Df_e = pi DF T in the eigenbasis. No T^{-1} on the left: T is the identity except
            // row 2, and row/column 2 (E) is dropped below, as is y; L_inv does not mix row 2
            // into the others.
            D = L_inv * D;
            
            IMatrix D_reduced = matrix_erase_cord(D,1); // drop the y-coordinate (always 0 on this section)
            D_total = (i == 0) ? D_reduced : intervalHull(D_reduced,D_total);

            tau_i += tau_delta;
        }
        
        IMatrix M = matrix_erase_cord(D_total,1);

        IMatrix M_T = M;
        M_T.transpose();
        IMatrix V_h = M_T * Q_h * M - m_h * Q_h;
        IMatrix V_v = M_T * Q_v * M - m_v * Q_v;


        // Gershgorin circle criterion: a sufficient condition for
        // positive-definiteness -- each diagonal entry must exceed the sum of
        // absolute values of the rest of its row.
        auto is_positive_definite = [](const IMatrix& Q) {
            for(int i = 0; i < 4; i++) {
                Interval R = 0;
                for(int j = 0; j < 4; j++) if(j != i) R += capd::abs(Q[i][j]);
                if(!(Q[i][i] - R > 0)) return false;
            }
            return true;
        };

        if(is_positive_definite(V_h) && is_positive_definite(V_v))
            return true;
        return false;
    }

public:
    enum Status {NEW, COMPLETED, FAILED};

    IVector result; // final (dx,dz)-residual hull, read by main.cpp after task->join()
    Status status = NEW; // outcome reported back to main.cpp's task loop
    bool compute_boundary;

    // shared across all tasks/threads; used by main.cpp to report progress
    static std::atomic<int> completed_count;

    // Calibrates the input box: widens x,z,dy by their own remainder (plus a tiny
    // inflation first if E is degenerate) so the interval-Newton checks below have
    // a genuine box to work with, rather than a single point.
    ComputeImageTask(Interval tau, Interval E, IVector V_in,
            const IMatrix& J, const IMatrix& L, const IMatrix& L_inv, bool compute_boundary) :
            tau(tau), E(E), V(V_in),J(J),L(L), L_inv(L_inv),
            compute_boundary(compute_boundary) {

        result = IVector(2);
        E.split(E0,E_remainder);

        // A degenerate (point) energy interval implies x,z,dy are points too,
        // which would zero out their remainders below; inflate them first so they don't.
        bool degenerate_E = (E.left() == E.right());

        auto inflate_coordinate = [&](int idx) {
            Interval v = V[idx];
            if(degenerate_E) v += 1e-14 * Interval(-1,1);
            Interval v0, v_remainder;
            v.split(v0,v_remainder);
            V[idx] = v + 0.3 * v_remainder;
        };
        inflate_coordinate(0); // x
        inflate_coordinate(2); // z
        inflate_coordinate(4); // dy

        V0 = V; V_remainder = IVector(6);
        split(V0,V_remainder);
    }

    // the calibrated box V (after the inflation above), read by main.cpp for the hull check
    const IVector& get_V() const { return V; }
    

    // Entry point for one piece of Omega: prove the periodic orbit, 
    // check the cone conditions (strips only), compute the rigorous image of this
    // box along the unstable manifold, and (boundary pieces only) test for a zero via the
    // Brouwer-degree argument.
    void run(unsigned threadID) {
        Interval returnTime;
        IVector zero_vector{0,0,0,0,0,0};
        IVector zero{0,0};

        try {
            pm_x = &CR3BPPool::pool_X->getSolver(threadID).pm;
            pm_y = &CR3BPPool::pool_Y->getSolver(threadID).pm;
            pm_z = &CR3BPPool::pool_Z->getSolver(threadID).pm;
            E_map = &CR3BPPool::pool_E->getSolver(threadID);
            LDE_map = &CR3BPPool::pool_LDE->getSolver(threadID);

            // check if V encloses v(e) over the energy interval E
            if(!prove_periodic_orbit(V)) throw std::runtime_error("Proper enclosure of fixed points is not proven");

            // check if F = (P_y)^2 satisfies proper cone conditions on an explicit region;
            // they depend on E only, so they are checked on the interior strips, which cover E
            // (and whose boxes V contain the boxes V of the boundary pieces)
            if(!compute_boundary)
                if(!cone_conditions(5)) throw std::runtime_error("Some quadratic forms are not positive definite");
            
            // same cone-width formula as inside cone_conditions, but here for this box's own
            // (narrow) tau-interval -- cone_conditions itself checks a wider neighborhood (tau_max)
            Interval gamma = sqrt(constants::ALPHA_H) * tau * Interval(-1,1);

            IVector eta{tau,0,0,gamma,gamma,gamma};
            auto [center, T, eta_remainder, remainder] = initial_set(eta, 5e-9);
            capd::C0HOTripletonSet S(center, T * L, eta_remainder, remainder);
            
            // J scales the Poincare image so it lands in the +-1 square
            (*pm_x)(S);
            IVector res = (*pm_y)(S,zero_vector,J,returnTime);
            

            // B rescales Omega onto [-1,1]^2 (Sec. 5)
            IMatrix B = {{1/constants::DELTA_TAU,0},{0,1/constants::DELTA_E}};

            // res1 is the image of the affine map B(zeta - zeta~), res2 the enclosure of G;
            // their hull contains H([0,1], piece). On boundary pieces we check that it misses zero.
            IVector res1 = B * IVector{tau - constants::TAU_TILDE, E - constants::E_TILDE};
            IVector res2 = IVector{res[3],res[5]};
            result = capd::vectalg::intervalHull(res1,res2);
            
            if(compute_boundary)
                if(subset(zero,result)){
                    std::cout << result << std::endl;
                    throw std::runtime_error("Solution contains zero");
                } 

            status = COMPLETED;
        } catch(const std::exception& exc) {
            std::cout << exc.what() << std::endl;
            status = FAILED;
        }
        ++completed_count;
    }
};

std::atomic<int> ComputeImageTask::completed_count{0};

/**
 Checks that the interval matrix [D_(x,z,dy) F] of the Step 1 equation (y1, dx1, E - e) is
 regular on the box U, i.e. every matrix in it is invertible. Same matrix as M in
 prove_periodic_orbit (it does not depend on e). Regularity is checked by the interval
 Gauss elimination. Then F(., e) is injective on U (U is convex), so it has at most one
 zero in U for every e.
 @param [in] U - a box (x,y,z,dx,dy,dz), here the hull of the boxes V of all tasks
 @return - true if the interval matrix is regular
 */
inline bool is_jacobian_regular(const capd::IVector& U, capd::IPoincareMap& pm_z, capd::IMap& E_map) {
    capd::IMatrix D(6,6);
    capd::C1HORect2Set S(U);
    capd::IVector Y = pm_z(S,D);
    D = pm_z.computeDP(Y,D);
    capd::IVector e_grad = E_map.derivative(U)[0];

    capd::IMatrix M{{D[1][0],D[1][2],D[1][4]},
                    {D[3][0],D[3][2],D[3][4]},
                    {e_grad[0],e_grad[2],e_grad[4]}};
    try {
        capd::matrixAlgorithms::gauss(M,capd::IVector(3));
    } catch(const std::runtime_error&) {
        return false;
    }
    return true;
}

#endif