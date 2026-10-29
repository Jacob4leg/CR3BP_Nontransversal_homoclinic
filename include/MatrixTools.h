#include<iostream>
#include<tuple>
#include "capd/capdlib.h"


#ifndef _MATRIX_TOOLS_H_
#define _MATRIX_TOOLS_H_

/**
 * Erases row and column corresponding to k-th coordinate. Used mostly to erase y-coordinate
 * in proper places of the code
 * @param [in] M - N x N matrix
 * @param [in] k - coordinate index
 * @return (N-1) x (N-1) matrix without k-th row and column 
 */
template<typename T>
inline capd::vectalg::Matrix<T,0,0> matrix_erase_cord(capd::vectalg::Matrix<T,0,0> M, int k) {
    auto dim = M.dimension();
    int n = std::get<0>(dim);
    assert(n == std::get<1>(dim));

    capd::vectalg::Matrix<T,0,0> M_res(n-1,n-1);

    for(int i = 0; i < n-1; i++) {
        int l = i < k ? i : i+1;
        for(int j = 0; j < n-1; j++) {
            int m = j < k ? j : j+1;
            M_res[i][j] = M[l][m];
        }
    }
    return M_res;
}

/**
 * Adds row and column between k-th and (k+1)-th coordinates. All added entries are equal to 0
 * @param [in] M - N x N matrix
 * @param [in] k - coordinate index
 * @return (N+1) x (N+1) matrix with added additional empty column and row 
 */
template<typename T>
inline capd::vectalg::Matrix<T,0,0> matrix_add_cord(capd::vectalg::Matrix<T,0,0> M, int k) {
    auto dim = M.dimension();
    int n = std::get<0>(dim);
    assert(n == std::get<1>(dim));

    capd::vectalg::Matrix<T,0,0> M_res(n+1,n+1);

    for(int i = 0; i < n+1; i++) {
        if(i == k) continue;
        int l = i < k ? i : i-1;
        for(int j = 0; j < n+1; j++) {
            if(j == k) continue;
            int m = j < k ? j : j-1;
            M_res[i][j] = M[l][m];
        }
    }
    return M_res;
}

/**
 * Computes change of basis L in which A is almost block-diagonal
 * @param [in] A - N x N matrix
 * @return change to eigenbasis
 */
inline capd::LDMatrix eigenbasis(capd::LDMatrix A) {
    std::pair<int,int> dim = A.dimension();
    int n = std::get<0>(dim);
    int m = std::get<1>(dim);
    assert(n == m);

    capd::LDVector rV(n), iV(n);
    capd::LDMatrix rVec(n,n), iVec(n,n);
    capd::alglib::computeEigenvaluesAndEigenvectors(A,rV,iV,rVec,iVec);

    for(int i = 0; i < n; i++) {
        if(iV[i] != 0) {
            rVec.column(i+1) = iVec.column(i);
            i++;
        }
    }
    return rVec;
}


#endif
