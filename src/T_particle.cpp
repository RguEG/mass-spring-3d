#include <T_particle.h>
#include <iostream>

void T_particle(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, double mass) {
    T = 0.0;
    typedef Eigen::SparseMatrix<double> SparseMatrixType;
    SparseMatrixType M;

    M.resize(qdot.size(), qdot.size());
    M.setZero();
    for (int i = 0; i < M.rows(); i++) {
        M.insert(i, i) = mass;
    }

    T = 0.5 * qdot.transpose() * M * qdot;

}