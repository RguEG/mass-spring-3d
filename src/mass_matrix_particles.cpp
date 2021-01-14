#include <mass_matrix_particles.h>

void mass_matrix_particles(Eigen::SparseMatrixd& M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {
    //Input:
//  q - generalized coordiantes for the mass-spring system
//  mass - the mass of each particle in the mass-spring system.
//Output:
//  M - sparse mass matrix for mass-spring system

    M.resize(q.size(),q.size());
    M.setZero();
    for (int i = 0; i < M.rows(); i++) {
                M.insert(i, i) = mass;
    }
}
