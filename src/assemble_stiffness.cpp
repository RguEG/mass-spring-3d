#include <assemble_stiffness.h>
#include <iostream>
typedef Eigen::Triplet<double> T;

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
    //Input:
//  q - generalized coordinates for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  E - the mx2 spring connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  l0 - the mx1 vector of undeformed spring lengths
//  k - the stiffness of each spring in the mass-spring system
//Output:
//  K - the 3nx3n sparse stiffness matrix which is the negative hessian of the potential energy function. 
    Eigen::Matrix66d H;
    Eigen::Vector3d  q0,q1;
    K.resize(q.size(),q.size());
    K.setZero();
    std::vector<T> tripletList;
    tripletList.reserve(9 * V.rows() * V.rows());

    for (int i = 0; i < E.rows(); i++) {
        q0[0] = q[3 * E(i, 0)];
        q0[1] = q[3 * E(i, 0) + 1];
        q0[2] = q[3 * E(i, 0) + 2];
        q1[0] = q[3 * E(i, 1)];
        q1[1] = q[3 * E(i, 1) + 1];
        q1[2] = q[3 * E(i, 1) + 2];
        d2V_spring_particle_particle_dq2(H, q0, q1, l0[i], k);
        for (int m = 0; m < 3;m++) {
            for (int n = 0; n < 3; n++) {
                tripletList.push_back(T(3 * E(i, 0) + m, 3 * E(i, 0) + n, H(m, n)));
                tripletList.push_back(T(3 * E(i, 0) + m, 3 * E(i, 1) + n, H(m, n+3.0)));
                tripletList.push_back(T(3 * E(i, 1) + m, 3 * E(i, 0) + n, H(m+3.0, n)));
                tripletList.push_back(T(3 * E(i, 1) + m, 3 * E(i, 1) + n, H(m+3.0, n+3.0)));
            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());
    };