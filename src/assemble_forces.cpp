#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) {
//Input:
//  q - generalized coordinates for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  E - the mx2 spring connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  l0 - the mx1 vector of undeformed spring lengths
//  m - the mass of each particle in the mass-spring system
//  k - the stiffness of each spring in the mass-spring system
//Output:
//  f - the vector 3xn vector of forces acting on each node of the mass-spring system
    Eigen::Vector6d flocal;
    Eigen::Vector3d q0,q1;
    double q_0, q_1, q_2, q_3, q_4, q_5;
    f.resize(q.size());
    f.setZero();
    for (int i = 0; i < E.rows();i++) {
        q0[0] = q[3.0 * E(i, 0)];
        q0[1] = q[3.0 * E(i, 0) + 1];
        q0[2] = q[3.0 * E(i, 0) + 2];
        q1[0] = q[3.0 * E(i, 1)];
        q1[1] = q[3.0 * E(i, 1) + 1];
        q1[2] = q[3.0 * E(i, 1) + 2];
        dV_spring_particle_particle_dq(flocal, q0, q1, l0[i], k);
        f[3.0 * E(i, 0)] += flocal[0];
        f[3.0 * E(i, 0) + 1] += flocal[1];
        f[3.0 * E(i, 0) + 2] += flocal[2];
        f[3.0 * E(i, 1)] += flocal[3];
        f[3.0 * E(i, 1) + 1] += flocal[4];
        f[3.0 * E(i, 1) + 2] += flocal[5];
    }
    f = -1.0 * f;
    }