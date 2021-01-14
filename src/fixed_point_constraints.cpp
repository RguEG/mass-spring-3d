#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    //Input:
//  q_size - total number of scalar generalized coordinaes (3 times number of vertices in the mesh) 3n
//  indices - indices (row ids in V) for fixed vertices indices.size() = l
//Output:
//  P - mxn sparse matrix which projects out fixed vertices 3(n-l) x 3n

    P.resize(q_size-3*indices.size(),q_size);
	P.setZero();
    int count = 0;


	for (int i = 0; i < q_size - 3 * indices.size(); i++) {
		if (count < indices.size() && i == 3.0 * indices[count] - 3.0 * count) { count++; }
		P.insert(i, i + 3.0 * count) = 1.0;
	}

}