Ryan Brossy, Melvin Laux
(Group 4)

Initialise Alignment

We compute the center of gravity for both the source and target meshes by iterating over all points in each mesh and summing their coordinates, and then dividing by the number of points. We then iterate over each mesh once again and find the point with the greatest distance from the center of gravity by continually updating a variable with the max of its current value and the distance for the current point. The scaling is then calculated by dividing the maximum distance of the source mesh by that of the target mesh. The translations are then assigned component-wise by taking the difference between the center points for x,y, and z and multiplying them by the scale.

Bad Pair Rejection

We loop through the correspondences vector, at each iteration extracting the appropriate vertices and normals according to the indices found in the ith pair of the vector. The angle between the normals is then calculated as the arccos of the dot product of the source and target divided by the product of their norms. An if statement is then used to check whether the distance between the corresponding vertices (norm of their difference) is less than 0.05, and if the calculated angle is less than 30 degrees (~0.5236 radians), and if the conditions are met the pair is pushed to the filtered correspondences vector.

Incremental Transformation

At each call to the function, we resize the b vector and J matrix according to the size of the correspondence vector. The centers of the corresponding points for both the source and target mesh are then found by taking the average of all of these points. The M and N matrices needed to contruct A are then built according to the provided formula (M = sum of the products of the barycentered source and transposed target vertices, N = sum of product of the barycentered target vertices with their transposes). The matrix A is then calculated by multiplying M with the inverse of N.

JacobiSVD from Eigen is then used to perform single value decomposition of the matrix A, and the rotation matrix constructed by taking the resultant U and V* matrices and multiplying them. The translation vector is then given by the center of the source mesh points minus the rotation matrix times the target mesh point center.

The derivatives of the component rotation matrices are then defined as below, by taking the partial derivatives of the respective rotation matrices and applying the small angle approximation:

  Rx' = [0.0, 0.0, 0.0,
  		 0.0, 0.0, -1.0, 
  		 0.0, 1.0, 0.0]

  Ry' = [0.0, 0.0, 1.0,
  		 0.0, 0.0, 0.0,
  		 -1.0, 0.0, 0.0]

  Rz' = [0.0, -1.0, 0.0,
  		 1.0, 0.0, 0.0,
  		 0.0, 0.0, 0.0]

 While the rotation matrices themselves become the identity matrix under the same approximation. We then construct the Jacobian for each point using the following definitions (per variable), obtained by evaluating the partial derivatives of the function for all 7 variables and setting scaling s to 0:

 Rotation: target_normal_transpose * Rx' * source_vertex (same for y and z)
 Translation: target_normal_transpose * x_unit_vector (same for y and z)
 Scaling: 2 * target_normal_transpose * rotation * source_vertex

 The b vector is given as follows: target_normal_transpose * (source_vertex - target_vertex).

 We initially had problems with divergence of the alignment by setting the b vector equal to normals*(rotation*source + translation - target)), as seemingly implied by the notes (b vector = f(x_t - 1)), and eventually after additional research found that the definition above caused the meshes to align properly - as this is the correct point-to-plane error estimator.

 This Jacobian and b vector are then passed to the solver to perform the least squares analysis and rigid alignment.