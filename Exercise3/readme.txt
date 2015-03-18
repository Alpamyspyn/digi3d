Ryan Brossy, Melvin Laux
(Group 4)

Exercise 3.1 Uniform Laplace

Use a vertex iterator and vertex around vertex circulator to go through the neighbourhood of each vertex. At each iteration we sum the vectors given by the difference between the vertex and its neighbours, then average this resultant vector and store half its norm.

Exercise 3.2 Laplace-Beltrami

We use a vertex iterator and halfedge around vertex circulator to iterate through the set of half edges for each vertex in the mesh. For each halfedge, of the current vertex we determine the neighbour vertex it is pointing at. We then use the precalculated (edge)weights to produce the sum of the Laplaceian and eventually multiply this sum with the precalculated (vertex)weights, where the (edge)weights represent the the sum of the angles of the sum within the Laplacian and the (vertex)weight represent the area around the current vertex. Finally, we store the half of the resulting vector as the approximation of the Laplace-Beltrami curvature.

Exercise 3.3 Triangle Shape

We loop over all triangles in the mesh using a face iterator, and within this loop we store the positions of the 3 vertices making up the triangle. These values are used to construct the 3 edge vectors, and the minimum edge length retrieved. For negative values of the denominator (a x b), the triangle shape is assigned FLT_MAX, else we store the value of half the product over the edge lengths divided by the denominator (value stored in circum_radius_sq) divided by the minimum edge length.



Exercise 3.4 Gaussian Curvature
