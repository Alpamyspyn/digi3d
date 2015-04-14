Ryan Brossy, Melvin Laux
(Group 4)

Splitting Long Edges

We iterate over all edges in the mesh, and for each select the 2 vertices involved. We set the desired edge length as half of the sum of the target length for each vertex. We then check if the edge is greater than 4/3rds of the desired edge length, and if so add a new vertex at the midpoint of the edge. We then calculate a normal and target length for this vertex by taking an average of these properties for the surrounding vertices, and call the mesh.split() function to finish the procedure. If an edge is split we set the finished variable is false so that the loop continues.


Collapsing Small Edges

We iterate over all (non-deleted) edges, and select the 2 vertices and 2 half edges for each.

We check if 4/5ths of the mean length is greater than the desired length, and if so then check if the edge is collapsible. To do this, we define booleans for each half edge indicating if they connect boundary with non-boundary vertices, and combine these with the mesh.is_collapse_ok() function. Using these values and the valences of each vertex, we then choose the appropriate (if any) half edge to collapse, setting the finished variable to false is a change is made. FInally, we call garbage collection to actually delete the collapsed edges.


Equalising Valences

The 2 vertices currently forming the edge are obtained using mesh.vertex, and the 2 vertices involved are flipping are obtained by taking the next_halfedge for the 2 halfedges currently forming the edge. The valances for these vertices are then stored.

The optimal valence values are assigned through ternary logic operators checking for boundary vertices. The deviation is then calculated as the sum of the squares of the valencies, and the edge flipped. Following this, the deviation from optimal valence is calculated again and the edge flipped back if this new value is not smaller than the old.


Tangential Smoothing

We calculate the laplacian for each vertex as in previous exercises. We then get the normal for the vertex and calculate the tangential component of the laplacian by subtracting the component of the the laplacian projected onto the normal vector. We obtain the normal component by taking the dot product of the laplacian with the normal vector, and multiplying by the normal. The update shift is stored, and we then loop through each vertex in the mesh and shift each using the update[] value.


Adaptive Remeshing

We store the mean curvature an gaussian curvature for each vertex, and compute the desired edge length using the provided formula (target length / H + sqrt(H^2 - K)), but set the desired length as target_length/K in the case of H^2 - K being negative to avoid NaN issues taking the square root of a negative number.

Uniform laplacian smoothing is applied over 5 iterations to the set of target lengths using the same procedure as for vertices but with scalars, adding half of the laplacian to the target length for each vertex per iteration.

In the last step we rescale the mean of the target lengths by calculating the current mean, and then dividing each length by this value and multiplying by the desired mean.


Results

We are aware, that the provided code does not produce the desired results. Many hours of debugging and double-checking the code unfortunately did not help to find our mistake(s). This could stem from the fact that our team consists of only two people with relatively little experience in C++ coding.