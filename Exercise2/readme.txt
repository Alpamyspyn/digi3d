Ryan Brossy, Melvin Laux
(Group 4)

Exercise 2.1
We stored the indices returned by the knn_search() function for point i and (n_neighbours + 1) in a temporary vector, and summed all the elements of this vector. This was used to get p_bar. The matrix D was then constructed by taking the difference of each of the neighbourhood points and p_bar and storing them in the columns. M was then constructed by taking the product of D and its transpose, and the normals obtained by passing this matrix to the smallest_eigenvector() function.


Exercise 2.2

Using Hoppe's distance gave a sphere with a highly distorted surface compared to what would be desired given the point data.


Exercise 2.3

Bounding box dimensions are found through the use of the maxCoeff() and minCoeff() functions on rows of the points_ matrix, and the diagonal given by the L2 norm of the difference between the maximum and minimum obtained.

We then define the centers by storing the coordinates of each available point as well as the offset centre obtained by adding 0.001*BoxDiagonal multiplied by the corresponding normal vector. These are then used to setup the matrix of coefficients M through the use of the kernel() function. 

The distance constraint vector is obtained by setting the value for all centres to 0 and for the offset centres to 0.001*BoxDiagonal. 

Finally, the coefficient matrix and distance constraint vector are passed to the linear solver.

Exercise 2.4

In the case of using 1000 neighbours to compute the normals (sphere_1000Neighbours.png) we see that there is a severe distortion in the final sphere, and all normals are facing the in the same direction.

For the minimal number of neighbours (sphere_3Neighbours.png) the sphere still renders normally, We can therefore say that for the sphere this is a suitable value.

In changing the offset value, even multiplying the box diagonal by 0.0000001 (sphere_0000001.png) still resulted in a perfectly rendered sphere.

Using an offset of 5 times the box diagonal, we also get a normally rendered sphere. This indicated that at least in the case of a sphere the final result is not very depedant on the offset. Other objects took extremely long times to render and we were unable to test enough cases in time.
