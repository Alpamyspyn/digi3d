Exercise 2.1
We stored the indices returned by the knn_search() function for point i and (n_neighbours + 1) in a temporary vector, and summed all the elements of this vector. This was used to get p_bar. The matrix D was then constructed by taking the difference of each of the neighbourhood points and p_bar and storing them in the columns. M was then constructed by taking the product of D and its transpose, and the normals obtained by passing this matrix to the smallest_eigenvector() function.


Exercise 2.2

<observations>


Exercise 2.3

Bounding box dimensions are found through the use of the maxCoeff() and minCoeff() functions on rows of the points_ matrix, and the diagonal given by the L2 norm of the difference between the maximum and minimum obtained.

We then define the matrix of centers by ... and use these to setup the matrix of coefficients M.

The distance constraint vector and coefficient matrix are then passed to the linear solver.

Exercise 2.4

