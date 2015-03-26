Ryan Brossy, Melvin Laux
(Group 4)

Exercise 3.1 Uniform Laplacian Smoothing

We calculate the uniform Laplacian operator as in the previous exercise. In the first iteration we initialise v_old_pos_ to the initial positions. After each iteration we store the new vertex positions in v_new_pos, and after all calculations are complete apply the shift given by the uniform Laplacian smoothing formula.

Exercise 3.2 Laplace-Beltrami Smoothing

Implemented in an identical manner to 3.1 above, but using the Laplace-Beltrami operator.

Exercise 3.3 Feature Enhancement

We first call either the uniform or Laplace-Beltrami smoothing function, and then shift each vertex in the smoothed mesh by the difference between the initial and smoothed meshes, multiplied by the enhancement factor.

Enhanced Images:

Image1.png:

Applied feature enhanced with 10 uniform smoothing iterations and enhancement coefficient of 2.0, 5 times.

We observe the mesh forming jagged edges - this is due to the fact that after each iteration the gap between the initial and smoothed mesh is larger, compounding the effect of enhancement and forming points at the extremes. The effect is more noticeable with a higher enhancement coefficient.  From a signal processing point of view this can be seen as applying a high pass filter multiple times.

image2.png:

To the already processed image above we applied a further 3 Laplace-Beltrami smoothing iterations with the same settings.

We regain much of the original appearance of the model, however with a less refined overall look. Smoothing is equivalent to a low pass filter, and so we are cutting out the more exaggerated features and begin to return to the original model, but some detail is lost in the process.


image3.png:

We apply 2 uniform smoothings followed by 3 uniform Laplacian enhancements, with 10 iterations and an enhancement coefficient of 3.0

We observe sharper/more exaggerated edges and a "squaring" of the features. By first smoothing we filter out the fine details with a low pass filter, and then only the largest features remain to be enhanced by the high pass filter.

image4.png

Compared to image3, we apply an additional 5 smoothings and 6 enhancements, 2 smoothings, and a final enhancement with the same settings.

We obtain an even more exaggerated analogue of image 3, with only the strongest features being grossly enhanced.