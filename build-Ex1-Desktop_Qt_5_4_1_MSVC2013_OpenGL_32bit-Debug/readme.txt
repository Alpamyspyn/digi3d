Ryan Brossy, Melvin Laux, Solomon Kamugasa

Exercise 1.2
We implemented the calc_valences function by iterating over all vertices in the mesh using Mesh::Vertex_iterator, 
and at each vertex used the Mesh::Vertex_around_vertex_circulator to iterate through the neighbourhood, counting the vertices.
No changes were made to the valence viewer class.

Exercise 1.3



We scaled the valence values of all vertices (excluding the maximum and minimum to get a better distribution) and then 
mapped the colours in the following way:

blue: 	min(max(0.75 - normalized_valence, 0), 1);
red: 	min(max(normalized_valence - 0.25, 0), 1);
green:  min(max(abs(normalized_valence - 0.5) - 1, 0), 1);

With this continuous mapping we receive blue vertices for low values and yellow values for high values.