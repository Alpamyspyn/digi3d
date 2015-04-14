//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Mark Pauly, Mario Botsch, Balint Miklos, Hao Li, Sofien Bouaziz,
//   Yuliy Schwartzburg, Duygu Ceylan, Mario Deuss, Bailin Deng,
//   Alexandru Ichim, Anastasia Tkach
//
//   Copyright (C) 2007-2015 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne, Applied Geometry Group and Computer Graphics
//         Laboratory, ETH Zurich
//
//-----------------------------------------------------------------------------
//
//                                License
//
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License
//   as published by the Free Software Foundation; either version 2
//   of the License, or (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin Street, Fifth Floor,
//   Boston, MA  02110-1301, USA.
//
//=============================================================================

//== INCLUDES =================================================================

#include <surface_mesh/IO.h>
#include <surface_mesh/Surface_mesh.h>
#include <iostream>
#include <set>
#include <cmath>
#include <float.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

//== IMPLEMENTATION ===========================================================

using namespace surface_mesh;

typedef Surface_mesh            Mesh;
Mesh                            mesh;
Mesh::Vertex_property<Point>    update;

//desired edge length
Scalar                           target_length;

// Curvatures
Mesh::Vertex_property<Scalar>          vweight_, vunicurvature_, vcurvature_, vgausscurvature_, vtargetlength_, vnewtargetlength_;
Mesh::Edge_property<Scalar>            eweight_;
Mesh::Vertex_property<Normal>          vnormal_;

void  split_long_edges();
void  collapse_short_edges();
void  equalize_valences();
void  tangential_relaxation();

void  calc_weights();
void  calc_mean_curvature();
void  calc_uniform_mean_curvature();
void  calc_gauss_curvature();
void  calc_target_length();

void  set_normal(Mesh::Vertex v_, const Normal &normal_);
Normal get_normal(Mesh::Vertex v_);


int main(int argc, char **argv) {
	if (argc < 4) {
		std::cerr << "Usage: \n"
			<< argv[0] << " <edge-length>  <input_mesh>  <output_mesh>\n\n";
		exit(1);
	}

	target_length = atof(argv[1]);
	if (target_length <= 0) {
		std::cerr << "Error: edge-length needs to be positive." << std::endl;
		exit(1);
	}

	// read mesh
	if (!surface_mesh::read_mesh(mesh, argv[2])) {
		std::cerr << "Error: unable to parse the input mesh file." << std::endl;
		exit(1);
	}

	mesh.update_vertex_normals();
	mesh.update_face_normals();

	// initialize properties
	vcurvature_ = mesh.vertex_property<Scalar>("v:meancurvature", 0);
	vunicurvature_ = mesh.vertex_property<Scalar>("v:unimeancurvature", 0);
	vweight_ = mesh.vertex_property<Scalar>("v:weight", 0);
	eweight_ = mesh.edge_property<Scalar>("e:weight", 0);
	vgausscurvature_ = mesh.vertex_property<Scalar>("v:gausscurvature", 0);
    vtargetlength_ = mesh.vertex_property<Scalar>("v:length", 0);
    vnewtargetlength_ = mesh.vertex_property<Scalar>("v:newlength", 0);
	vnormal_ = mesh.vertex_property<Normal>("v:normal");
	update = mesh.vertex_property<Point>("v:update");


	calc_weights();
	calc_mean_curvature();
	calc_uniform_mean_curvature();
	calc_gauss_curvature();
    calc_target_length();


	// main remeshing loop
	for (int i = 0; i < 5; ++i) {
		split_long_edges();
		collapse_short_edges();
        equalize_valences();
        tangential_relaxation();
	}


	// write mesh
	if (!surface_mesh::write_mesh(mesh, argv[3])) {
		std::cerr << "Error: unable to write to the output mesh file." << std::endl;
		exit(1);
	}
}


void  split_long_edges() {
	Mesh::Edge_iterator     e_it, e_end(mesh.edges_end());
	Mesh::Vertex   v0, v1, v;
	bool            finished;
	int             i;

    std::printf("Splitting edges...\n");

    for (finished = false, i = 0; !finished && i < 100; ++i) {
		finished = true;
        std::printf("Iteration: %d; num_vertices: %d\n", i, mesh.n_vertices());

		for (e_it = mesh.edges_begin(); e_it != e_end; ++e_it) {
            // -----------------------------------------------------------
			// INSERT CODE:
			//  1) compute the desired length as the mean between the property vtargetlength_ of two vertices of the edge
			//  If the edge is longer than 4/3 * desired length
			//		2) add the midpoint to the mesh
			//		3) set the interpolated normal and interpolated vtargetlength_ property to the vertex
			//		4) split the edge with this vertex (use openMesh function split)
			// Leave the loop running until no splits are done (use the finished variable)
			// -----------------------------------------------------------

            // fetch vertices of the edge
            v0 = mesh.vertex(*e_it, 0);
            v1 = mesh.vertex(*e_it, 1);

            // calculate desired edge length
            float desired_length = 0.5f * (vtargetlength_[v0] + vtargetlength_[v1]);
            //printf("(target_v1, target_v2) = (%f, %f)\n", vtargetlength_[v0], vtargetlength_[v1]);
            //printf("desired: %f - length: %f \n",desired_length,mesh.edge_length(*e_it));
            if (mesh.edge_length(*e_it) > (4.0f / 3.0f) * desired_length)
            {
                //printf("desired: %f - length: %f \n",dlength,mesh.edge_length(*e_it));

                // calculate position of the new vertex
                Point mid = mesh.position(v0) + (mesh.position(v1) - mesh.position(v0)) / 2.0f;
                // add it to the mesh
                Mesh::Vertex v = mesh.add_vertex(mid);

                // interpolate normal andtargetlenth of the new vertex
                set_normal(v, (get_normal(v0) + get_normal(v1)) / 2.0f);//.normalize());
                vtargetlength_[v] = (vtargetlength_[v0] + vtargetlength_[v1]) / 2.0f;

                //printf("new target length: %f \n", get_normal(v));

                // split the edge
                mesh.split(*e_it, v);


                //printf("edge length: %f \n", mesh.edge_length(mesh.find_edge(v, v0)));
                //printf("edge length: %f \n", mesh.edge_length(mesh.find_edge(v, v1)));
                //printf("edge length: %f \n", mesh.edge_length(mesh.find_edge(v, v2)));
                //printf("edge length: %f \n", mesh.edge_length(mesh.find_edge(v, v3)));

                finished = false;
            }


        }

       // mesh.garbage_collection();
	}
}


void  collapse_short_edges() {
	Mesh::Edge_iterator     e_it, e_end(mesh.edges_end());
	Mesh::Vertex   v0, v1;
	Mesh::Halfedge  h01, h10;
	bool            finished, b0, b1;
	int             i;
	bool            hcol01, hcol10;

    std::printf("Collapsing short edges...\n");

    for (finished = false, i = 0; !finished && i < 100; ++i)
    {
        finished = true;

        std::printf("Iteration: %d; num_vertices: %d\n", i, mesh.n_vertices());

        for (e_it = mesh.edges_begin(); e_it != e_end; ++e_it)
        {
            if (!mesh.is_deleted(*e_it)) // might already be deleted
			{
                // -----------------------------------------------------------
				// INSERT CODE:
				//  1) Compute the desired length as the mean between the property vtargetlength_ of two vertices of the edge
				// If the edge is shorter than 4/5 of the desired length
				//		2) Check if halfedge connects a boundary vertex with a non-boundary vertex. If so, don't collapse. 
				//		3) Check if halfedges collapsible
				//		4) Select the halfedge to be collapsed if at least one halfedge can be collapsed
				//		5) Collapse the halfedge
				// Leave the loop running until no collapse has been done (use the finished variable)
                // -----------------------------------------------------------

                // fetch both halfedges of the current edge
                h01 = mesh.halfedge(*e_it,0);
                h10 = mesh.halfedge(*e_it,1);

                // fetch both vertices of the current edge
                v0 = mesh.from_vertex(h01);
                v1 = mesh.from_vertex(h10);

              //  printf("edge length: %f\n", mesh.edge_length(*e_it));

                // calculate desired edge length
                Scalar desired_length = (vtargetlength_[v0] + vtargetlength_[v1]) / 2.0f;

                // collapse edge if condition is met
                if(0.8f * desired_length > mesh.edge_length(*e_it))
                {
                    // check if halfedge is from boundary vertex to non-boundary vertex
                    b0 = mesh.is_boundary(v0) && !mesh.is_boundary(v1);
                    b1 = mesh.is_boundary(v1) && !mesh.is_boundary(v0);

                    // check if halfedge is collapsable
                    hcol01 = mesh.is_collapse_ok(h01) && !b0;
                    hcol10 = mesh.is_collapse_ok(h10) && !b1;

                    if(hcol01)
                    {
                        if(hcol10){
                            if(mesh.valence(v0) > mesh.valence(v1))
                            {
                                mesh.collapse(h10);
                                finished = false;
                            }
                            else
                            {
                                mesh.collapse(h01);
                                finished = false;
                            }
                        }
                        else
                        {
                            mesh.collapse(h01);
                            finished = false;
                        }
                    }
                    else if(hcol10)
                    {
                        mesh.collapse(h10);
                        finished = false;
                    }
                }
			}
		}
	}

	mesh.garbage_collection();

	if (i == 100) std::cerr << "collapse break\n";
}


void  equalize_valences() {
	Mesh::Edge_iterator e_it, e_end(mesh.edges_end());
	Mesh::Vertex v0, v1, v2, v3;
	Mesh::Halfedge h;
	int val0, val1, val2, val3;
	int val_opt0, val_opt1, val_opt2, val_opt3;
	int ve0, ve1, ve2, ve3, ve_before, ve_after;
	bool finished;
	int i;


    std::printf("Equalising valences...\n");

	// flip all edges
    for (finished = false, i = 0; !finished && i < 100; ++i)
    {
		finished = true;

        std::printf("Iteration: %d\n", i);

        for (e_it = mesh.edges_begin(); e_it != e_end; ++e_it)
        {
            if (!mesh.is_boundary(*e_it))
            {
                // -----------------------------------------------------------
				// INSERT CODE:
				//  1) Extract valences of the four vertices involved to an eventual flip.
				//  2) Compute the sum of the squared valence deviances before flip
				//  3) Compute the sum of the squared valence deviances after and eventual flip
				//  4) If valence deviance is decreased and flip is possible, flip the vertex
				// Leave the loop running until no collapse has been done (use the finished variable)
				// -----------------------------------------------------------

                if(!mesh.is_flip_ok(*e_it))
                    continue;

                h = mesh.halfedge(*e_it, 0);
                v0 = mesh.from_vertex(h);
                v1 = mesh.to_vertex(h);

                v2 = mesh.to_vertex(mesh.next_halfedge(h));

                h = mesh.opposite_halfedge(h);
                v3 = mesh.to_vertex(mesh.next_halfedge(h));

                val0 = mesh.valence(v0);
                val1 = mesh.valence(v1);
                val2 = mesh.valence(v2);
                val3 = mesh.valence(v3);

                val_opt0 = mesh.is_boundary(v0) ? 4 : 6;
                val_opt1 = mesh.is_boundary(v1) ? 4 : 6;
                val_opt2 = mesh.is_boundary(v2) ? 4 : 6;
                val_opt3 = mesh.is_boundary(v3) ? 4 : 6;

                ve0 = val0 - val_opt0;
                ve1 = val1 - val_opt1;
                ve2 = val2 - val_opt2;
                ve3 = val3 - val_opt3;

                ve_before = pow(ve0,2) + pow(ve1,2) + pow(ve2,2) + pow(ve3,2);

                ve0 = val0 - 1 - val_opt0;
                ve1 = val1 - 1 - val_opt1;
                ve2 = val2 + 1 - val_opt2;
                ve3 = val3 + 1 - val_opt3;

                ve_after = pow(ve0,2) + pow(ve1,2) + pow(ve2,2) + pow(ve3,2);

                if(ve_after < ve_before)
                {
                    mesh.flip(*e_it);
                    //printf("before: %d, after: %d ---> flip\n", ve_before, ve_after);
                    finished = false;
                }

                    //printf("before: %d, after: %d ---> no flip\n", ve_before, ve_after);


			}
		}
	}

    if (i == 100) std::cerr << "flip break\n";
}


void  tangential_relaxation() {
	Mesh::Vertex_iterator     v_it, v_end(mesh.vertices_end());
	Mesh::Vertex_around_vertex_circulator   vv_c, vv_end;
	int    valence;
	Point     u, n;
	Point     laplace;


    // smooth
    std::printf("Tangential relaxation...\n");
    for (int iters = 0; iters < 10; ++iters)
    {
        std::printf("Iteration: %d\n", iters);
        for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
        {
            if (!mesh.is_boundary(*v_it))
            {
                // -----------------------------------------------------------
                // INSERT CODE:
                //  1) Compute uniform laplacian curvature approximation vector
                //  2) Compute the tangential component of the laplacian vector and move the vertex
                //  3) Store smoothed vertex location in the update vertex property.
                //     (you don't have to use 1/2 attenuation in this case, it's fine without attenuation)
                // -----------------------------------------------------------

                laplace = Point(0.0);
                valence = 0;
                vv_c = mesh.vertices(*v_it);
                vv_end = vv_c;

                do{
                    laplace += (mesh.position(*vv_c) - mesh.position(*v_it));
                    ++valence;
                }while(++vv_c != vv_end);

                laplace /= valence;

                n = get_normal(*v_it);

                u = laplace - dot(laplace,n)*n;

                update[*v_it] = u;
            }
        }

        for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
            if (!mesh.is_boundary(*v_it))
                mesh.position(*v_it) += update[*v_it];
    }
}


void calc_target_length() {
	Mesh::Vertex_iterator v_it, v_end(mesh.vertices_end());
	Mesh::Vertex_around_vertex_circulator vv_c, vv_end;
	Scalar length;
	Scalar mean_length;
    Scalar H;
	Scalar K;

    std::printf("calculating target length...\n");


    // -----------------------------------------------------------
	// INSERT CODE:
	//  1) Get the maximal curvature at each vertex (use the precomputed mean (property: vcurvature_) and gaussian curvature (property: vgausscurvature_))
	//  2) Calculate the desired edge length as the target_length divided by the maximal curvature at each vertex, and assign it to the property vtargetlength_
	//  2) Smooth the maximal curvature uniformly, use the property vnewtargetlength_ to store the smoothed values intermediately
    //  3) Rescale the property vtargetlength_ such that its mean equals the user specified target_length
	// -----------------------------------------------------------

	// calculate desired length
    for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
    {
		length = 1.0;
        if (!mesh.is_boundary(*v_it))
        {
            // 1)
            H = vcurvature_[*v_it];
            K = vgausscurvature_[*v_it];
            length = target_length;

            if (pow(H, 2) - K >= 0.0)
                    length /= H + sqrt(pow(H ,2) - K);
            else
                    length /= H;
		}


        vtargetlength_[*v_it] = length;
	}

	// smooth desired length
    for (int i = 0; i < 5; i++)
    {

        // 2)
        for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
        {
            if(mesh.is_boundary(*v_it))
                continue;

            Scalar laplacian = 0.0;
            int n = 0;

            vv_c = mesh.vertices(*v_it);
            vv_end = vv_c;
            do
            {
                laplacian += vtargetlength_[*vv_c] - vtargetlength_[*v_it];
                ++n;
            }
            while(++vv_c != vv_end);

            laplacian /= n;
            vnewtargetlength_[*v_it] = vtargetlength_[*v_it] + 0.5 * laplacian;

        }


        for(v_it=mesh.vertices_begin(); v_it!=v_end;++v_it)
            vtargetlength_[*v_it] = vnewtargetlength_[*v_it];

	}


	// rescale desired length:
	// 3)

    mean_length = 0.0;

    for(v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
        mean_length += vtargetlength_[*v_it];

    mean_length /= mesh.n_vertices();

    for(v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
    {
        //vtargetlength_[*v_it] *= target_length / mean_length;
        vtargetlength_[*v_it] = target_length;
    }
}


void calc_weights() {
	Mesh::Vertex_iterator  v_it, v_end(mesh.vertices_end());
	Mesh::Edge_iterator          e_it, e_end(mesh.edges_end());
	Mesh::Face_around_vertex_circulator    vf_c, vf_end;
	Mesh::Vertex_around_face_circulator    fv_c;
	Mesh::Halfedge    h0, h1, h2;
	Mesh::Vertex      v0, v1;
	Point             p0, p1, p2, d0, d1;
	Scalar            w, area;



    for (e_it = mesh.edges_begin(); e_it != e_end; ++e_it)
    {
		w = 0.0;

		h0 = mesh.halfedge(*e_it, 0);
		v0 = mesh.to_vertex(h0);
		p0 = mesh.position(v0);

		h1 = mesh.halfedge(*e_it, 1);
		v1 = mesh.to_vertex(h1);
		p1 = mesh.position(v1);

		h2 = mesh.next_halfedge(h0);
		p2 = mesh.position(mesh.to_vertex(h2));
		d0 = normalize(p0 - p2);
		d1 = normalize(p1 - p2);
		w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, dot(d0, d1)))));

		h2 = mesh.next_halfedge(h1);
		p2 = mesh.position(mesh.to_vertex(h2));
		d0 = normalize(p0 - p2);
		d1 = normalize(p1 - p2);
		w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, dot(d0, d1)))));

		w = std::max(0.0f, w);
		eweight_[*e_it] = w * 0.5;
	}


    for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
    {
		area = 0.0;
		vf_c = mesh.faces(*v_it);
		vf_end = vf_c;

		do {
			fv_c = mesh.vertices(*vf_c);

			const Point& P = mesh.position(*fv_c);  ++fv_c;
			const Point& Q = mesh.position(*fv_c);  ++fv_c;
			const Point& R = mesh.position(*fv_c);

			area += norm(cross(Q - P, R - P)) * 0.5f * 0.3333f;

		} while (++vf_c != vf_end);

		vweight_[*v_it] = area;
	}
}


void calc_mean_curvature() {
	Mesh::Vertex_iterator        v_it, v_end(mesh.vertices_end());
	Mesh::Halfedge_around_vertex_circulator   vh_c, vh_end;
	Mesh::Vertex      neighbor_v;
	Mesh::Edge        e;
	Point             laplace(0.0);


	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
		Scalar curv = 0;

		if (!mesh.is_boundary(*v_it)) {
			laplace = Point(0.0);

			vh_c = mesh.halfedges(*v_it);
			vh_end = vh_c;

			do {
				e = mesh.edge(*vh_c);
				neighbor_v = mesh.to_vertex(*vh_c);
				laplace += eweight_[e] * (mesh.position(neighbor_v) - mesh.position(*v_it));

			} while (++vh_c != vh_end);

			laplace *= vweight_[*v_it];
			curv = 0.5 * norm(laplace);
		}

		vcurvature_[*v_it] = curv;
	}
}

void calc_uniform_mean_curvature() {
	Mesh::Vertex_iterator        v_it, v_end(mesh.vertices_end());
	Mesh::Vertex_around_vertex_circulator   vv_c, vv_end;
	Point             laplace(0.0);


	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
		Scalar curv = 0;

		if (!mesh.is_boundary(*v_it)) {
			laplace = Point(0.0);
			double n = 0;
			vv_c = mesh.vertices(*v_it);
			vv_end = vv_c;

			do {
				laplace += (mesh.position(*vv_c) - mesh.position(*v_it));
				++n;
			} while (++vv_c != vv_end);

			laplace /= n;

			curv = 0.5 * norm(laplace);
		}

		vunicurvature_[*v_it] = curv;
	}
}

void calc_gauss_curvature() {
	Mesh::Vertex_iterator        v_it, v_end(mesh.vertices_end());
	Mesh::Vertex_around_vertex_circulator   vv_c, vv_c2, vv_end;
	Point             d0, d1;
	Scalar            angles, cos_angle;
	Scalar            lb(-1.0f), ub(1.0f);


	// compute for all non-boundary vertices
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
		Scalar curv = 0;

		if (!mesh.is_boundary(*v_it)) {
			angles = 0.0;

			vv_c = mesh.vertices(*v_it);
			vv_end = vv_c;

			do {
				vv_c2 = vv_c;
				++vv_c2;
				d0 = normalize(mesh.position(*vv_c) - mesh.position(*v_it));
				d1 = normalize(mesh.position(*vv_c2) - mesh.position(*v_it));
				cos_angle = std::max(lb, std::min(ub, dot(d0, d1)));
				angles += acos(cos_angle);
			} while (++vv_c != vv_end);

			curv = (2 * M_PI - angles) * 2.0f * vweight_[*v_it];
		}

		vgausscurvature_[*v_it] = curv;
	}

}

void set_normal(Mesh::Vertex v_, const Normal &normal_) {
	vnormal_[v_] = normal_;
}

Normal get_normal(Mesh::Vertex v_) {
	return vnormal_[v_];
}
