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
//   Copyright (C) 2015 by Bailin Deng <bldeng@gmail.com>
//   Copyright (C) 2006-2015 by Computer Graphics and Geometry Laboratory, EPFL
//                              Computer Graphics Laboratory, ETH Zurich,
//                              and Computer Graphics Group, RWTH Aachen
//
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
//=============================================================================
//
//  CLASS QualityViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "QualityViewer.hh"
#include <vector>
#include <float.h>
#include <math.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif


//== IMPLEMENTATION ========================================================== 


QualityViewer::
QualityViewer(const char* _title, int _width, int _height)
: MeshViewer(_title, _width, _height)
{ 
	add_draw_mode("Uniform Mean Curvature");
	add_draw_mode("Mean Curvature");
	add_draw_mode("Gaussian Curvature");
	add_draw_mode("Triangle Shape");
	add_draw_mode("Reflection Lines");

	init();
}


//-----------------------------------------------------------------------------


QualityViewer::
~QualityViewer()
{
	if (glIsTexture(textureID_))  
		glDeleteTextures( 1, &textureID_);
}

//-----------------------------------------------------------------------------


void
QualityViewer::
init()
{
	// base class first
	MeshViewer::init();


	// generate checkerboard-like image
	GLubyte tex[256*256*3], *tp=tex;
    for (int x=0; x<256; ++x)
		for (int y=0; y<256; ++y)
			if (((x+2)/4 % 10) == 0 || ((y+2)/4 % 10) == 0)
			{
				*(tp++) = 0;
				*(tp++) = 0;
				*(tp++) = 0;
			}
			else
			{
				*(tp++) = 255;
				*(tp++) = 255;
				*(tp++) = 255;
			}


    // generate texture
    if (!glIsTexture(textureID_))
        glGenTextures(1, &textureID_);
    glBindTexture(GL_TEXTURE_2D, textureID_);


    // copy texture to GL
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, 3, 256, 256,
        0, GL_RGB, GL_UNSIGNED_BYTE, tex);
}


//-----------------------------------------------------------------------------

void
QualityViewer::
init_properties()
{
    vcurvature_ = mesh_.vertex_property<Scalar>("v:meancurvature", 0);
    vunicurvature_ = mesh_.vertex_property<Scalar>("v:unimeancurvature", 0);
    vweight_ = mesh_.vertex_property<Scalar>("v:weight", 0);
    eweight_ = mesh_.edge_property<Scalar>("e:weight", 0);
    tshape_ = mesh_.face_property<Scalar>("f:tshape", 0);
    vgausscurvature_ = mesh_.vertex_property<Scalar>("v:gausscurvature", 0);
}


//-----------------------------------------------------------------------------


bool
QualityViewer::
open_mesh(const char* _filename)
{
	// load mesh
	if (MeshViewer::open_mesh(_filename))
	{
        // initialize the properties we use
        init_properties();

		// compute curvature stuff
		calc_weights();
		calc_mean_curvature();
		calc_uniform_mean_curvature();
		calc_gauss_curvature();
		calc_triangle_quality();
		face_color_coding();

		glutPostRedisplay();
		return true;
	}
	return false;
}


//-----------------------------------------------------------------------------



void
QualityViewer::
calc_weights()
{
    calc_edge_weights();
    calc_vertex_weights();
}


//-----------------------------------------------------------------------------


void
QualityViewer::
calc_edge_weights()
{
    Mesh::Edge_iterator          e_it, e_end(mesh_.edges_end());
    Mesh::Halfedge    h0, h1, h2;
    Mesh::Vertex      v0, v1;
    Point             p0, p1, p2, d0, d1;
    Scalar            w;



    for (e_it=mesh_.edges_begin(); e_it!=e_end; ++e_it)
    {
        w  = 0.0;

        h0 = mesh_.halfedge(*e_it, 0);
        v0 = mesh_.to_vertex(h0);
        p0 = mesh_.position(v0);

        h1 = mesh_.halfedge(*e_it, 1);
        v1 = mesh_.to_vertex(h1);
        p1 = mesh_.position(v1);

        h2 = mesh_.next_halfedge(h0);
        p2 = mesh_.position(mesh_.to_vertex(h2));
        d0 = normalize(p0 - p2);
        d1 = normalize(p1 - p2);
        w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, dot(d0, d1)))));

        h2 = mesh_.next_halfedge(h1);
        p2 = mesh_.position(mesh_.to_vertex(h2));
        d0 = normalize(p0 - p2);
        d1 = normalize(p1 - p2);
        w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, dot(d0, d1)))));

        w = std::max(0.0f, w);
        eweight_[*e_it] = w * 0.5;
    }
}



//-----------------------------------------------------------------------------


void
QualityViewer::
calc_vertex_weights()
{
    Mesh::Vertex_iterator  v_it, v_end(mesh_.vertices_end());
    Mesh::Face_around_vertex_circulator    vf_c, vf_end;
    Mesh::Vertex_around_face_circulator    fv_c;
    Scalar            area;

    for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
    {
        area = 0.0;
        vf_c = mesh_.faces(*v_it);

        // Check for isolated vertices
        if(!vf_c){
            continue;
        }

        vf_end = vf_c;

        do{
            fv_c = mesh_.vertices(*vf_c);

            const Point& P = mesh_.position(*fv_c);  ++fv_c;
            const Point& Q = mesh_.position(*fv_c);  ++fv_c;
            const Point& R = mesh_.position(*fv_c);

            area += norm(cross(Q-P, R-P)) * 0.5f * 0.3333f;

        }while(++vf_c != vf_end);

        vweight_[*v_it] = 0.5 / area;
    }
}


//-----------------------------------------------------------------------------


void 
QualityViewer::
calc_mean_curvature()
{
    Mesh::Vertex_iterator        v_it, v_end(mesh_.vertices_end());
    Mesh::Halfedge_around_vertex_circulator   vh_c, vh_end;
    Mesh::Vertex      neighbor_v;
    Mesh::Edge        e;
    Point             laplace(0.0);


	for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
	{
        Scalar curv = 0;

        if (!mesh_.is_boundary(*v_it))
		{
            laplace = Point(0.0);

            vh_c = mesh_.halfedges(*v_it);
            vh_end = vh_c;

            do{
                e = mesh_.edge(*vh_c);
                neighbor_v = mesh_.to_vertex(*vh_c);
                laplace += eweight_[e] * (mesh_.position(neighbor_v) - mesh_.position(*v_it));

            }while(++vh_c != vh_end);

            laplace *= vweight_[*v_it];
            curv = 0.5 * norm(laplace);
		}

        vcurvature_[*v_it] = curv;
	}
}

void 
QualityViewer::
calc_uniform_mean_curvature()
{
    Mesh::Vertex_iterator        v_it, v_end(mesh_.vertices_end());
    Mesh::Vertex_around_vertex_circulator   vv_c, vv_end;
    Point             laplace(0.0);


	for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
	{
        Scalar curv = 0;

        if (!mesh_.is_boundary(*v_it))
		{
            laplace = Point(0.0);
			double n = 0;
            vv_c = mesh_.vertices(*v_it);        
            vv_end = vv_c;

            do{
                laplace += (mesh_.position(*vv_c) - mesh_.position(*v_it));
                ++n;
            }while(++vv_c != vv_end);

            laplace /= n;

            curv = 0.5 * norm(laplace);
		}

        vunicurvature_[*v_it] = curv;
	}
}

void 
QualityViewer::
calc_gauss_curvature()
{
    Mesh::Vertex_iterator        v_it, v_end(mesh_.vertices_end());
    Mesh::Vertex_around_vertex_circulator   vv_c, vv_c2, vv_end;
    Point             d0, d1;
    Scalar            angles, cos_angle;
    Scalar            lb(-1.0f), ub(1.0f);


	// compute for all non-boundary vertices
	for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
	{
        Scalar curv = 0;

        if (!mesh_.is_boundary(*v_it))
		{
			angles = 0.0;

            vv_c = mesh_.vertices(*v_it);
            vv_end = vv_c;

            do{
                vv_c2 = vv_c;
                ++ vv_c2;
                d0 = normalize(mesh_.position(*vv_c) - mesh_.position(*v_it));
                d1 = normalize(mesh_.position(*vv_c2) - mesh_.position(*v_it));
                cos_angle = std::max(lb, std::min(ub, dot(d0, d1)));
                angles += acos(cos_angle);
            }while(++vv_c != vv_end);

            curv = (2 * M_PI - angles) * 2.0f * vweight_[*v_it];
		}

        vgausscurvature_[*v_it] = curv;
	}

}

//-----------------------------------------------------------------------------


void 
QualityViewer::
calc_triangle_quality()
{
    Mesh::Face_iterator				f_it, f_end(mesh_.faces_end());
    Mesh::Vertex_around_face_circulator	fv_c;
    Point				v0,v1,v2;
    Point				v0v1,v0v2,v1v2;
    Scalar				denom, circum_radius_sq, min_length_sq;

	for (f_it=mesh_.faces_begin(); f_it!=f_end; ++f_it)
	{
        fv_c = mesh_.vertices(*f_it);

        v0 = mesh_.position(*fv_c); ++fv_c;
        v1 = mesh_.position(*fv_c); ++fv_c;
        v2 = mesh_.position(*fv_c);

		v0v1 = v1-v0;
		v0v2 = v2-v0;
		v1v2 = v2-v1;

        denom = 4.0*sqrnorm(cross(v0v1,v0v2));
		if (fabs(denom) < FLT_MIN)  circum_radius_sq = FLT_MAX;
        else circum_radius_sq = sqrnorm(v0v1) * sqrnorm(v0v2) * sqrnorm(v1v2) / denom;
		
        min_length_sq = sqrnorm(v0v1);
        if (sqrnorm(v0v2) < min_length_sq) min_length_sq = sqrnorm(v0v2);
        if (sqrnorm(v1v2) < min_length_sq) min_length_sq = sqrnorm(v1v2);
		
        tshape_[*f_it] = sqrt(circum_radius_sq/min_length_sq);
	}
}

//-----------------------------------------------------------------------------

void 
QualityViewer::
face_color_coding()
{
	face_colors_.clear();
	face_colors_.reserve(mesh_.n_faces()*3);

    Scalar min_shape = 0.6f, max_shape = 2.0f;

    // map values to colors
    Mesh::Face_iterator f_it, f_end(mesh_.faces_end());
    for (f_it = mesh_.faces_begin(); f_it!=f_end; ++f_it)
	{
        Color col = value_to_color(tshape_[*f_it], min_shape, max_shape);
        face_colors_.push_back(col[0]);
        face_colors_.push_back(col[1]);
        face_colors_.push_back(col[2]);
	}
}



//-----------------------------------------------------------------------------


void 
QualityViewer::
color_coding(Mesh::Vertex_property<Scalar> prop)
{	
    // Get the value array
    std::vector<Scalar> values = prop.vector();

    // discard upper and lower 5%
	unsigned int n = values.size()-1;
	unsigned int i = n / 20;
	std::sort(values.begin(), values.end());
    Scalar min_value = values[i], max_value = values[n-1-i];

    // map values to colors
    Mesh::Vertex_iterator v_it, v_end(mesh_.vertices_end());
	for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
	{
        set_color(*v_it, value_to_color(prop[*v_it], min_value, max_value));
	}
}


Color
QualityViewer::
value_to_color(Scalar value, Scalar min_value, Scalar max_value) {
    Scalar v0, v1, v2, v3, v4;
    v0 = min_value + 0.0/4.0 * (max_value - min_value);
    v1 = min_value + 1.0/4.0 * (max_value - min_value);
    v2 = min_value + 2.0/4.0 * (max_value - min_value);
    v3 = min_value + 3.0/4.0 * (max_value - min_value);
    v4 = min_value + 4.0/4.0 * (max_value - min_value);

    Color col(1.0);

    if (value < v0) col = Color(0, 0, 1);
    else if (value > v4) col = Color(1, 0, 0);

	else if (value <= v2) 
	{
		if (value <= v1) // [v0, v1]
		{
            Scalar u =  (value - v0) / (v1 - v0);
            col = Color(0, u, 1);
		}      
		else // ]v1, v2]
		{
            Scalar u = (value - v1) / (v2 - v1);
            col = Color(0, 1, 1-u);
		}
	}
	else 
	{
		if (value <= v3) // ]v2, v3]
		{
            Scalar u = (value - v2) / (v3 - v2);
            col = Color(u, 1, 0);
		}
		else // ]v3, v4]
		{
            Scalar u = (value - v3) / (v4 - v3);
            col = Color(1, 1-u, 0);
		}
	}

	return col;
}


//-----------------------------------------------------------------------------


void 
QualityViewer::
draw(const std::string& _draw_mode)
{

	if (indices_.empty())
	{
		MeshViewer::draw(_draw_mode);
		return;
	}

	if (_draw_mode == "Mean Curvature") color_coding(vcurvature_);
	if (_draw_mode == "Gaussian Curvature") color_coding(vgausscurvature_);
	if (_draw_mode == "Uniform Mean Curvature") color_coding(vunicurvature_);

	if (_draw_mode == "Mean Curvature" || _draw_mode == "Gaussian Curvature" || _draw_mode == "Uniform Mean Curvature")
	{

		glDisable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);
        GL::glVertexPointer(vertex_positions_pointer());
        GL::glNormalPointer(vertex_normals_pointer());
        GL::glColorPointer(vertex_colors_pointer());
		
		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);

	}

	if (_draw_mode == "Triangle Shape")
	{

		glDisable(GL_LIGHTING);
		glShadeModel(GL_FLAT);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
        GL::glVertexPointer(vertex_positions_pointer());
        GL::glNormalPointer(vertex_normals_pointer());


		glDepthRange(0.01, 1.0);
		glBegin(GL_TRIANGLES);
		for (unsigned i=0; i<indices_.size(); i++)
		{
            if (i%3==0) glColor3f(face_colors_[i], face_colors_[i+1], face_colors_[i+2]);
			glArrayElement(indices_[i]);
		}
		glEnd();


		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);

		glColor3f(0.3, 0.3, 0.3);

		glEnableClientState(GL_VERTEX_ARRAY);
        GL::glVertexPointer(vertex_positions_pointer());

		glDrawBuffer(GL_BACK);
		glDepthRange(0.0, 1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDepthFunc(GL_LEQUAL);
		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDepthFunc(GL_LESS);
	}

	else if (_draw_mode == "Reflection Lines")
	{
		glTexGeni( GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP );
		glTexGeni( GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP );
		glEnable( GL_TEXTURE_GEN_S );
		glEnable( GL_TEXTURE_GEN_T );
		glEnable( GL_TEXTURE_2D );    
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
        GL::glVertexPointer(vertex_positions_pointer());
        GL::glNormalPointer(vertex_normals_pointer());

		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);

		glDisable( GL_TEXTURE_GEN_S );
		glDisable( GL_TEXTURE_GEN_T );
		glDisable( GL_TEXTURE_2D );
	}



	else MeshViewer::draw(_draw_mode);
}


//=============================================================================