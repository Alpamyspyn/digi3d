//=============================================================================
//
//   Code framework for the lecture
//
//   "Surface Representation and Geometric Modeling"
//
//   Mark Pauly, Mario Botsch, Balint Miklos, and Hao Li
//
//   Copyright (C) 2007 by  Applied Geometry Group and
//							Computer Graphics Laboratory, ETH Zurich
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
//  CLASS RegistrationViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include <surface_mesh/IO.h>
#include <surface_mesh/Surface_mesh.h>
#include "RegistrationViewer.hh"
#include "Transformation.hh"
#include "Fitting.hh"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>


//== IMPLEMENTATION ==========================================================


RegistrationViewer::
RegistrationViewer(const char* _title, int _width, int _height)
    : GlutExaminer(_title, _width, _height)
{
    clear_draw_modes();

    src_index_ = 0;
    numProcessed_ = 0;

    mode_ = VIEW;
}


//-----------------------------------------------------------------------------


RegistrationViewer::
~RegistrationViewer()
{
}

//-----------------------------------------------------------------------------


/// set output filename
void
RegistrationViewer::
set_output( const std::string & filename )
{
    outputFilename_ = filename;
}

//-----------------------------------------------------------------------------

bool
RegistrationViewer::
open_meshes(const std::vector<std::string> & _filenames,
            const std::vector<size_t> &indices_template)
{
    indices_template_ = indices_template;

    for(int i = 0; i < (int) _filenames.size(); i++)
    {
        Mesh mesh;
        // load mesh
        if (surface_mesh::read_mesh(mesh, _filenames[i]))
        {

            // compute face & vertex normals
            mesh.update_face_normals();
            mesh.update_vertex_normals();


            // update face indices for faster rendering
            update_face_indices();


            // info
            std::cerr << mesh.n_vertices() << " vertices, "
                      << mesh.n_faces()    << " faces\n";

            meshes_.push_back( mesh );

            transformations_.push_back( Transformation() );
        }
        else
        {
            std::cerr << "Problems reading the meshes\n";
            return false;
        }
    }

    return true;
}

bool RegistrationViewer::initialize() {

    initialize_alignment (meshes_[0], meshes_[1], transformations_[1]);

    Vec3f bbMin, bbMax;
    bbMin = bbMax = meshes_[0].position(*meshes_[0].vertices_begin ());
    for (Mesh::Vertex_iterator v_it = meshes_[0].vertices_begin ();
         v_it != meshes_[0].vertices_end (); ++v_it)
    {
        bbMin.minimize(meshes_[0].position(*v_it));
        bbMax.maximize(meshes_[0].position(*v_it));
    }

    set_scene((bbMin + bbMax)*0.5, 0.5* norm (bbMin - bbMax));


    // update face indices for faster rendering
    update_face_indices();
    glutPostRedisplay();

    numProcessed_ = std::min( 2, int(meshes_.size()) );
    src_index_ = 1; //std::max( 0, numProcessed_-1 );

    d_previous = Eigen::VectorXd::Zero (3 * meshes_[0].n_vertices ());
    return true;

}


//-----------------------------------------------------------------------------


void
RegistrationViewer::
update_face_indices()
{
    indices_.clear();

    for(int i = 0; i < (int) meshes_.size(); i++)
    {
        indices_.push_back( std::vector<unsigned int>() );

        Mesh::Face_iterator        f_it(meshes_[i].faces_begin()),
                f_end(meshes_[i].faces_end());
        Mesh::Vertex_around_face_circulator  fv_c;

        indices_[i].clear();
        indices_[i].reserve(meshes_[i].n_faces()*3);

        for (; f_it!=f_end; ++f_it)
        {
            fv_c=meshes_[i].vertices(*f_it);
            indices_[i].push_back((*fv_c).idx());

            ++ fv_c;
            indices_[i].push_back((*fv_c).idx());

            ++ fv_c;
            indices_[i].push_back((*fv_c).idx());
        }
    }
}


//-----------------------------------------------------------------------------


void
RegistrationViewer::
draw(const std::string& _draw_mode)
{

    if (indices_.empty())
    {
        GlutExaminer::draw(_draw_mode);
        return;
    }

    // display scans
    for(int i = 0; i < numProcessed_; i++)
    {
        if( i == src_index_ )
            draw(i, Vec3f(0.1,0.5,0.1) );
        else
            draw(i, Vec3f(0.5,0.5,0.5) );
    }
}



//-----------------------------------------------------------------------------


void
RegistrationViewer::
draw(int index, const Vec3f & color)
{
    /// TODO - have a drawWireframe (target) and a drawSmooth (source)


    glPushMatrix();
    // apply transformation matrix of scan
    transformations_[index].apply_gl();

    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);
    //glEnable(GL_RESCALE_NORMAL);
    glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);
    glColor3f(color[0], color[1], color[2]);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    GL::glVertexPointer(&meshes_[index].points()[0]);
    GL::glNormalPointer(&meshes_[index].vertex_property<Normal>("v:normal").vector ()[0]);

    glDrawElements(GL_TRIANGLES, indices_[index].size(), GL_UNSIGNED_INT, &(indices_[index][0]));

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_COLOR_MATERIAL);

    glPopMatrix();
}



//-----------------------------------------------------------------------------

void
RegistrationViewer::
keyboard(int key, int x, int y)
{

    switch (key)
    {

    case 'r':
    {
        std::cout << std::endl << "RIGID REGISTRATION" << std::endl;
        rigid_registration (meshes_[0], meshes_[1],
                indices_template_,
                transformations_[0], transformations_[1]);
        glutPostRedisplay();
        break;
    }
    case 'n':
    {
        std::cout << std::endl << "NONRIGID REGISTRATION" << std::endl;
        Eigen::VectorXd d;
        nonrigid_registration (meshes_[0], meshes_[1],
                indices_template_,
                transformations_[0], transformations_[1], d_previous, d);
        d_previous = d;
        glutPostRedisplay();
        break;
    }
    case 's':
    {
        // Save the fitted mesh
        surface_mesh::write_obj (meshes_[0], "fitted_template.obj");

        // Transform the mesh
        Mesh mesh_src_tr = meshes_[0];
        for (Mesh::Vertex_iterator v_it = mesh_src_tr.vertices_begin ();
             v_it != mesh_src_tr.vertices_end (); ++v_it)
        {
            const Vec3 &p = mesh_src_tr.position (*v_it);
            mesh_src_tr.position (*v_it) = transformations_[1].inverse ().transformPoint (p);

            const Vec3 &n = mesh_src_tr.get_vertex_property<Normal> ("v:normal")[*v_it];
            mesh_src_tr.vertex_property<Normal> ("v:normal")[*v_it] = transformations_[1].inverse ().transformNormal (n);
        }
        surface_mesh::write_obj (mesh_src_tr, "fitted_template_tr.obj");


        // Transform the point cloud
        Mesh mesh_tgt_tr = meshes_[1];
        for (Mesh::Vertex_iterator v_it = mesh_tgt_tr.vertices_begin ();
             v_it != mesh_tgt_tr.vertices_end (); ++v_it)
        {
            const Vec3 &p = mesh_tgt_tr.position (*v_it);
            mesh_tgt_tr.position (*v_it) = transformations_[1].transformPoint (p);

            const Vec3 &n = mesh_tgt_tr.get_vertex_property<Normal> ("v:normal")[*v_it];
            mesh_tgt_tr.vertex_property<Normal> ("v:normal")[*v_it] = transformations_[1].transformNormal (n);
        }

        // And write it to hdd
        surface_mesh::write_obj (mesh_tgt_tr, "pointcloud_transformed.obj");

        printf ("Wrote meshes to fitted_template.obj and pointcloud_transformed.obj\n");
        break;
    }
    case 'h':
    {
        printf("Help:\n");
        printf("SHIFT and move mouse: manual alignment\n");
        printf("'h'\t-\thelp\n");
        printf("'n'\t-\tnext mesh\n");
        printf("'r'\t-\tregister current mesh selected mesh using point-2-point optimization\n");
        printf("' '\t-\tregister current mesh selected mesh using point-2-surface optimization\n");
        printf("'s'\t-\tsave points to output\n");
        break;
    }
    default:
    {
        GlutExaminer::keyboard(key, x, y);
        break;
    }

    }

}



//-----------------------------------------------------------------------------


void
RegistrationViewer::
mouse(int button, int state, int x, int y)
{
    // manual object transformation when pressing SHIFT
    if( glutGetModifiers() & GLUT_ACTIVE_SHIFT ) mode_ = MOVE;
    else mode_ = VIEW;

    GlutExaminer::mouse(button, state, x, y);
}





//-----------------------------------------------------------------------------
/// called during mouse motion while button is pressed
void
RegistrationViewer::
motion(int x, int y)
{
    switch (mode_)
    {
    default:
    case VIEW:
    {
        GlutExaminer::motion(x, y);
        break;
    }
    case MOVE:
    {
        // manual object transformation when pressing SHIFT

        // zoom
        if (button_down_[0] && button_down_[1])
        {
            float dy = y - last_point_2D_[1];
            float h  = height_;
            Transformation mv_tr = Transformation::retrieve_gl();
            Transformation tr;
            tr.scale_ = 1. - dy / h * 5.;
            transformations_[src_index_] =  tr * transformations_[src_index_];
        }
        // rotation
        else if (button_down_[0])
        {
            if (last_point_ok_)
            {
                Vec2i  new_point_2D;
                Vec3f  new_point_3D;
                bool   new_point_ok;

                new_point_2D = Vec2i(x, y);
                new_point_ok = map_to_sphere(new_point_2D, new_point_3D);

                if (new_point_ok)
                {
                    Vec3f axis      = cross (last_point_3D_, new_point_3D);
                    float cos_angle = dot (last_point_3D_, new_point_3D);

                    if (fabs(cos_angle) < 1.0)
                    {
                        float angle = 2.0*acos(cos_angle);
                        //rotate(axis, angle);
                        Transformation mv_tr = Transformation::retrieve_gl();
                        mv_tr.translation_.fill(0);
                        Transformation tr(angle, Vec3f(axis[0],axis[1],axis[2]));

                        transformations_[src_index_] = mv_tr.inverse() * tr * mv_tr * transformations_[src_index_];
                    }
                }
            }
        }
        // translation
        else if (button_down_[1])
        {
            float dx = x - last_point_2D_[0];
            float dy = y - last_point_2D_[1];
            printf ("translate dx %f, dy %f --- x %d y %d last_point %d %d\n", dx, dy, x, y, last_point_2D_[0], last_point_2D_[1]);

            float z = - ((modelview_matrix_[ 2]*center_[0] +
                    modelview_matrix_[ 6]*center_[1] +
                    modelview_matrix_[10]*center_[2] +
                    modelview_matrix_[14]) /
                    (modelview_matrix_[ 3]*center_[0] +
                    modelview_matrix_[ 7]*center_[1] +
                    modelview_matrix_[11]*center_[2] +
                    modelview_matrix_[15]));

            float aspect = (float)width_ / (float)height_;
            float up     = tan(fovy_/2.0f*M_PI/180.f) * near_;
            float right  = aspect*up;

            Transformation mv_tr = Transformation::retrieve_gl();
            Transformation tr(2.0*dx/width_*right/near_*z, -2.0*dy/height_*up/near_*z, 0.0f);
            transformations_[src_index_] = mv_tr.inverse() * tr * mv_tr * transformations_[src_index_];
        }


        // remeber points
        last_point_2D_ = Vec2i(x, y);
        last_point_ok_ = map_to_sphere(last_point_2D_, last_point_3D_);

        glutPostRedisplay();
        break;
    }
    }
}
