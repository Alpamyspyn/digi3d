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
//=============================================================================
//
//  CLASS ValenceViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "ValenceViewer.hh"
#include <vector>
#include <algorithm>
#include <stdio.h>



//== IMPLEMENTATION ==========================================================


ValenceViewer::
ValenceViewer(const char* _title, int _width, int _height)
  : MeshViewer(_title, _width, _height)
{
  add_draw_mode("Vertex Valences");
}


//-----------------------------------------------------------------------------


ValenceViewer::
~ValenceViewer()
{
}

//-----------------------------------------------------------------------------

bool
ValenceViewer::
open_mesh(const char* _filename)
{
  // load mesh
  if (MeshViewer::open_mesh(_filename))
  {
    // compute vertex valence and color coding
    calc_valences();
    color_coding();

    glutPostRedisplay();
    return true;
  }
  return false;
}


//-----------------------------------------------------------------------------


void
ValenceViewer::
calc_valences()
{
    // initialize the valence property of the mesh
    vvalence_ = mesh_.vertex_property<unsigned int>("v:valence", 0);
    int vcount;

    // EXERCISE 1.2 /////////////////////////////////////////////////////////////
    // Compute valence of every vertex of "mesh_" and store them in vertex property vvalence_
    // (hint: use Mesh::Vertex_iterator)

    //iterate over all vertices
    Mesh::Vertex_iterator v_it, v_begin, v_end;
    v_begin = mesh_.vertices_begin();
    v_end = mesh_.vertices_end();
    for( v_it = v_begin ; v_it != v_end; ++ v_it )
    {

         Mesh::Vertex v = *v_it;
         vcount = 0;
         // iterate through neighbourhood
         //Mesh::Vertex v0;
         Mesh::Vertex_around_vertex_circulator vc, vc_end;
         vc = mesh_.vertices(v);
         vc_end = vc;
         do {
          vcount++;
         } while(++vc != vc_end);
         vvalence_[v] = vcount;
    }

    /////////////////////////////////////////////////////////////////////////////
}


//-----------------------------------------------------------------------------


void
ValenceViewer::
color_coding()
{
    // EXERCISE 1.3 /////////////////////////////////////////////////////////////
    // Implement a color visualization of your choice that shows the valence of
    // each vertex of "mesh_". Store the color of each vertex using the method
    // set_color(Vertex, Color)
    Mesh::Vertex_property<Color> colors;
    colors = mesh_.vertex_property<Color>("v:color");
    Mesh::Vertex_iterator v_it, v_begin, v_end;
    v_begin = mesh_.vertices_begin();
    v_end = mesh_.vertices_end();
    unsigned int valmin = MAXINT;
    unsigned int valmax = 0;
    unsigned int valsec = 0;
    for (v_it = v_begin; v_it != v_end; v_it++)
    {
        Mesh::Vertex v = *v_it;
        if (vvalence_[v] < valmin && vvalence_[v] > 0 ) { valmin = vvalence_[v]; }
        if (vvalence_[v] > valmax) {
            valsec = valmax;
            valmax = vvalence_[v];
        }
    }
    float fractv;

    for( v_it = v_begin ; v_it != v_end; ++ v_it )
    {
         Mesh::Vertex v = *v_it;
         fractv = ((float)vvalence_[v] - (float)valmin)/float(valsec-valmin);
         float blue = std::min(std::max((0.75-fractv), 0.0), 1.0);
         float red = std::min(std::max((fractv-0.25), 0.0), 1.0);
         float green = std::min(std::max((std::abs(fractv-0.5))-1.0, 0.0), 1.0);
         colors[v] = Color((int)(red*255),(int)(green*255),(int)(blue*255));
    }


    /////////////////////////////////////////////////////////////////////////////
}


//-----------------------------------------------------------------------------

//color ValenceViewer::valenceColor(unsigned int valence)
//{

//}


//-----------------------------------------------------------------------------

void
ValenceViewer::
draw(const std::string& _draw_mode)
{
    if (indices_.empty())
    {
        MeshViewer::draw(_draw_mode);
        return;
    }

    if (_draw_mode == "Vertex Valences")
    {
        glDisable(GL_LIGHTING);
        glShadeModel(GL_SMOOTH);
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
        GL::glVertexPointer(vertex_positions_pointer());
        GL::glNormalPointer(vertex_normals_pointer());
        GL::glColorPointer(vertex_colors_pointer());
        glDepthRange(0.01, 1.0);
        glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);
        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);
        glColor3f(0.1, 0.1, 0.1);
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

    else MeshViewer::draw(_draw_mode);
}


//=============================================================================
