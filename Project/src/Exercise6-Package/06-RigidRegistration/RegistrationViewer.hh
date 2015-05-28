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
//  CLASS RegistrationViewer
//
//=============================================================================


#ifndef REGISTRATIONVIEWERWIDGET_HH
#define REGISTRATIONVIEWERWIDGET_HH


//== INCLUDES =================================================================


#include "GLUTViewer/GlutExaminer.hh"
#include "Transformation.hh"
#include "surface_mesh/Surface_mesh.h"

//== CLASS DEFINITION =========================================================



class RegistrationViewer : public GlutExaminer
{

public:

  /// default constructor
  RegistrationViewer (const char* _title, int _width, int _height);

  // destructor
  ~RegistrationViewer ();

  /// set output filename
  void set_output (const std::string & filename);

  /// open meshes
  bool open_meshes (const std::vector<std::string> & _filenames,
                    const std::vector<size_t> &indices_template);

  void initialize ();


protected:
  typedef surface_mesh::Surface_mesh  Mesh;

  virtual void draw (const std::string& _draw_mode);
  virtual void keyboard (int key, int x, int y);
  virtual void motion (int x, int y);
  virtual void mouse (int button, int state, int x, int y);

private:

  /// update buffer with face indices
  void update_face_indices ();

  /// draw the mesh in scene
  virtual void draw (int index, const Vec3f & color);


protected:


  enum Mode { VIEW, MOVE } mode_;

protected:

	std::string								outputFilename_;

  float										averageVertexDistance_;
  int src_index_ = 0;
  int 										numProcessed_;
  std::vector< Mesh >                       meshes_;
  std::vector< std::vector<unsigned int> >  indices_;
  std::vector< Transformation >				transformations_;

  std::vector<size_t> indices_template_;
};


//=============================================================================
#endif // REGISTRATIONVIEWERWIDGET_HH defined
//=============================================================================

