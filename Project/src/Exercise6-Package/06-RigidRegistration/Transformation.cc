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
//  CLASS Transformation - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "Transformation.hh"
#include "../GLUTViewer/gl.hh"


//== IMPLEMENTATION ==========================================================


Transformation::
Transformation()
{
	set_identity();
}


//=============================================================================
Transformation::
Transformation(float tx, float ty, float tz)
{
	set_identity();
  translation_ = Eigen::Vector3d(tx,ty,tz);
}


//=============================================================================
Transformation::
Transformation(float angle, Vec3f axis)
{
	set_identity();

  double l = norm( axis );
	if (l > 0) {
		double l1 = 1.0/l;
		double x = axis[0]*l1;
		double y = axis[1]*l1;
		double z = axis[2]*l1;
		double s = sin(angle), c = cos(angle);
		double xs = x*s, ys = y*s, zs = z*s, c1 = 1.0-c;
		double xx = c1*x*x, yy = c1*y*y, zz = c1*z*z;
		double xy = c1*x*y, xz = c1*x*z, yz = c1*y*z;
    rotation_(0, 0) = xx+c;
    rotation_(0, 1) = xy-zs;
    rotation_(0, 2) = xz+ys;
    rotation_(1, 0) = xy+zs;
    rotation_(1, 1) = yy+c;
    rotation_(1, 2) = yz-xs;
    rotation_(2, 0) = xz-ys;
    rotation_(2, 1) = yz+xs;
    rotation_(2, 2) = zz+c;
	}
}



//=============================================================================

void
Transformation::
set_identity()
{
  scale_ = 1.;
  rotation_.setIdentity();
  translation_.setZero ();
}


//=============================================================================

Transformation
Transformation::
operator*( const Transformation & o )
{
	Transformation t;

  t.scale_ = scale_ * o.scale_;
	t.rotation_ = rotation_ * o.rotation_;
  t.translation_ = scale_ * rotation_ * o.translation_ + translation_;

	return t;
}


//=============================================================================

// inverse rigid motion
Transformation
Transformation::
inverse() {
	Transformation t;
  t.scale_ = 1. / scale_;
	t.rotation_ = rotation_.transpose();
  t.translation_ = - 1. / scale_ * t.rotation_ * translation_;
	return t;
}

//=============================================================================


// Transform point
Vec3f Transformation::transformPoint( const Vec3f & p ) const
{
  Eigen::Vector3d p_eigen (p[0], p[1], p[2]);
  Eigen::Vector3d aux = scale_ * rotation_ * p_eigen + translation_;
  return (Vec3f (aux(0), aux(1), aux(2)));
}


Eigen::Vector3d Transformation::transformPoint( const Eigen::Vector3d & p ) const
{
  Eigen::Vector3d p_tr = scale_ * rotation_ * p + translation_;
  return (p_tr);
}


// Transform normal
Vec3f Transformation::transformNormal( const Vec3f & n ) const
{
  Eigen::Vector3d n_eigen (n[0], n[1], n[2]);
  Eigen::Vector3d aux = rotation_ * n_eigen;
  return (Vec3f (aux(0), aux(1), aux(2)));
}


Eigen::Vector3d Transformation::transformNormal( const Eigen::Vector3d & n ) const
{
  Eigen::Vector3d n_tr = rotation_ * n;
  return (n_tr);
}


//=============================================================================

/// apply transformation to current OpenGL matrix
void
Transformation::
apply_gl()
{
	float data[16];
	memset(data,0,sizeof(float)*16);
	data[15] = 1;

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
      data[4*j+i] = scale_ * rotation_(i, j);
	for(int i = 0; i < 3; i++)
    data[12+i] = translation_(i);

	glMultMatrixf(data);
}


//=============================================================================

/// retrieve current OpenGL transformation
Transformation
Transformation::
retrieve_gl()
{
	Transformation tr;

	double data[16];
	glGetDoublev( GL_MODELVIEW_MATRIX, data);

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
      tr.rotation_(i, j) = data[4*j+i];
	for(int i = 0; i < 3; i++)
    tr.translation_(i) = data[12+i];


  /// Extract the scaling
  tr.scale_ = tr.rotation_.col (0).norm ();
  tr.rotation_ /= tr.scale_;

	return tr;
}

//=============================================================================
