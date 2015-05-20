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
//  CLASS Transformation
//
//=============================================================================


#ifndef TRANSFORMATION_HPP_
#define TRANSFORMATION_HPP_

#include <vector>
#include "surface_mesh/Vector.h"
#include "../eigen/Eigen/Core"

using namespace surface_mesh;

/**
 * Transformation class
 *
 * contains a rotation and translation defining a rigid transformation
 */
class Transformation {

public:
    /// constructor: identity transformation
    Transformation();

    /// constructor: translation
    Transformation(float tx, float ty, float tz);

    /// constructor: rotation around axis
    Transformation(float angle, Vec3f axis);

    /// set identity transformation
    void set_identity();

    /// apply transformation to current OpenGL Matrix
    void apply_gl();

    /// retrieve curren OpenGL transformation
    static Transformation retrieve_gl();

    /// concatenate two transformations
    Transformation operator*( const Transformation & o );

    /// return inverse transformation
    Transformation inverse();

    /// Transform point
    Vec3f transformPoint( const Vec3f & p ) const;

    Eigen::Vector3d transformPoint (const Eigen::Vector3d &p) const;

    /// Transform a normal
    Vec3f transformNormal (const Vec3f &n) const;

    Eigen::Vector3d transformNormal (const Eigen::Vector3d &n) const;

    double scale_;
    Eigen::Matrix3d rotation_;
    Eigen::Vector3d translation_;
};

#endif /* TRANSFORMATION_HPP_ */
