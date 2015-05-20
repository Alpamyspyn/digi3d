#pragma once

#include "../eigen/Eigen/Core"
#include "../eigen/Eigen/Geometry"
#include <QColor>

#ifdef __APPLE__
#elif _WIN32
# include <GL/glew.h>
#endif

#include "../GLUTViewer/gl.hh"
#include "surface_mesh/Surface_mesh.h"

#include <QImage>

typedef surface_mesh::Surface_mesh  Mesh;

struct CameraParams
{
  double fx, fy, cx, cy;
  Eigen::Affine3d pose;
  int width, height;

  void
  print ()
  {
      std::cout << "Focals: " << fx << ", " << fy << "\nCenter: " << cx << ", " << cy << "\n" << "width " << width << ", " << height << std::endl;
      std::cout << "Pose:\n" << pose.matrix () << std::endl;
  }

  Eigen::Matrix3d
  getK ()
  {
      Eigen::Matrix3d K (Eigen::Matrix3d::Identity ());
      K (0, 0) = fx;
      K (1, 1) = fy;
      K (2, 2) = 1.;
      K (0, 2) = cx;
      K (1, 2) = cy;
      return (K);
  }
};


class SelectionBufferGeneratorFBO
{
public:
  SelectionBufferGeneratorFBO ()
  {}

  void
  setMesh (const Mesh &mesh)
  { mesh_ = mesh; }

  void
  setCameraParams (const CameraParams &params);

  void
  draw ();

  void
  getProjectedIds (Eigen::MatrixXi &buffer);

  void
  getRenderImage (QImage &render);

  void
  getVisibleTris (std::vector<bool> &visible_tris);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
  Mesh mesh_;
  Eigen::Matrix4f proj_matrix_;
  Eigen::Matrix4f model_matrix_;

  static GLuint fbo_id_, rb_id_, tex_id_;
  static int win_id_;
  int width_, height_;
};
