#pragma once


#include "surface_mesh/Surface_mesh.h"
#include "Transformation.hh"
#include "eigen/Eigen/LU"
#include "eigen/Eigen/Geometry"
#include "../KdTree/KdTree.h"
#include "../eigen/Eigen/SparseCholesky"

typedef surface_mesh::Surface_mesh  Mesh;
typedef Eigen::Matrix<double, 3, Eigen::Dynamic> PointCloudMatrix;


void
initialize_alignment (const Mesh &mesh_src,
                      const Mesh &mesh_tgt,
                      Transformation &transf_tgt)
{
  /*
   * Perform the initial alignment of the two meshes as explained in the exercise sheet.
   */


  printf ("Initial target transformation: scale %f, translation %f %f %f\n",
          transf_tgt.scale_, transf_tgt.translation_(0), transf_tgt.translation_(1), transf_tgt.translation_(2));
}


void
compute_correspondences (const PointCloudMatrix &verts_src,
                         const PointCloudMatrix &verts_tgt,
                         const std::vector<size_t> &indices_src,
                         std::vector<std::pair<size_t, size_t> > &corresps)
{
  /*
  /// Brute-force
  for (size_t i = 0; i < indices_src.size (); ++i)
  {
    size_t i_src = indices_src[i];
    double closest_dist = std::numeric_limits<double>::max ();
    size_t closest_index = 0;
    for (size_t i_tgt = 0; i_tgt < verts_tgt.cols (); i_tgt += 100)
    {
      double dist = (verts_src.col (i_src) - verts_tgt.col (i_tgt)).norm ();
      if (dist < closest_dist)
      {
        closest_dist = dist;
        closest_index = i_tgt;
      }
    }

    corresps.push_back (std::make_pair (i_src, closest_index));
  }
  */

  /// Use a Kd-tree -> much faster
  ClosestPoints<double> closest_points (verts_tgt);
  for (size_t i = 0; i < indices_src.size (); ++i)
  {
    size_t i_src = indices_src[i];
    std::vector<size_t> nn_indices = closest_points.knn_search (verts_src.col (i_src), 1);
    corresps.push_back (std::make_pair (i_src, nn_indices.front ()));
  }
}


void
filter_correspondences (const PointCloudMatrix &verts_src,
                        const PointCloudMatrix &verts_tgt,
                        const PointCloudMatrix &normals_src,
                        const PointCloudMatrix &normals_tgt,
                        const std::vector<std::pair<size_t, size_t> > &corresps,
                        std::vector<std::pair<size_t, size_t> > &corresps_filtered)
{
  /*
   * This is where the correspondence filtering code should go in.
   * Filter the correspondences by their distnace and by normal compatibility.
   */


  printf ("--> filtering corresps: before %d, after %d\n", corresps.size (), corresps_filtered.size ());
}



void
compute_incremental_transformation (const PointCloudMatrix &vertices_src,
                                    const PointCloudMatrix &vertices_tgt,
                                    const PointCloudMatrix &normals_tgt,
                                    const std::vector<std::pair<size_t, size_t> > &corresps,
                                    Transformation &transform_inc)
{
  Eigen::Matrix<double, Eigen::Dynamic, 7> J;
  Eigen::VectorXd b;


  /*
   * Build the jacobian and the residual of the linearized system using the point to plane
   * metric on the given correspondences.
   * The linear system is built and solved for you, and the incremental transformation is
   * extracted from the solution vector.
   */



  Eigen::Matrix<double, 7, 7> JtJ = J.transpose () * J;
  Eigen::Matrix<double, 7, 1> Jtb = J.transpose () * b;
  Eigen::Matrix<double, 7, 1> var_result = JtJ.fullPivLu ().solve (-Jtb);

  transform_inc.scale_ = (var_result (6) + 1.) * (var_result (6) + 1.);
  transform_inc.rotation_ = Eigen::AngleAxisd (var_result (0), Eigen::Vector3d::UnitX ()) *
                            Eigen::AngleAxisd (var_result (1), Eigen::Vector3d::UnitY ()) *
                            Eigen::AngleAxisd (var_result (2), Eigen::Vector3d::UnitZ ()).matrix ();
  transform_inc.translation_ = var_result.block<3, 1> (3, 0);
}


void
extract_vertices_and_normals_from_mesh (const Mesh &mesh,
                                        PointCloudMatrix &vertices,
                                        PointCloudMatrix &normals)
{
  vertices.resize (3, mesh.n_vertices ());
  normals.resize (3, mesh.n_vertices ());

  for (Mesh::Vertex_iterator v_it = mesh.vertices_begin ();
       v_it != mesh.vertices_end (); ++v_it)
  {
    const Vec3f p = mesh.position (*v_it);
    vertices.col ((*v_it).idx ()) = Eigen::Vector3d (p[0], p[1], p[2]);

    const Vec3f n = mesh.get_vertex_property<Normal>("v:normal") [*v_it];
    normals.col ((*v_it).idx ()) = Eigen::Vector3d (n[0], n[1], n[2]);
  }
}


void
apply_transformation (const Transformation &tr,
                      PointCloudMatrix &vertices,
                      PointCloudMatrix &normals)
{
  for (size_t i = 0; i < vertices.cols (); ++i)
  {
    vertices.col (i) = tr.transformPoint (vertices.col (i));
    normals.col (i) = tr.transformNormal (normals.col (i));
  }
}


bool
rigid_registration (const Mesh &mesh_src,
                    const Mesh &mesh_tgt,
                    const std::vector<size_t> &indices_src,
                    const Transformation &transform_src,
                    Transformation &transform_tgt)
{
  printf ("Transforming the points\n");
  /// Transform the two sets of points
  PointCloudMatrix vertices_src (3, mesh_src.n_vertices ()),
                   vertices_tgt (3, mesh_tgt.n_vertices ()),
                   normals_src (3, mesh_src.n_vertices ()),
                   normals_tgt (3, mesh_tgt.n_vertices ());

  extract_vertices_and_normals_from_mesh (mesh_src, vertices_src, normals_src);
  extract_vertices_and_normals_from_mesh (mesh_tgt, vertices_tgt, normals_tgt);
  apply_transformation (transform_src, vertices_src, normals_src);
  apply_transformation (transform_tgt, vertices_tgt, normals_tgt);

  printf ("Computing correspondences\n");
  /// Compute correspondences
  std::vector<std::pair<size_t, size_t> > corresps;
  compute_correspondences (vertices_src, vertices_tgt,
                           indices_src,
                           corresps);

  printf ("Filtering correspondences\n");
  std::vector<std::pair<size_t, size_t> > corresps_filtered;
  filter_correspondences (vertices_src, vertices_tgt,
                          normals_src, normals_tgt,
                          corresps, corresps_filtered);

  Transformation transform_inc;
  compute_incremental_transformation (vertices_src, vertices_tgt, normals_tgt,
                                      corresps_filtered, transform_inc);

  transform_tgt = transform_inc.inverse () * transform_tgt;

  return (true);
}
