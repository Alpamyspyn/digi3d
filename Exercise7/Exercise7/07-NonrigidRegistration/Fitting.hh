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

    Mesh::Vertex_iterator v_it, v_end(mesh_src.vertices_end());

    //centers
    Point src_center(0.0, 0.0, 0.0);
    for(v_it = mesh_src.vertices_begin(); v_it != v_end; ++v_it)
        src_center += mesh_src.position(*v_it);

    src_center /= mesh_src.n_vertices();

    Mesh::Vertex_iterator v_it2, v_end2(mesh_tgt.vertices_end());

    Point tgt_center(0.0, 0.0, 0.0);
    for(v_it2 = mesh_tgt.vertices_begin(); v_it2 != v_end2; ++v_it2)
        tgt_center += mesh_tgt.position(*v_it2);

    tgt_center /= mesh_tgt.n_vertices();

    //scaling
    double src_dist = 0.0;
    for(v_it = mesh_src.vertices_begin(); v_it != v_end; ++v_it)
        src_dist = std::max(src_dist, (double) norm((mesh_src.position(*v_it) - src_center)));

    double tgt_dist = 0.0;
    for(v_it2 = mesh_tgt.vertices_begin(); v_it2 != v_end2; ++v_it2)
        tgt_dist = std::max(tgt_dist, (double) norm((mesh_tgt.position(*v_it2) - tgt_center)));


    transf_tgt.scale_ = src_dist / tgt_dist;
    transf_tgt.translation_(0) = transf_tgt.scale_ * (src_center[0] - tgt_center[0]);
    transf_tgt.translation_(1) = transf_tgt.scale_ * (src_center[1] - tgt_center[1]);
    transf_tgt.translation_(2) = transf_tgt.scale_ * (src_center[2] - tgt_center[2]);

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

    for(size_t i = 0; i < corresps.size(); ++i)
    {
        Vector3d v_src = verts_src.col(corresps.at(i).first);
        Vector3d v_tgt = verts_tgt.col(corresps.at(i).second);
        Vector3d n_src = normals_src.col(corresps.at(i).first);
        Vector3d n_tgt = normals_tgt.col(corresps.at(i).second);

        double angle = acos(n_src.dot(n_tgt) / n_src.norm() * n_tgt.norm());

        if((v_src - v_tgt).norm() < 0.05 && angle < 0.5236)
            corresps_filtered.push_back(corresps.at(i));
    }


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

  Eigen::Vector3d src_center, tgt_center, translation;
  Eigen::Matrix3d rotation, A, M, N;
  double scaling = 0.0;
  size_t n_corresps = corresps.size();

  b.resize(n_corresps);
  J.resize(n_corresps, 7);

  for(size_t i = 0; i < n_corresps; ++i)
  {
      src_center += vertices_src.col(corresps.at(i).first);
      tgt_center += vertices_tgt.col(corresps.at(i).second);
  }

  src_center /= n_corresps;
  tgt_center /= n_corresps;

  for(size_t i = 0; i < n_corresps; ++i)
  {
      M += (vertices_src.col(corresps.at(i).first) - src_center) * (vertices_tgt.col(corresps.at(i).second) - tgt_center).transpose();
      N += (vertices_tgt.col(corresps.at(i).second) - tgt_center) * (vertices_tgt.col(corresps.at(i).second) - tgt_center).transpose();
  }

  A = M * N.inverse();

  Eigen::JacobiSVD<Eigen::Matrix3d> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

  rotation = svd.matrixU() * svd.matrixV();
  translation = src_center  - rotation * tgt_center;

  Eigen::Matrix3d rot_x, rot_y, rot_z;

  rot_x << 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0;
  rot_y << 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0;
  rot_z << 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;

  for(size_t i = 0; i < n_corresps; ++i)
  {
      b(i) = normals_tgt.col(corresps.at(i).second).transpose() * (vertices_src.col(corresps.at(i).first) - vertices_tgt.col(corresps.at(i).second));

      // rotation
      J(i,0) = normals_tgt.col(corresps.at(i).second).transpose() * rot_x * vertices_src.col(corresps.at(i).first);
      J(i,1) = normals_tgt.col(corresps.at(i).second).transpose() * rot_y * vertices_src.col(corresps.at(i).first);
      J(i,2) = normals_tgt.col(corresps.at(i).second).transpose() * rot_z * vertices_src.col(corresps.at(i).first);

      // translation
      J(i,3) = normals_tgt.col(corresps.at(i).second).transpose() * Eigen::Vector3d::UnitX();
      J(i,4) = normals_tgt.col(corresps.at(i).second).transpose() * Eigen::Vector3d::UnitY();
      J(i,5) = normals_tgt.col(corresps.at(i).second).transpose() * Eigen::Vector3d::UnitZ();

      // scaling
      J(i,6) = 2.0 * normals_tgt.col(corresps.at(i).second).transpose() * rotation * vertices_src.col(corresps.at(i).first);
  }

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



void
compute_displacements (const Mesh &mesh_src,
                       const PointCloudMatrix &vertices_src,
                       const PointCloudMatrix &vertices_tgt,
                       const PointCloudMatrix &normals_tgt,
                       const std::vector<std::pair<size_t, size_t> > &corresps,
                       const Eigen::VectorXd &d_previous,
                       Eigen::VectorXd &d)
{


    const double weight_closeness = 1e-2;
    const double weight_smoothness = 1.0;

    /*
     * Build the jacobian and the residual of the linearized system using the point to plane
     * metric on the given correspondences.
     * The linear system is built and solved for you, and the incremental transformation is
     * extracted from the solution vector.
     */

    std::cout << "  Distance term"<< std::endl;
    Eigen::SparseMatrix<double> J_distance (corresps.size (), 3 * vertices_src.cols ());
    Eigen::VectorXd b_distance (corresps.size ());
    std::vector<Eigen::Triplet<double> > triplets;
    for (size_t i = 0; i < corresps.size(); ++i)
    {
        size_t c1 = corresps.at(i).first;
        size_t c2 = corresps.at(i).second;
        Vector3d dist = vertices_src.col(c1) - vertices_tgt.col(c2);
        Vector3d jx(1,dist(1),dist(2));
        Vector3d jy(dist(0),1,dist(2));
        Vector3d jz(dist(0),dist(1),1);
      //  Vector3d v(normals_tgt.col(c2).transpose()*jx,normals_tgt.col(c2).transpose()*jy,normals_tgt.col(c2).transpose()*jz);
        triplets.push_back(Eigen::Triplet<double>(i, c1*3, normals_tgt.col(c2).transpose()*jx));
        triplets.push_back(Eigen::Triplet<double>(i,c1*3+1,normals_tgt.col(c2).transpose()*jy));
        triplets.push_back(Eigen::Triplet<double>(i,c1*3+2,normals_tgt.col(c2).transpose()*jz));

    }
    J_distance.setFromTriplets(triplets.begin(), triplets.end());

    for (size_t i = 0; i<corresps.size(); ++i)
    {
        b_distance(i) = normals_tgt.col(corresps.at(i).second).transpose() * (vertices_src.col(corresps.at(i).first) - vertices_tgt.col(corresps.at(i).second));
    }


    std::cout << "  Closeness term"<< std::endl;
    Eigen::SparseMatrix<double> J_closeness (3 * vertices_src.cols (), 3 * vertices_src.cols ());

    J_closeness = Eigen::MatrixXd::Identity(3 * vertices_src.cols (), 3 * vertices_src.cols ());

    std::cout << "  Smoothness term"<< std::endl;
    Eigen::SparseMatrix<double> J_smoothness (3 * vertices_src.cols (), 3 * vertices_src.cols ());

    Eigen::Matrix3d submatrix = MatrixXd::Constant(3, 3, 1.0);

    Mesh::Vertex_around_vertex_circulator vv_c, vv_end;
    for(size_t i=0; i<vertices_src.cols(); ++i)
    {
        for(size_t j=0; j<vertices_src.cols(); ++i)
        {
            if (i==j)
            {
                size_t n_size = 0;
                vv_c = vertices_src.col(i);
                vv_end = vv_c;
                do{
                        n_size++;
                  }while(++vv_c != vv_end);
                J_smoothness.block(i,j,3,3) = n_size*submatrix;
                break;
            }

            vv_c = vertices_src.col(i);
            vv_end = vv_c;
            do{
                if(vv_c == vertices_src.col(j))
                {
                    J_smoothness.block(i,j,3,3) = submatrix;
                    break;
                }
            }while(++vv_c != vv_end);
            J_smoothness.block(i,j,3,3) = 0*submatrix;
        }
    }


    std::cout << "  Building the linear system"<< std::endl;
    Eigen::SparseMatrix<double> JtJ = J_distance.transpose () * J_distance
            + weight_closeness * J_closeness
            + weight_smoothness * J_smoothness.transpose () * J_smoothness;

    Eigen::VectorXd Jtb =  - J_distance.transpose () * b_distance;
    // At previous iteration J' * b_previous =  J' * J * d_previous
    Eigen::VectorXd Jtb_previous = - J_distance.transpose () * J_distance * d_previous;

    Jtb -= Jtb_previous;

    std::cout << "  Solving the linear system"<< std::endl;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute (JtJ);
    if (solver.info () != Eigen::Success)
        printf ("Solver decomposition failed.\n");
    d = solver.solve (Jtb);
}


bool
nonrigid_registration (Mesh &mesh_src,
                       const Mesh &mesh_tgt,
                       const std::vector<size_t> &indices_src,
                       const Transformation &transform_src,
                       const Transformation &transform_tgt,
                       const Eigen::VectorXd &d_previous,
                       Eigen::VectorXd &d)
{
    /// Transform the two sets of points
    PointCloudMatrix vertices_src (3, mesh_src.n_vertices ()),
            vertices_tgt (3, mesh_tgt.n_vertices ()),
            normals_src (3, mesh_src.n_vertices ()),
            normals_tgt (3, mesh_tgt.n_vertices ());

    extract_vertices_and_normals_from_mesh (mesh_src, vertices_src, normals_src);
    extract_vertices_and_normals_from_mesh (mesh_tgt, vertices_tgt, normals_tgt);
    apply_transformation (transform_src, vertices_src, normals_src);
    apply_transformation (transform_tgt, vertices_tgt, normals_tgt);

    std::cout << "Computing correspondences"<< std::endl;
    std::vector<std::pair<size_t, size_t> > corresps;
    compute_correspondences (vertices_src, vertices_tgt,
                             indices_src,
                             corresps);

    std::cout << "Filtering correspondences"<< std::endl;
    std::vector<std::pair<size_t, size_t> > corresps_filtered;
    filter_correspondences (vertices_src, vertices_tgt,
                            normals_src, normals_tgt,
                            corresps, corresps_filtered);

    std::cout << "Computing displacements"<< std::endl;
    compute_displacements (mesh_src, vertices_src, vertices_tgt, normals_tgt, corresps_filtered,
                           d_previous, d);

    std::cout << "Updating the source mesh" << std::endl << std::endl;
    for (Mesh::Vertex_iterator v_it = mesh_src.vertices_begin ();
         v_it != mesh_src.vertices_end (); ++v_it)
    {
        Eigen::Vector3d v_eig = vertices_src.col ((*v_it).idx ()) -
                d_previous.block<3, 1> (3 * (*v_it).idx (), 0) +
                d.block<3, 1> (3 * (*v_it).idx (), 0);
        mesh_src.position (*v_it) = Vec3f (v_eig (0), v_eig (1), v_eig (2));
    }

    return true;
}
