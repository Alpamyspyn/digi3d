#include "../surface_mesh/Surface_mesh.h"
#include "../surface_mesh/IO.h"

#include "../eigen/eigen_include.h"
#include <iostream>

typedef surface_mesh::Surface_mesh  Mesh;


void
copyPositions (const Mesh &src,
               Mesh &tgt)
{
  for (Mesh::Vertex_iterator src_it = src.vertices_begin (),
       tgt_it = tgt.vertices_begin (); src_it != src.vertices_end ();
       ++src_it, ++tgt_it)
    tgt.position (*tgt_it) = src.position (*src_it);
}

int main(int argc, char **argv)
{

  /// argv[1] -> macaw neutral
  /// argv[2] -> macaw expression
  /// argv[3] -> new person neutral


  /// Read in the meshes, triangulate the template, copy the topology to the rest
  Mesh template_neutral, template_exp, person_neutral, aux;
  surface_mesh::read_mesh (template_neutral, argv[1]);
  template_neutral.triangulate ();
  template_neutral.update_face_normals ();
  template_neutral.update_vertex_normals ();
  surface_mesh::read_mesh (aux, argv[2]);
  template_exp = template_neutral;
  copyPositions (aux, template_exp);
  surface_mesh::read_mesh (aux, argv[3]);
  person_neutral = template_neutral;
  copyPositions (aux, person_neutral);

  // Save the meshes for debug reasons
  surface_mesh::write_obj (template_neutral, "debug_template_neutral.obj");
  surface_mesh::write_obj (template_exp, "debug_template_exp.obj");
  surface_mesh::write_obj (person_neutral, "debug_person_neutral.obj");


  // Set up the linear system to solve
  int num_vertices = template_neutral.n_vertices ();
  int num_tris = template_neutral.n_faces ();
  Eigen::Matrix<double, 4, 3> M;
  M << -1., -1., -1.,
        1.,  0.,  0.,
        0.,  1.,  0.,
        0.,  0.,  1.;


  // Accumulate the template vertices for the neutral
  Eigen::Matrix<double, Eigen::Dynamic, 4> V_template_neutral (3 * num_tris, 4),
                                           V_template_exp (3 * num_tris, 4),
                                           V_person_neutral (3 * num_tris, 4);
  for (Mesh::Face_iterator f_it = template_neutral.faces_begin (); f_it != template_neutral.faces_end (); ++f_it)
  {
    std::vector<Eigen::Vector3d> vertices;
    for (Mesh::Vertex_around_face_circulator v_it = Mesh::Vertex_around_face_circulator (&template_neutral, *f_it).begin ();
         v_it != Mesh::Vertex_around_face_circulator (&template_neutral, *f_it).end (); ++v_it)
    {
      surface_mesh::Vec3f v = template_neutral.position (*v_it);
      vertices.push_back (Eigen::Vector3d (v[0], v[1], v[2]));
    }

    for (size_t i = 0; i < 3; ++i)
      V_template_neutral.block<3, 1> (3 * (*f_it).idx (), i) = vertices[i];

    Eigen::Vector3d v_cross = (vertices[1] - vertices[0]).cross (vertices[2] - vertices[0]);
    Eigen::Vector3d v3 = vertices[0] + v_cross / v_cross.norm ();
    V_template_neutral.block<3, 1> (3 * (*f_it).idx (), 3) = v3;
  }

  // Accumulate the template vertices for the expression
  for (Mesh::Face_iterator f_it = template_exp.faces_begin (); f_it != template_exp.faces_end (); ++f_it)
  {
    std::vector<Eigen::Vector3d> vertices;
    for (Mesh::Vertex_around_face_circulator v_it = Mesh::Vertex_around_face_circulator (&template_exp, *f_it).begin ();
         v_it != Mesh::Vertex_around_face_circulator (&template_exp, *f_it).end (); ++v_it)
    {
      surface_mesh::Vec3f v = template_exp.position (*v_it);
      vertices.push_back (Eigen::Vector3d (v[0], v[1], v[2]));
    }

    for (size_t i = 0; i < 3; ++i)
      V_template_exp.block<3, 1> (3 * (*f_it).idx (), i) = vertices[i];

    Eigen::Vector3d v_cross = (vertices[1] - vertices[0]).cross (vertices[2] - vertices[0]);
    Eigen::Vector3d v3 = vertices[0] + v_cross / v_cross.norm ();
    V_template_exp.block<3, 1> (3 * (*f_it).idx (), 3) = v3;
  }

  // Accumulate the actor vertices for the neutral
  for (Mesh::Face_iterator f_it = person_neutral.faces_begin (); f_it != person_neutral.faces_end (); ++f_it)
  {
    std::vector<Eigen::Vector3d> vertices;
    for (Mesh::Vertex_around_face_circulator v_it = Mesh::Vertex_around_face_circulator (&person_neutral, *f_it).begin ();
         v_it != Mesh::Vertex_around_face_circulator (&person_neutral, *f_it).end (); ++v_it)
    {
      surface_mesh::Vec3f v = person_neutral.position (*v_it);
      vertices.push_back (Eigen::Vector3d (v[0], v[1], v[2]));
    }

    for (size_t i = 0; i < 3; ++i)
      V_person_neutral.block<3, 1> (3 * (*f_it).idx (), i) = vertices[i];

    Eigen::Vector3d v_cross = (vertices[1] - vertices[0]).cross (vertices[2] - vertices[0]);
    Eigen::Vector3d v3 = vertices[0] + v_cross / v_cross.norm ();
    V_person_neutral.block<3, 1> (3 * (*f_it).idx (), 3) = v3;
  }

  /// Compute the residual for the deformation energy.
  Eigen::VectorXd b (9 * num_tris);
  for (size_t f_i = 0; f_i < num_tris; ++f_i)
  {
    Eigen::Matrix3d b_block = -V_template_exp.block<3, 4> (3 * f_i, 0) * M *
                             (V_template_neutral.block<3, 4> (3 * f_i, 0) * M).inverse ();
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        b (9 * f_i + 3 * i + j) = b_block (i, j);
  }


  // Build the Jacobian and right-hand side of the regularization energy linearization
  Eigen::SparseMatrix<double> J (9 * num_tris, 3 * (num_vertices + num_tris));
  std::vector<Eigen::Triplet<double> > triplets;
  Eigen::VectorXd b_reg (Eigen::VectorXd::Zero (3 * (num_vertices + num_tris)));

  ///
  /// Your solution goes here
  ///

  // Make the b
      for (Mesh::Vertex_iterator v_it = person_neutral.vertices_begin();
           v_it != person_neutral.vertices_end (); ++v_it)
      {
          surface_mesh::Point pos = person_neutral.position(*v_it);
          b_reg[3*(*v_it).idx()] = -1.0*pos[0];
          b_reg[3*(*v_it).idx()+1] = -1.0*pos[1];
          b_reg[3*(*v_it).idx()+2] = -1.0*pos[2];
      }

      size_t N = person_neutral.n_vertices()+1;

      for (size_t i = N; i < num_tris; ++i)
      {
          b_reg[3*i] = -1.0*V_person_neutral(i-N,3);
          b_reg[3*i+1] = -1.0*V_person_neutral(i-N+1,3);
          b_reg[3*i+2] = -1.0*V_person_neutral(i-N+2,3);
      }

  // constructing the jacobian
  for(size_t f_i = 0; f_i != num_tris; ++f_i)
  {
      Eigen::MatrixXd X = M * (V_person_neutral.block<3, 4> (3 * f_i, 0) * M).inverse();

      // pushing ALL the triplets, block insertion would clearly make things too easy
      for(size_t x = 0; x < 2; x++)
      {
          for(size_t i = 0; i < 3; i++)
          {
              for(size_t j = 0; j < 4; j++)
              {
                 triplets.push_back(Eigen::Triplet<double>(f_i*9 + x*3 + i, x*4 + j, X(j, i)));
              }
          }
      }
  }

  J.setFromTriplets (triplets.begin (), triplets.end ());

  Eigen::SparseMatrix<double> Id (3 * (num_vertices + num_tris), 3 * (num_vertices + num_tris));
  Id.setIdentity ();

  printf ("Solving the sparse linear system ...\n");
  Eigen::SparseMatrix<double> JtJ = J.transpose () * J + Id;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  solver.compute (JtJ);
  printf ("Decomposing ...\n");
  if (solver.info () != Eigen::Success)
    printf ("Solver decomposition failed.\n");
  printf ("Solving ...\n");
  Eigen::VectorXd sol = solver.solve (J.transpose () * (-b) - b_reg);

  /// Put the solution back
  Mesh person_exp = person_neutral;
  for (Mesh::Vertex_iterator v_it = person_exp.vertices_begin (); v_it != person_exp.vertices_end (); ++v_it)
    person_exp.position (*v_it) = surface_mesh::Vec3f (sol (3 * (*v_it).idx () + 0),
                                                       sol (3 * (*v_it).idx () + 1),
                                                       sol (3 * (*v_it).idx () + 2));
  person_exp.update_face_normals ();
  person_exp.update_vertex_normals ();

  surface_mesh::write_obj (person_exp, "debug_person_exp.obj");

  return (0);
}
