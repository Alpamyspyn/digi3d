#include "animation_loader.h"
#include <iostream>
#include <string>
#include "surface_mesh/Surface_mesh.h"
#include "surface_mesh/IO.h"
#include "RegistrationViewer.hh"

bool
readIndicesFile (const std::string &filename, std::vector<size_t> &indices)
{
  std::ifstream file (filename.c_str ());
  if (!file.is_open ())
  {
    printf ("Error reading indices file %s.\n", filename.c_str ());
    return (false);
  }
  size_t points_size;
  file >> points_size;
  for (size_t i = 0; i < points_size; ++i)
  {
    size_t v;
    file >> v;
    indices.push_back (v);
  }
  file.close ();
  return (true);
}

int
main (int argc,
      char **argv)
{
  typedef surface_mesh::Surface_mesh  Mesh;

  /// Read the fsb file given in argv[1]
  std::vector<fs::fsMsgTrackingState> tracking_data;
  loadFaceshiftAnimation (argv[1], tracking_data);
  /// Read folder containing blendshapes models in argv[2]
  std::vector<Mesh> bs;
  Mesh neutral, result;
  std::string path(argv[2]);
  surface_mesh::read_mesh (neutral, path+"/neutral.obj");

  //printf ("%d vertices\n",neutral.n_vertices());


  for (int i = 0; i < 49 ; ++i)
  {
      std::string current_bs=path+"/"+std::to_string(i)+".obj";

      //std::cout<<current_bs<<std::endl;
      Mesh bs_read;
      surface_mesh::read_mesh (bs_read, current_bs);
      bs.push_back(bs_read);
  }



  /// Example of using the tracking data
  const size_t num_frames = tracking_data.size ();
  for (size_t f_i = 0; f_i < num_frames; ++f_i)
  {

      surface_mesh::read_mesh (result, path+"/neutral.obj");

    printf ("Frame %ld / %ld\n", f_i, num_frames);

   /* /// Rigid transformation
    printf ("rotation as a quaternion: %f %f %f %f\n",
            tracking_data[f_i].tracking_data ().m_headRotation.x,
            tracking_data[f_i].tracking_data ().m_headRotation.y,
            tracking_data[f_i].tracking_data ().m_headRotation.z,
            tracking_data[f_i].tracking_data ().m_headRotation.w);
    printf ("translation: %f %f %f\n",
            tracking_data[f_i].tracking_data ().m_headTranslation.x,
            tracking_data[f_i].tracking_data ().m_headTranslation.y,
            tracking_data[f_i].tracking_data ().m_headTranslation.z);
*/
    /// Blendshape weights
    const int num_bs = tracking_data[f_i].tracking_data ().m_coeffs.size ();
    for (size_t i = 0; i < num_bs; ++i)
    {
        const int num_bs = tracking_data[0].tracking_data ().m_coeffs.size ();
        for(size_t i = 0; i < num_bs; ++i)
        {
            Surface_mesh::Vertex_iterator v_it_neutral = neutral.vertices_begin();
            Surface_mesh::Vertex_iterator v_it_bs = bs[i].vertices_begin();
            Surface_mesh::Vertex_iterator v_it_result = result.vertices_begin();

            for(v_it_bs = bs[i].vertices_begin(); v_it_bs != bs[i].vertices_end(); ++v_it_bs )
            {
                Point x = tracking_data[0].tracking_data().m_coeffs[i]* (bs[i].position(*v_it_bs) - neutral.position(*v_it_neutral));
                result.position(*v_it_result) = result.position(*v_it_result) + x;
                ++v_it_neutral;
                ++v_it_result;
            }
        }
    }
    printf ("\n\n");

    surface_mesh::write_obj(result, "test"+std::to_string(f_i)+".obj");

  }





  glutInit (&argc, argv);

  RegistrationViewer window ("Registration", 512, 512);

  std::vector<std::string> files;
  /// Template mesh
  files.push_back ("test0.obj");

  std::vector<size_t> indices;
  readIndicesFile (argv[3], indices);

  if (argc > 1)
  {
    window.open_meshes (files, indices);
    window.initialize ();
  }

  glutMainLoop ();

  /*

  */

  return (0);
}
