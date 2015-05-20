#include <iostream>

#include <QDomDocument>
#include <QFile>
#include <QStringList>
#include <QDir>
#include <QImage>
#include <QPainter>

#include "../surface_mesh/Surface_mesh.h"
#include "../surface_mesh/IO.h"

#include "../eigen/eigen_include.h"

#include "VertexVisbility.hh"

#define TEXTURE_SIZE 1024


/*
 * reads the cameras.xml file from Photoscan
 * The file can be exported from Photoscan by: Tools->Export->Export Cameras
 */
bool
readPhotoscanData (const std::string &xml_filename,
                   const std::string &folder_name,
                   std::vector<CameraParams> &camera_params,
                   std::vector<QImage> &images)
{
  QDomDocument doc;
  QFile file (QString (xml_filename.c_str ()));
  if (!file.open(QIODevice::ReadOnly) || !doc.setContent(&file))
  {
    std::cerr << "Could not find the xml file at: " << xml_filename << std::endl;
    return (false);
  }

  CameraParams cam_params_init;
  QDomElement sensor_xml = doc.firstChildElement ("document").firstChildElement ("chunk").firstChildElement ("sensors").
      firstChildElement ("sensor").firstChildElement ("calibration");
  cam_params_init.fx = sensor_xml.firstChildElement ("fx").text ().toDouble ();
  cam_params_init.fy = sensor_xml.firstChildElement ("fy").text ().toDouble ();
  cam_params_init.cx = sensor_xml.firstChildElement ("cx").text ().toDouble ();
  cam_params_init.cy = sensor_xml.firstChildElement ("cy").text ().toDouble ();
  cam_params_init.width = sensor_xml.firstChildElement ("resolution").attributes ().namedItem ("width").nodeValue ().toInt ();
  cam_params_init.height = sensor_xml.firstChildElement ("resolution").attributes ().namedItem ("height").nodeValue ().toInt ();

  QDomNodeList camera_nodes = doc.firstChildElement ("document").firstChildElement ("chunk").firstChildElement ("cameras").elementsByTagName ("camera");
  std::cout << "Camera nodes count: " << camera_nodes.size () << std::endl;
  std::vector<QString> camera_image_names;
  for (size_t cam_i = 0; cam_i < camera_nodes.size (); ++cam_i)
  {
    QDomNode node = camera_nodes.item (cam_i);
    /// Check if it is enabled
    if (node.attributes ().namedItem ("enabled").nodeValue () == "false")
      continue;
    camera_image_names.push_back (node.attributes ().namedItem ("label").nodeValue ());

    QStringList transf_list = node.firstChildElement ("transform").text ().split (" ");
    Eigen::Matrix4d transform (Eigen::Matrix4d::Identity ());
    for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 4; ++j)
        transform (i, j) = transf_list[4 * i +j].toDouble ();

    CameraParams aux_params = cam_params_init;
    aux_params.pose = Eigen::Affine3d (transform);
    camera_params.push_back (aux_params);
  }

  std::cout << "Found " << camera_params.size () << " valid cameras. Proceeding to reading in the images.\n";

  QDir folder (QString (folder_name.c_str ()));
  for (size_t im_i = 0; im_i < camera_image_names.size (); ++im_i)
  {
    if (!folder.exists (camera_image_names[im_i]))
    {
      std::cerr << "Error: could not find image file " << camera_image_names[im_i].toStdString () << std::endl;
      return (false);
    }

    QImage image (folder.absoluteFilePath (camera_image_names[im_i]));
//    QImage image;
    images.push_back (image);
    std::cout << "   read image: " << camera_image_names[im_i].toStdString () << "\n";
  }

  return (true);
}


/*
 * Computes the texture space location and barycentric coordinates of each texture pixel corresponding
 * to each triangle of the input mesh
 */
void
computeTextureTrianglePixels (const Mesh &mesh,
                              std::vector<std::vector<Eigen::Vector2i> > &texture_triangle_pixels,
                              std::vector<std::vector<Eigen::Vector3d> > &texture_triangle_barys)
{
  /// For each triangle
  QImage texture_colors (TEXTURE_SIZE, TEXTURE_SIZE, QImage::Format_RGB888);
  texture_colors.fill (QColor (0, 0, 0));
  QPainter painter;
  for (Mesh::Face_iterator f_it = mesh.faces_begin (); f_it != mesh.faces_end (); ++f_it)
  {
    std::vector<Eigen::Vector2d> tex_coords;
    for (Mesh::Halfedge_around_face_circulator h_it = Mesh::Halfedge_around_face_circulator (&mesh, *f_it).begin ();
         h_it != Mesh::Halfedge_around_face_circulator (&mesh, *f_it).end (); ++h_it)
    {
      surface_mesh::Vec3f t = mesh.get_halfedge_property<surface_mesh::Texture_coordinate> ("h:texcoord") [*h_it];
      tex_coords.push_back (Eigen::Vector2d (t[0], t[1]));
    }

    std::vector<QPoint> tri_points;
    for (size_t i = 0; i < 3; ++i)
      tri_points.push_back (QPoint (tex_coords[i] (0) * static_cast<double> (TEXTURE_SIZE),
                                    tex_coords[i] (1) * static_cast<double> (TEXTURE_SIZE)));
    unsigned char r, g, b;
    size_t temp = (*f_it).idx () + 1; /// to avoid confusion with the black background
    r = temp / (256 * 256);
    g = (temp % (256 * 256)) / 256;
    b = temp % 256;

    painter.begin (&texture_colors);
    painter.setPen (QPen (QColor (r, g, b)));
    painter.setBrush (QBrush (QColor (r, g, b)));
    painter.drawPolygon (&tri_points[0], 3);
    painter.end ();
  }

  texture_triangle_pixels.resize (mesh.n_faces ());
  texture_triangle_barys.resize (mesh.n_faces ());
  for (int x = 0; x < texture_colors.width (); ++x)
    for (int y = 0; y < texture_colors.height (); ++y)
    {
      QRgb pix = texture_colors.pixel (x, y);
      int tri_id = qRed (pix) * 256 * 256 + qGreen (pix) * 256 + qBlue (pix);
      if (tri_id == 0)
        continue;

      tri_id --;
      Mesh::Face f_handle (tri_id);

      Eigen::Matrix3d U;
      int index = 0;
      for (Mesh::Halfedge_around_face_circulator h_it = Mesh::Halfedge_around_face_circulator (&mesh, f_handle).begin ();
           h_it != Mesh::Halfedge_around_face_circulator (&mesh, f_handle).end (); ++h_it)
      {
        surface_mesh::Vec3f t = mesh.get_halfedge_property<surface_mesh::Texture_coordinate> ("h:texcoord") [*h_it];
        for (size_t i = 0; i < 3; ++i)
          U (i, index) = t[i];
        index ++;
      }

      Eigen::Vector3d A = U.inverse () * Eigen::Vector3d (static_cast<double> (x) / static_cast<double> (TEXTURE_SIZE),
                                                          static_cast<double> (y) / static_cast<double> (TEXTURE_SIZE),
                                                          1.);

      texture_triangle_pixels[tri_id].push_back (Eigen::Vector2i (x, TEXTURE_SIZE - 1 - y));
      texture_triangle_barys[tri_id].push_back (A);
    }
}



int
main (int argc,
      char **argv)
{
  /// argv[1] - cameras.xml
  /// argv[2] - folder with undistorted images
  /// argv[3] - fitted mesh

  /// Read the camera parameters
  std::vector<CameraParams> camera_params;
  std::vector<QImage> images;
  if (!readPhotoscanData (argv[1], argv[2],
                          camera_params, images))
  {
    std::cerr << "Bad stuff happenened, bailing out.\n";
    return (-1);
  }


  /// Read the mesh
  Mesh mesh;
  surface_mesh::read_mesh (mesh, argv[3]);
  mesh.triangulate ();
  mesh.update_face_normals ();
  mesh.update_vertex_normals ();

  /// Prepare the texture info
  std::vector<std::vector<Eigen::Vector2i> > texture_triangle_pixels;
  std::vector<std::vector<Eigen::Vector3d> > texture_triangle_barys;
  computeTextureTrianglePixels (mesh, texture_triangle_pixels, texture_triangle_barys);


  /*
   * Data structure initialization to suggest a good way of approaching the problem.
   */
  QImage texture_best_view (TEXTURE_SIZE, TEXTURE_SIZE, QImage::Format_RGB888);
  texture_best_view.fill (QColor (0, 0, 0));
  std::vector<double> best_view_angles (TEXTURE_SIZE * TEXTURE_SIZE, 0.);

  QImage texture_average (TEXTURE_SIZE, TEXTURE_SIZE, QImage::Format_RGB888);
  texture_average.fill (QColor (0, 0, 0));
  std::vector<Eigen::Vector3d> texture_average_vec (TEXTURE_SIZE * TEXTURE_SIZE, Eigen::Vector3d::Zero ());
  std::vector<double> average_weights (TEXTURE_SIZE * TEXTURE_SIZE, 0.);

  /// Create the FBO wrapper for rendering the mesh and getting the visible triangles
  SelectionBufferGeneratorFBO sbg;
  sbg.setMesh (mesh);

  /*
   * Testing if the FBO works on your machine. It should save an image with the view of each camera.
   * If the images are all black, then probably the xml file or the mesh is incorrect.
   * Or there are problems with GL.
   */
  for (size_t cam_i = 0; cam_i < camera_params.size (); cam_i += 1)
  {
    camera_params[cam_i].print ();
    sbg.setCameraParams (camera_params[cam_i]);

    /// Example on how to get the visible triangles
    std::vector<bool> visible_tris (mesh.n_faces (), false);
    sbg.getVisibleTris (visible_tris);

    /// Get the actual rendering
    QImage render_image;
    sbg.getRenderImage (render_image);
    char str[512];
    sprintf (str, "camera_%ld.png", cam_i);
    render_image.save (str);
   }


  /*
   * Your solution should go here.
   */
  
  
    Mesh::Face_iterator f_it, f_end(mesh.faces_end());
    Mesh::Vertex_around_face_circulator vc, vc_end;
    vc = mesh.vertices(*f_it);
    vc_end = vc;
    
    for(f_it = mesh.faces_begin(); i != f_end; ++ f_it)
    {
        vc = mesh.vertices(*f_it);
       do
        {

        }while(++vc != vc_end);
    }

    for (size_t cam_i = 0; cam_i < camera_params.size(); ++cam_i)
    {
        sbg.setCameraParams(camera_params[cam_i]);
        std::vector<bool> visible_tris(mesh.n_faces(), false);
        sbg.getVisibleTris(visible_tris);
        
        for (f_it = mesh.faces_begin(); f_it != f_end; ++f_it)
        {
            if(visible_tris.at((*f_it).idx()))
            {
                for(size_t tex_i = 0; tex_i < texture_triangle_pixels.at((*f_it).idx()).size(); ++tex_i)
                {
                    Eigen::Vector2i texel = texture_triangle_pixels.at((*f_it).idx()).at(tex_i);

                    Eigen::Vector3d betas = texture_triangle_barys.at((*f_it).idx()).at(tex_i);

                    Eigen::Affine3d T = camera_params[cam_i].pose;

                    Eigen::Matrix3d K = camera_params[cam_i].getK();
                    mesh.f
                }
            }
        }
    }


  /// Save the two textures
  texture_best_view.save ("texture_best_view.png");
  /// Transfer the colors to the texture_average image and save it
  for (size_t x = 0; x < TEXTURE_SIZE; ++x)
    for (size_t y = 0; y < TEXTURE_SIZE; ++y)
    {
      int vector_index = x * TEXTURE_SIZE + y;
      QColor color (texture_average_vec[vector_index] (0),
                    texture_average_vec[vector_index] (1),
                    texture_average_vec[vector_index] (2));
      texture_average.setPixel (x, y, color.rgb ());
    }
  texture_average.save ("texture_weighted_average.png");


  return (0);
}
