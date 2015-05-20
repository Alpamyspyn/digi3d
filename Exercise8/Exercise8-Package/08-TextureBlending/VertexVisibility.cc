#include "VertexVisbility.hh"


void
SelectionBufferGeneratorFBO::setCameraParams (const CameraParams &params)
{
  proj_matrix_ = Eigen::Matrix4f::Identity ();
  float zFar = 1000., zNear = 0.1;
  proj_matrix_ (0, 0) = 2. * params.fx / static_cast<float> (params.width);
  proj_matrix_ (0, 2) = (static_cast<float> (params.width) - 2. * params.cx) / static_cast<float> (params.width);
  proj_matrix_ (1, 1) = 2. * params.fy / static_cast<float> (params.height);
  proj_matrix_ (1, 2) = (-static_cast<float> (params.height) + 2. * params.cy) / static_cast<float> (params.height);
  proj_matrix_ (2, 2) = (-zFar - zNear) / (zFar - zNear);
  proj_matrix_ (2, 3) = -2. * zFar * zNear / (zFar - zNear);
  proj_matrix_ (3, 2) = -1.;
  proj_matrix_ (3, 3) = 0;

  model_matrix_ = params.pose.inverse ().matrix ().cast<float> ();
  model_matrix_.row (1) *= -1.;
  model_matrix_.row (2) *= -1.;

  width_ = params.width;
  height_ = params.height;

  std::cout << "set camera params\n";
}


void
SelectionBufferGeneratorFBO::draw ()
{
  std::cout << "drawing ...\n";
  GL::glCheckErrors ();

  /// Do the actual drawing
  glDisable (GL_TEXTURE_2D);
  glDisable (GL_BLEND);
  glDisable (GL_COLOR_MATERIAL);
  glClearColor (0., 0., 0., 0.);
  glClearDepth (1.);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);
  glDisable (GL_LIGHT0);
  glDisable (GL_LIGHTING);

  GL::glCheckErrors ();

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  glMultMatrixf (proj_matrix_.data ());

  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity ();
  glMultMatrixf (model_matrix_.data ());

  GL::glCheckErrors ();
  glShadeModel (GL_FLAT);

  GL::glCheckErrors ();

  // if (!wireframe_)
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
  // else
    // glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

  glBegin (GL_TRIANGLES);
  for (Mesh::Face_iterator f_it = mesh_.faces_begin (); f_it != mesh_.faces_end (); ++f_it)
  {
    unsigned char r, g, b;
    size_t temp = (*f_it).idx () + 1; /// to avoid confusion with the black background
    r = temp / (256 * 256);
    g = (temp % (256 * 256)) / 256;
    b = temp % 256;

    glColor3ub (r, g, b);
    for (Mesh::Vertex_around_face_circulator v_it = Mesh::Vertex_around_face_circulator (&mesh_, *f_it).begin ();
         v_it != Mesh::Vertex_around_face_circulator (&mesh_, *f_it).end (); ++v_it)
    {
      glVertex3f (mesh_.position (*v_it)[0],
                  mesh_.position (*v_it)[1],
                  mesh_.position (*v_it)[2]);
    }
  }
  glEnd ();

  GL::glCheckErrors ();

  glColor3ub (255, 255, 255);
}


GLuint SelectionBufferGeneratorFBO::fbo_id_ = 0;
GLuint SelectionBufferGeneratorFBO::rb_id_ = 0;
GLuint SelectionBufferGeneratorFBO::tex_id_ = 0;
int SelectionBufferGeneratorFBO::win_id_ = 0;
void
SelectionBufferGeneratorFBO::getProjectedIds (Eigen::MatrixXi &buffer)
{
  QImage img_frame;
  std::cerr << "before getrenderimage\n";
  getRenderImage (img_frame);
  std::cerr << "after getrenderimage\n";

  /// Now copy the cv::Mat into the Eigen matrix
  buffer = Eigen::MatrixXi::Ones (height_, width_) * (-1);
  for (int x = 0; x < width_; ++x)
    for (int y = 0; y < height_; ++y)
    {
      QRgb pix = img_frame.pixel (x, y);
      int tri_id = static_cast<int> (qRed (pix)) * 256 * 256 +
                   static_cast<int> (qGreen (pix)) * 256 +
                   static_cast<int> (qBlue (pix));

      if (tri_id != 0)
        buffer (y, x) = tri_id - 1;
    }
}


void
display ()
{}


void
SelectionBufferGeneratorFBO::getRenderImage (QImage &render)
{
  /// Check if GLUT was initialized already (by ourselves or by VTK)
  /// Don't do it a second time
  if (glutGetWindow () == 0)
  {
    std::cerr << "Creating glut window\n";
    int argc = 0;
    char** argv = NULL;
    glutInit (&argc, argv);
    glutInitWindowSize (320, 320);
    win_id_ = glutCreateWindow ("glut dummy window render_with_texture");
#ifndef __APPLE__
    glewInit ();
#endif
    glutDisplayFunc (display);
  }

  win_id_ = glutGetWindow ();
  glutSetWindow (win_id_);
  GL::glCheckErrors ();

  if (rb_id_ == 0)
  {
    glGenFramebuffers (1, &fbo_id_);
    glGenRenderbuffers (1, &rb_id_);
    glGenTextures (1, &tex_id_);
    glBindFramebuffer (GL_FRAMEBUFFER, fbo_id_);

    glBindTexture (GL_TEXTURE_2D, tex_id_);
    glTexImage2D (GL_TEXTURE_2D, 0, GL_RGB, width_, height_, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glFramebufferTexture2D (GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_id_, 0);

    glBindRenderbuffer (GL_RENDERBUFFER, rb_id_);
    glRenderbufferStorage (GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width_, height_);

    glFramebufferRenderbuffer (GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rb_id_);

    GL::glCheckErrors ();

    GLenum status = glCheckFramebufferStatus (GL_FRAMEBUFFER);
    if (status != GL_FRAMEBUFFER_COMPLETE)
      std::cerr << "Error with our framebuffer.\n";
  }

  GL::glCheckErrors ();
  glViewport (0, 0, width_, height_);

  GL::glCheckErrors ();
  draw ();
  std::cerr << "drawn\n";

  GL::glCheckErrors ();

  /// Get the data from the texture into a cv::Mat
  glReadBuffer (GL_COLOR_ATTACHMENT0);

  GL::glCheckErrors ();

  std::cerr << "getting the data\n";
  render = QImage (width_, height_, QImage::Format_RGB888);
  render.fill (QColor (0, 0, 0));
  glReadPixels (0, 0, width_, height_, GL_RGB, GL_UNSIGNED_BYTE, render.bits ());
  std::cerr << "got the data\n";

  /// Check how many non-zero pixels we have
  int count = 0;
  for (size_t x = 0; x < render.width (); ++x)
    for (size_t y = 0; y < render.height (); ++y)
    {
      QRgb pix = render.pixel (x, y);
      if (qRed (pix) != 0 ||
          qGreen (pix) != 0 ||
          qBlue (pix) != 0)
        count ++;
    }
  std::cerr << "non-zero pixels: " << count << std::endl;


  GL::glCheckErrors ();

  glDeleteTextures (1, &tex_id_);
  glDeleteRenderbuffers (1, &rb_id_);
  glBindFramebuffer (GL_FRAMEBUFFER, 0);
  glDeleteFramebuffers (1, &fbo_id_);
  tex_id_ = rb_id_ = fbo_id_ = 0;
  GL::glCheckErrors ();
}


void
SelectionBufferGeneratorFBO::getVisibleTris (std::vector<bool> &visible_tris)
{
  Eigen::MatrixXi projected_ids;
  getProjectedIds (projected_ids);
  for (int x = 0; x < projected_ids.rows (); ++x)
    for (int y = 0; y < projected_ids.cols (); ++y)
      if (projected_ids (x, y) != -1)
        visible_tris[projected_ids (x, y)] = true;
}



