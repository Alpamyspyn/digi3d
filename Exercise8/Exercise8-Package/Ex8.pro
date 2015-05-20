TEMPLATE = subdirs

SUBDIRS += \
    surface_mesh \
    GLUTViewer \
    08-TextureBlending

GLUTViewer.depends = surface_mesh
08-TextureBlending.depends = GLUTViewer
