TEMPLATE = subdirs

SUBDIRS += \
    surface_mesh \
    GLUTViewer \
    06-RigidRegistration

GLUTViewer.depends = surface_mesh
06-RigidRegistration.depends = GLUTViewer
