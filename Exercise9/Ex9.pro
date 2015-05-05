TEMPLATE = subdirs

SUBDIRS += \
    surface_mesh \
    09-DeformationTransfer

09-DeformationTransfer.depends = surface_mesh
