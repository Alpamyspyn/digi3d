INCLUDEPATH += $${_PRO_FILE_PWD_}/../
INCLUDEPATH += $${_PRO_FILE_PWD_}/../external
INCLUDEPATH += $${_PRO_FILE_PWD_}/../external/eigen

#Generate binary executables at the root of the build folder
DESTDIR = $${OUT_PWD}/../

win32{
    LIBS += "$${OUT_PWD}/../surface_mesh/surface_mesh.lib"

    
}

unix{
    INCLUDEPATH *= /usr/include/eigen3
    INCLUDEPATH *= /usr/local/include/eigen3
    INCLUDEPATH *= /usr/include
    INCLUDEPATH *= /usr/local/include


    LIBS += -L$${OUT_PWD}/../surface_mesh -lsurface_mesh
}
