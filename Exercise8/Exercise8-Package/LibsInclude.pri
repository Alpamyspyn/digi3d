INCLUDEPATH += $${_PRO_FILE_PWD_}/../
INCLUDEPATH += $${_PRO_FILE_PWD_}/../external
INCLUDEPATH += $${_PRO_FILE_PWD_}/../external/eigen

#Generate binary executables at the root of the build folder
DESTDIR = $${OUT_PWD}/../

win32{
    INCLUDEPATH *= $${_PRO_FILE_PWD_}/../freeglut/include
    LIBS += "$${OUT_PWD}/../GLUTViewer/GLUTViewer.lib"
    LIBS += "$${OUT_PWD}/../surface_mesh/surface_mesh.lib"

    # For glut, link against the lib file and copy the dll file to the same folder as the executable
    contains(QMAKE_TARGET.arch, x86_64){
        message("win64 build GLUT")
        LIBS += "$${_PRO_FILE_PWD_}/../freeglut/lib/win64/freeglut.lib"
        QMAKE_POST_LINK += copy $$shell_path($${_PRO_FILE_PWD_}/../freeglut/bin/win64/freeglut.dll) $$shell_path($${OUT_PWD}/../) /Y
    }

    # For glew, link against the lib file and copy the dll file to the same folder as the executable
    INCLUDEPATH *= $${_PRO_FILE_PWD_}/../glew/include
    contains(QMAKE_TARGET.arch, x86_64){
        message("win64 build GLEW")
        LIBS += "$${_PRO_FILE_PWD_}/../glew/lib/Release/x64/glew32.lib"
        QMAKE_POST_LINK += & copy $$shell_path($${_PRO_FILE_PWD_}/../glew/bin/Release/x64/glew32.dll) $$shell_path($${OUT_PWD}/../) /Y
    }
}

unix{
    INCLUDEPATH *= /usr/include/eigen3
    INCLUDEPATH *= /usr/local/include/eigen3
    INCLUDEPATH *= /usr/include
    INCLUDEPATH *= /usr/local/include

    LIBS += -L$${OUT_PWD}/../GLUTViewer -L$${OUT_PWD}/../surface_mesh -lGLUTViewer -lsurface_mesh

    INCLUDEPATH += $${_PRO_FILE_PWD_}/../glew/include
    LIBS += -L$${_PRO_FILE_PWD_}/../glew/build/lib/

    macx{
        LIBS += -framework GLUT -framework OpenGL -lGLEW
    }
    else{
        LIBS += -lGL -lGLU -lglut -lGLEW
    }
}
