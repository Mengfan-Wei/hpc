PROJECT(vtk-demo)

CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

SET(VTK_DIR /work/mae-weimf/lib/VTK-8.2.0/lib64/cmake/vtk-8.2)

FIND_PACKAGE(VTK COMPONENTS vtkCommonCore vtkCommonSystem
  vtkCommonMisc vtkCommonMath vtkIOCore vtkIOLegacy vtkIOXML REQUIRED)

INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
ADD_EXECUTABLE( vtkdemo pro5.cpp )
TARGET_LINK_LIBRARIES( vtkdemo ${VTK_LIBRARIES} )

# EOF
