# --------------------------------------------------
# Filename: CMakeLists.txt for icnsVonEconomo
# Initial Author: Rene Werner
# Copyright (C) 2016 Department of Computational 
#   Neuroscience, UKE
# --------------------------------------------------

PROJECT(icnsVonEconomo)

# Based on VTK:

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#LINK_LIBRARIES(
#  ${VTK_LIBRARIES}
#)

# Set executable:
# VE3D core executables:

ADD_EXECUTABLE(icnsVonEconomo MACOSX_BUNDLE icnsVonEconomo.cxx)
ADD_EXECUTABLE(icnsVE3D_defineLandmarks MACOSX_BUNDLE icnsVE3D_defineLandmarks.cxx vtkTexturingHelper.h vtkTexturingHelper.cpp)

# Surface registration functionality:

ADD_EXECUTABLE(icnsManualSurfaceRegistration MACOSX_BUNDLE icnsManualSurfaceRegistration.cxx vtkTexturingHelper.h vtkTexturingHelper.cpp)
ADD_EXECUTABLE(icnsLandmarkBasedSurfaceRegistration MACOSX_BUNDLE icnsLandmarkBasedSurfaceRegistration.cxx vtkTexturingHelper.h vtkTexturingHelper.cpp)
ADD_EXECUTABLE(icnsICPBasedSurfaceRegistration MACOSX_BUNDLE icnsICPBasedSurfaceRegistration.cxx vtkTexturingHelper.h vtkTexturingHelper.cpp)

# Surface model viewer functionality:

ADD_EXECUTABLE(icnsViewVTKPolyData MACOSX_BUNDLE icnsViewVTKPolyData.cxx vtkTexturingHelper.h vtkTexturingHelper.cpp)

# Surface model manipulation functionality:

ADD_EXECUTABLE(icnsFlipVTKPolyData MACOSX_BUNDLE icnsFlipVTKPolyData.cxx vtkTexturingHelper.h vtkTexturingHelper.cpp vtkOBJWriter.h vtkOBJWriter.cxx)
 
# Define libraries to link:

if(VTK_LIBRARIES)
  target_link_libraries(icnsVonEconomo ${VTK_LIBRARIES})
  target_link_libraries(icnsVE3D_defineLandmarks ${VTK_LIBRARIES})
  target_link_libraries(icnsManualSurfaceRegistration ${VTK_LIBRARIES})
  target_link_libraries(icnsLandmarkBasedSurfaceRegistration ${VTK_LIBRARIES})
  target_link_libraries(icnsICPBasedSurfaceRegistration ${VTK_LIBRARIES})
  target_link_libraries(icnsViewVTKPolyData ${VTK_LIBRARIES})
  target_link_libraries(icnsFlipVTKPolyData ${VTK_LIBRARIES})
else()
  target_link_libraries(icnsVonEconomo vtkHybrid vtkWidgets)
  target_link_libraries(icnsVE3D_defineLandmarks vtkHybrid vtkWidgets)
  target_link_libraries(icnsICPBasedSurfaceRegistration vtkHybrid vtkWidgets)
  target_link_libraries(icnsLandmarkBasedSurfaceRegistration vtkHybrid vtkWidgets)
  target_link_libraries(icnsViewVTKPolyData vtkHybrid vtkWidgets)
  target_link_libraries(icnsFlipVTKPolyData vtkHybrid vtkWidgets)
endif()
