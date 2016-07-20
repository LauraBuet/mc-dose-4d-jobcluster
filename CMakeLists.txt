# --------------------------------------------------
# Filename: CMakeLists.txt
# Initial Author: Rene Werner
# Copyright (C) 2016 Department of Computational 
#   Neuroscience, UKE
# --------------------------------------------------

# ---------- GENERAL STUFF:

# IF(UNIX)
#   set(CMAKE_C_COMPILER "gcc-4.5")
#   set(CMAKE_CXX_COMPILER "g++-4.5")
# ENDIF()

PROJECT( icnsWorking )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0012 NEW)
ENDIF(COMMAND CMAKE_POLICY)

SUBDIRS(
  icnsApplications
)

# ---------- Set the place where the lib and binary files should be built

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# ---------- Find where VTK was built or installed.

# FIND_PACKAGE(VTK)
 
# IF(VTK_FOUND)
#  INCLUDE(${VTK_USE_FILE})
# ELSE(VTK_FOUND)
#  MESSAGE(FATAL_ERROR "Cannot build without VTK. Please set VTK_DIR.")
# ENDIF(VTK_FOUND)

# SET(VTK_LIBS 
#   vtkRendering
#   vtkVolumeRendering 
#   vtkGraphics 
#   vtkHybrid 
#   vtkImaging 
#   vtkIO 
#   vtkFiltering 
#   vtkCommon
#   vtkWidgets 
# )

# ---------- Find where ITK was built or installed.

FIND_PACKAGE(ITK)

IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build without ITK. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

# ---------- Link libraries.

LINK_LIBRARIES(
  ${ITK_LIBRARIES}
 # ${VTK_LIBRARIES}
)

# SET(ITK_LIBS 
#   ITKAlgorithms
#   ITKBasicFilters 
#   ITKCommon 
#   ITKIO 
# )        
# 
# # Link ITK and VTK libraries
# LINK_LIBRARIES(
#   ${VTK_LIBS}
#   ${ITK_LIBS}
# )
# 
# LINK_DIRECTORIES(
#   ${LIBRARY_OUTPUT_PATH}
# )
