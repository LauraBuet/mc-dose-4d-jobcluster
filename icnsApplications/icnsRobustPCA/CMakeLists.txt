# --------------------------------------------------
# Filename: CMakeLists.txt for icnsRobustPCA
# Initial Author: Rene Werner
# Copyright (C) 2016 Department of Computational 
#   Neuroscience, UKE
# --------------------------------------------------

PROJECT(icnsRobustPCA)

# The project depends on Armadillo:

FIND_PACKAGE(Armadillo REQUIRED)
INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})

# Set executable:

ADD_EXECUTABLE( icnsRobustPCADecomposition  
  icnsRobustPCADecomposition.cxx
  icnsMatrixFunctions.h
  icnsMatrixFunctions.cpp
  icnsRobustPCA.h  
  icnsRobustPCA.cpp
)
TARGET_LINK_LIBRARIES( icnsRobustPCADecomposition ${ARMADILLO_LIBRARIES})
ADD_EXECUTABLE( icnsRobustPCADecompositionEvaluation icnsRobustPCADecompositionEvaluation.cxx )