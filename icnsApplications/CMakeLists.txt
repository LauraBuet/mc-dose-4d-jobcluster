# --------------------------------------------------
# Filename: CMakeLists.txt.h
# Initial Author: Rene Werner
# Copyright (C) 2016 Department of Computational
#   Neuroscience, UKE
# --------------------------------------------------
cmake_minimum_required(VERSION 3.0)
SUBDIRS(
  icnsImageTools
  icnsRobustPCA
  icnsSurrogateBasedMotionPrediction
  icnsVariationalRegistration
  icnsVectorFieldTools
 # icnsVisualizationTools
 # icnsVonEconomo3D
)
