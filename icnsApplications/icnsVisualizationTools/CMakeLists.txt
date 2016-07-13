# --------------------------------------------------
# Filename: CMakeLists.txt for icnsVisualizationTools
# Initial Author: Rene Werner
# Copyright (C) 2016 Department of Computational 
#   Neuroscience, UKE
# --------------------------------------------------

PROJECT( icnsVisualizationTools )

# Based on VTK:

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

LINK_LIBRARIES(
  ${VTK_LIBRARIES}
)

# Set executable:

ADD_EXECUTABLE(icnsImageSliceScreenshot icnsImageSliceScreenshot.cxx)
ADD_EXECUTABLE(icnsImageSliceScreenshot_2D icnsImageSliceScreenshot_2D.cxx)

#target_link_libraries(icnsVonEconmo ${VTK_LIBRARIES})
