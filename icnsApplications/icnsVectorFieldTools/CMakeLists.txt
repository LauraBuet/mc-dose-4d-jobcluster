# --------------------------------------------------
# Filename: CMakeLists.txt for icnsVectorFieldTools
# Initial Author: Rene Werner
# Copyright (C) 2016 Department of Computational 
#   Neuroscience, UKE
# --------------------------------------------------

PROJECT( icnsVectorFieldTools )

# Set executable:

ADD_EXECUTABLE( icnsComposeVectorFields icnsComposeVectorFields.cxx )
ADD_EXECUTABLE( icnsValidateLandmarks icnsValidateLandmarks.cxx )
ADD_EXECUTABLE( icnsVectorFieldDeformationMetricComputation icnsVectorFieldDeformationMetricComputation.cxx )
ADD_EXECUTABLE( icnsVectorFieldStatistics icnsVectorFieldStatistics.cxx )