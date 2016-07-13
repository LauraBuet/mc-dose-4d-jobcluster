# --------------------------------------------------
# Filename: CMakeLists.txt for icnsImageTools
# Initial Author: Rene Werner
# Copyright (C) 2016 Department of Computational 
#   Neuroscience, UKE
# --------------------------------------------------

PROJECT( icnsImageTools )

# Set executable:

ADD_EXECUTABLE( icnsAffineImageRegistration icnsAffineImageRegistration.cxx )
ADD_EXECUTABLE( icnsAffineImageTransformation icnsAffineImageTransformation.cxx )
ADD_EXECUTABLE( icnsImageCropping icnsImageCropping.cxx )
ADD_EXECUTABLE( icnsImageHistogramMatching icnsImageHistogramMatching.cxx )
ADD_EXECUTABLE( icnsImageResampling icnsImageResampling.cxx )
ADD_EXECUTABLE( icnsImageToImageDistanceComputation icnsImageToImageDistanceComputation.cxx )
ADD_EXECUTABLE( icnsImageThresholding icnsImageThresholding.cxx )
ADD_EXECUTABLE( icnsPythagoreanImageGeneration icnsPythagoreanImageGeneration.cxx )
ADD_EXECUTABLE( icnsSegmentationAccuracyEvaluation icnsSegmentationAccuracyEvaluation.cxx )