# --------------------------------------------------
# Filename: CMakeLists.txt for icnsSurrogateBasedMotionPrediction
# Initial Author: Rene Werner
# Copyright (C) 2016 Department of Computational 
#   Neuroscience, UKE
# --------------------------------------------------

PROJECT( icnsSurrogateBasedMotionPrediction )

SET(imiLOOCV_SRCS
#  imiPCRMotionPrediction.cpp
#  imiPLSMotionPrediction.cpp
  icnsMotionPrediction.cpp
  icnsMotionPredictionLOOCV.cpp
#  imiCCAMotionPrediction.cpp
#  imiKernelMLRMotionPrediction.cpp
#  imiKernelPCA.cpp
#  imiKernel.cpp
#  imiLinearKernel.cpp
#  imiGaussianKernel.cpp
)


#SET( icnsSurrogateBasedMotionPrediction_SRCS
#  icnsMLRMotionPrediction.cpp
#)

SET( icnsGenericMLRTraining_SRCS
  icnsMLRMotionPrediction.h
  icnsMLRMotionPrediction.cpp
)

SET( icnsGenericMLR_SRCS
  icnsMLRMotionPrediction.cpp
  icnsMotionPredictionSFS.cpp
)

#SET(imiGenericPCR_SRCS
#  imiPCRMotionPrediction.cpp
#  imiPLSMotionPrediction.cpp
#  imiMotionPrediction.cpp
#  imiMotionPredictionLOOCV.cpp
#  imiCCAMotionPrediction.cpp
#  imiKernelPCA.cpp
#  imiKernel.cpp
#  imiLinearKernel.cpp
#)

#SET(imiGenericPLS_SRCS
#  imiPLSMotionPrediction.cpp
#  imiMotionPrediction.cpp
#  imiMotionPredictionLOOCV.cpp
#  imiCCAMotionPrediction.cpp
#  imiPCRMotionPrediction.cpp
#  imiKernelPCA.cpp
#  imiKernel.cpp
#  imiLinearKernel.cpp
#)

#SET(imiGenericCCA_SRCS
#  imiCCAMotionPrediction.cpp
#  imiMotionPrediction.cpp
#  imiMotionPredictionLOOCV.cpp
#  imiPLSMotionPrediction.cpp
#  imiPCRMotionPrediction.cpp
#  imiKernelPCA.cpp
#  imiKernel.cpp
#  imiLinearKernel.cpp
#)

#SET(imiKernelPCA_SRCS
#  imiKernelPCA.cpp
#  imiKernel.cpp
#  imiGaussianKernel.cpp
#  imiLinearKernel.cpp
#)

#SET(imiKernelMLR_SRCS
#  imiKernelMLRMotionPrediction.cpp
#  imiKernel.cpp
#  imiGaussianKernel.cpp
#  imiLinearKernel.cpp
#  imiMotionPredictionSFS.cpp
#)

#SET(imiMethodsMLR_SRCS
#  imiMLRMotionPrediction.cpp
#  imiMotionPredictionSFS.cpp
#)

#SET(imiKernelPCR_SRCS
#  imiKernelPCRMotionPrediction.cpp
#  ${imiKernelPCA_SRCS}
#)

#SET(imiGenerateBeltSurrogate_SRCS  
#  imiBeltSurrogate.cpp
#)

#SET(imiGenerateRangeimageSurrogate_SRCS  
#  imiRangeImageSurrogate.cpp
#)

#SET(imiGenerateDiaphragmSurrogate_SRCS  
#  imiDiaphragmSurrogate.cpp
#)

#SET(imiGenerateDiaphragmSurrogate_SRCS  
#  imiDiaphragmSurrogate.cpp
#)

ADD_EXECUTABLE(icnsMATLABMatrixCompare icnsMATLABMatrixCompare.cxx)
ADD_EXECUTABLE(icnsGenerateMatrixFromVectorField icnsGenerateMatrixFromVectorFieldsMain.cxx )
ADD_EXECUTABLE(icnsGenerateVectorFieldFromVector icnsGenerateVectorFieldFromVectorMain.cxx )

ADD_EXECUTABLE(icnsGenericMLRTraining icnsGenericMLRTrainingMain.cxx ${icnsGenericMLRTraining_SRCS})
ADD_EXECUTABLE(icnsGenericMLRPrediction icnsGenericMLRPredictionMain.cxx ${icnsGenericMLR_SRCS})
ADD_EXECUTABLE(icnsGenericMLRTrainingAndPrediction icnsGenericMLRTrainingAndPredictionMain.cxx ${icnsGenericMLR_SRCS})

#ADD_EXECUTABLE(imiGenericPCRTrainingAndPrediction imiGenericPCRTrainingAndPredictionMain.cxx ${imiGenericPCR_SRCS} ${imiLOOCV_SRCS})
#ADD_EXECUTABLE(imiGenericPLSTrainingAndPrediction imiGenericPLSTrainingAndPredictionMain.cxx ${imiGenericPLS_SRCS} ${imiLOOCV_SRCS})
#ADD_EXECUTABLE(imiGenericCCATrainingAndPrediction imiGenericCCATrainingAndPredictionMain.cxx ${imiGenericCCA_SRCS} ${imiLOOCV_SRCS})
#ADD_EXECUTABLE(imiGenericKernelMLRTrainingAndPrediction imiGenericKernelMLRTrainingAndPredictionMain.cxx ${imiKernelMLR_SRCS} ${imiLOOCV_SRCS} ${imiGenericMLR_SRCS})
#ADD_EXECUTABLE(imiGenericKernelPCRTrainingAndPrediction imiGenericKernelPCRTrainingAndPredictionMain.cxx ${imiKernelPCR_SRCS})
#ADD_EXECUTABLE(imiGenerateBeltSurrogate imiGenerateBeltSurrogateMain.cxx ${imiGenerateBeltSurrogate_SRCS})
#ADD_EXECUTABLE(imiGenerateDiaphragmSurrogate imiGenerateDiaphragmSurrogateMain.cxx ${imiGenerateDiaphragmSurrogate_SRCS})
#ADD_EXECUTABLE(imiGenerateRangeImageSurrogate imiGenerateRangeImageSurrogateMain.cxx ${imiGenerateRangeimageSurrogate_SRCS})
#ADD_EXECUTABLE(imiGenerateAirVolumeSurrogate imiGenerateAirVolumeSurrogateMain.cxx )
#ADD_EXECUTABLE(imiSurrogateSignalDifferentiation imiSurrogateSignalDifferentiation.cxx )

#ADD_EXECUTABLE(imiKernelPCAProject imiKernelPCAProject.cxx ${imiKernelPCA_SRCS})
#ADD_EXECUTABLE(imiKernelPCAComputePreImage imiKernelPCAComputePreImage.cxx ${imiKernelPCA_SRCS})
#ADD_EXECUTABLE(imiKernelRegressionMotionPredictionTest imiKernelRegressionMotionPredictionTest.cxx ${imiLOOCV_SRCS})
#ADD_EXECUTABLE(imiVnlMatrixArithmetic imiVnlMatrixArithmetic.cxx)
#ADD_EXECUTABLE(imiKernelPCATest imiKernelPCATest.cxx ${imiKernelPCA_SRCS})
#ADD_EXECUTABLE(imiKernelPCADenoising imiKernelPCADenoising.cxx ${imiKernelPCA_SRCS})


