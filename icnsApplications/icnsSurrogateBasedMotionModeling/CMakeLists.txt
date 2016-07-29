##############################################################################
#
# Filename: CMakeLists.txt
# Initial Author: Rene Werner
# Copyright (C) 2012 Institute of Medical Informatics, University of Luebeck
#
##############################################################################

PROJECT(icnsSurrogateBasedMotionPrediction)

SET(imiLOOCV_SRCS
  imiPCRMotionPrediction.cpp
  imiPLSMotionPrediction.cpp
  imiMotionPrediction.cpp
  imiMotionPredictionLOOCV.cpp
  imiCCAMotionPrediction.cpp
  imiKernelMLRMotionPrediction.cpp
  imiKernelPCA.cpp
  imiKernel.cpp
  imiLinearKernel.cpp
  imiGaussianKernel.cpp
)

SET(imiSurrogateBasedMotionPrediction_SRCS
  imiMLRMotionPrediction.cpp
)

SET(imiGenericMLR_SRCS
  imiMLRMotionPrediction.cpp
  imiMotionPredictionSFS.cpp
)

SET(imiGenericPCR_SRCS
  imiPCRMotionPrediction.cpp
  imiPLSMotionPrediction.cpp
  imiMotionPrediction.cpp
  imiMotionPredictionLOOCV.cpp
  imiCCAMotionPrediction.cpp
  imiKernelPCA.cpp
  imiKernel.cpp
  imiLinearKernel.cpp
)

SET(imiGenericPLS_SRCS
  imiPLSMotionPrediction.cpp
  imiMotionPrediction.cpp
  imiMotionPredictionLOOCV.cpp
  imiCCAMotionPrediction.cpp
  imiPCRMotionPrediction.cpp
  imiKernelPCA.cpp
  imiKernel.cpp
  imiLinearKernel.cpp
)

SET(imiGenericCCA_SRCS
  imiCCAMotionPrediction.cpp
  imiMotionPrediction.cpp
  imiMotionPredictionLOOCV.cpp
  imiPLSMotionPrediction.cpp
  imiPCRMotionPrediction.cpp
  imiKernelPCA.cpp
  imiKernel.cpp
  imiLinearKernel.cpp
)

SET(imiKernelPCA_SRCS
  imiKernelPCA.cpp
  imiKernel.cpp
  imiGaussianKernel.cpp
  imiLinearKernel.cpp
)

SET(imiKernelMLR_SRCS
  imiKernelMLRMotionPrediction.cpp
  imiKernel.cpp
  imiGaussianKernel.cpp
  imiLinearKernel.cpp
  imiMotionPredictionSFS.cpp
)

SET(imiMethodsMLR_SRCS
  imiMLRMotionPrediction.cpp
  imiMotionPredictionSFS.cpp
)

SET(imiKernelPCR_SRCS
  imiKernelPCRMotionPrediction.cpp
  ${imiKernelPCA_SRCS}
)

SET(imiGenerateBeltSurrogate_SRCS  
  imiBeltSurrogate.cpp
)

SET(imiGenerateRangeimageSurrogate_SRCS  
  imiRangeImageSurrogate.cpp
)

SET(imiGenerateDiaphragmSurrogate_SRCS  
  imiDiaphragmSurrogate.cpp
)

SET(imiGenerateDiaphragmSurrogate_SRCS  
  imiDiaphragmSurrogate.cpp
)

#ADD_EXECUTABLE(imiMATLABMatrixCompare imiMATLABMatrixCompare.cxx)
ADD_EXECUTABLE(imiGenericMLRTraining imiGenericMLRTrainingMain.cxx ${imiGenericMLR_SRCS})
#ADD_EXECUTABLE(imiGenericMLRPrediction imiGenericMLRPredictionMain.cxx ${imiGenericMLR_SRCS})
#ADD_EXECUTABLE(imiGenericMLRTrainingAndPrediction imiGenericMLRTrainingAndPredictionMain.cxx ${imiGenericMLR_SRCS})
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
#ADD_EXECUTABLE(imiGenerateMatrixFromVectorField imiGenerateMatrixFromVectorFieldsMain.cxx )
#ADD_EXECUTABLE(imiGenerateVectorFieldFromVector imiGenerateVectorFieldFromVectorMain.cxx )

#ADD_EXECUTABLE(imiKernelPCAProject imiKernelPCAProject.cxx ${imiKernelPCA_SRCS})
#ADD_EXECUTABLE(imiKernelPCAComputePreImage imiKernelPCAComputePreImage.cxx ${imiKernelPCA_SRCS})
#ADD_EXECUTABLE(imiKernelRegressionMotionPredictionTest imiKernelRegressionMotionPredictionTest.cxx ${imiLOOCV_SRCS})
#ADD_EXECUTABLE(imiVnlMatrixArithmetic imiVnlMatrixArithmetic.cxx)
#ADD_EXECUTABLE(imiKernelPCATest imiKernelPCATest.cxx ${imiKernelPCA_SRCS})
#ADD_EXECUTABLE(imiKernelPCADenoising imiKernelPCADenoising.cxx ${imiKernelPCA_SRCS})


