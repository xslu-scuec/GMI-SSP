
ADD_ELXCOMPONENT( KNNGraphAlphaMutualInformationMetric
 elxKNNGraphAlphaMutualInformationMetric.h
 elxKNNGraphAlphaMutualInformationMetric.hxx
 elxKNNGraphAlphaMutualInformationMetric.cxx
 itkKNNGraphAlphaMutualInformationImageToImageMetric.h
 itkKNNGraphAlphaMutualInformationImageToImageMetric.hxx
)

if( USE_KNNGraphAlphaMutualInformationMetric )
  target_link_libraries( KNNGraphAlphaMutualInformationMetric KNNlib ANNlib )
endif()

