/*=========================================================================
 *
 *  Copyright UMC Utrecht and contributors
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkKNNGraphAlphaMutualInformationImageToImageMetric_h
#define __itkKNNGraphAlphaMutualInformationImageToImageMetric_h

/** Includes for the Superclass. */
#include "itkMultiInputImageToImageMetricBase.h"

/** Includes for the kNN trees. */
#include "itkArray.h"
#include "itkListSampleCArray.h"
#include "itkBinaryTreeBase.h"
#include "itkBinaryTreeSearchBase.h"

/** Supported trees. */
#include "itkANNkDTree.h"
#include "itkANNbdTree.h"
#include "itkANNBruteForceTree.h"

/** Supported tree searchers. */
#include "itkANNStandardTreeSearch.h"
#include "itkANNFixedRadiusTreeSearch.h"
#include "itkANNPriorityTreeSearch.h"

/** Include for the spatial derivatives. */
#include "itkArray2D.h"

#include "itkMultiThreader.h"

namespace itk
{
/**
 * \class KNNGraphAlphaMutualInformationImageToImageMetric
 *
 * \brief Computes similarity between two images to be registered.
 *
 * This metric computes the alpha-Mutual Information (aMI) between
 * two multi-channeled data sets. Said otherwise, given two sets of
 * features, the aMI between them is calculated.
 * Since for higher dimensional aMI it is infeasible to compute high
 * dimensional joint histograms, here we adopt a framework based on
 * the length of certain graphs, see Neemuchwala. Specifically, we use
 * the k-Nearest Neighbour (kNN) graph, using an implementation provided
 * by the Approximate Nearest Neighbour (ANN) software package.
 *
 * Note that the feature image are given beforehand, and that values
 * are calculated by interpolation on the transformed point. For some
 * features, it would be better (but slower) to first apply the transform
 * on the image and then recalculate the feature.
 *
 * All the technical details can be found in:\n
 * M. Staring, U.A. van der Heide, S. Klein, M.A. Viergever and J.P.W. Pluim,
 * "Registration of Cervical MRI Using Multifeature Mutual Information,"
 * IEEE Transactions on Medical Imaging, vol. 28, no. 9, pp. 1412 - 1421,
 * September 2009.
 *
 * \ingroup RegistrationMetrics
 */

template< class TFixedImage, class TMovingImage >
class KNNGraphAlphaMutualInformationImageToImageMetric :
  public MultiInputImageToImageMetricBase< TFixedImage, TMovingImage >
{
public:

  /** Standard itk. */
  typedef KNNGraphAlphaMutualInformationImageToImageMetric Self;
  typedef MultiInputImageToImageMetricBase<
    TFixedImage, TMovingImage >                             Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( KNNGraphAlphaMutualInformationImageToImageMetric,
    MultiInputImageToImageMetricBase );

  /** Typedefs from the superclass. */
  typedef typename
    Superclass::CoordinateRepresentationType CoordinateRepresentationType;
  typedef typename Superclass::MovingImageType            MovingImageType;
  typedef typename Superclass::MovingImagePixelType       MovingImagePixelType;
  typedef typename Superclass::MovingImageConstPointer    MovingImageConstPointer;
  typedef typename Superclass::FixedImageType             FixedImageType;
  typedef typename Superclass::FixedImageConstPointer     FixedImageConstPointer;
  typedef typename Superclass::FixedImageRegionType       FixedImageRegionType;
  typedef typename Superclass::TransformType              TransformType;
  typedef typename Superclass::TransformPointer           TransformPointer;
  typedef typename Superclass::InputPointType             InputPointType;
  typedef typename Superclass::OutputPointType            OutputPointType;
  typedef typename Superclass::TransformParametersType    TransformParametersType;
  typedef typename Superclass::TransformJacobianType      TransformJacobianType;
  typedef typename Superclass::InterpolatorType           InterpolatorType;
  typedef typename Superclass::InterpolatorPointer        InterpolatorPointer;
  typedef typename Superclass::RealType                   RealType;
  typedef typename Superclass::GradientPixelType          GradientPixelType;
  typedef typename Superclass::GradientImageType          GradientImageType;
  typedef typename Superclass::GradientImagePointer       GradientImagePointer;
  typedef typename Superclass::GradientImageFilterType    GradientImageFilterType;
  typedef typename Superclass::GradientImageFilterPointer GradientImageFilterPointer;
  typedef typename Superclass::FixedImageMaskType         FixedImageMaskType;
  typedef typename Superclass::FixedImageMaskPointer      FixedImageMaskPointer;
  typedef typename Superclass::MovingImageMaskType        MovingImageMaskType;
  typedef typename Superclass::MovingImageMaskPointer     MovingImageMaskPointer;
  typedef typename Superclass::MeasureType                MeasureType;
  typedef typename Superclass::DerivativeType             DerivativeType;
  typedef typename Superclass::ParametersType             ParametersType;
  typedef typename Superclass::FixedImagePixelType        FixedImagePixelType;
  typedef typename Superclass::MovingImageRegionType      MovingImageRegionType;
  typedef typename Superclass::ImageSamplerType           ImageSamplerType;
  typedef typename Superclass::ImageSamplerPointer        ImageSamplerPointer;
  typedef typename Superclass::ImageSampleContainerType   ImageSampleContainerType;
  typedef typename
    Superclass::ImageSampleContainerPointer ImageSampleContainerPointer;
  typedef typename Superclass::FixedImageLimiterType  FixedImageLimiterType;
  typedef typename Superclass::MovingImageLimiterType MovingImageLimiterType;
  typedef typename
    Superclass::FixedImageLimiterOutputType FixedImageLimiterOutputType;
  typedef typename
    Superclass::MovingImageLimiterOutputType MovingImageLimiterOutputType;
  typedef typename Superclass::NonZeroJacobianIndicesType NonZeroJacobianIndicesType;

  /** Typedef's for storing multiple inputs. */
  typedef typename Superclass::FixedImageVectorType             FixedImageVectorType;
  typedef typename Superclass::FixedImageMaskVectorType         FixedImageMaskVectorType;
  typedef typename Superclass::FixedImageRegionVectorType       FixedImageRegionVectorType;
  typedef typename Superclass::MovingImageVectorType            MovingImageVectorType;
  typedef typename Superclass::MovingImageMaskVectorType        MovingImageMaskVectorType;
  typedef typename Superclass::InterpolatorVectorType           InterpolatorVectorType;
  typedef typename Superclass::FixedImageInterpolatorVectorType FixedImageInterpolatorVectorType;

  /** The fixed image dimension. */
  itkStaticConstMacro( FixedImageDimension, unsigned int, FixedImageType::ImageDimension );
  itkStaticConstMacro( MovingImageDimension, unsigned int, MovingImageType::ImageDimension );

  /** Typedefs for the samples. */
  typedef Array< double >                           MeasurementVectorType;
  typedef typename MeasurementVectorType::ValueType MeasurementVectorValueType;
  typedef typename Statistics::ListSampleCArray<
    MeasurementVectorType, double >                   ListSampleType;
  typedef typename ListSampleType::Pointer ListSamplePointer;

  /** Typedefs for trees. */
  typedef BinaryTreeBase< ListSampleType >    BinaryKNNTreeType;
  typedef typename BinaryKNNTreeType::Pointer BinaryKNNTreePointer;
  typedef ANNkDTree< ListSampleType >         ANNkDTreeType;
  typedef ANNbdTree< ListSampleType >         ANNbdTreeType;
  typedef ANNBruteForceTree< ListSampleType > ANNBruteForceTreeType;

  /** Typedefs for tree searchers. */
  typedef BinaryTreeSearchBase< ListSampleType >     BinaryKNNTreeSearchType;
  typedef typename BinaryKNNTreeSearchType::Pointer  BinaryKNNTreeSearchPointer;
  typedef ANNStandardTreeSearch< ListSampleType >    ANNStandardTreeSearchType;
  typedef ANNFixedRadiusTreeSearch< ListSampleType > ANNFixedRadiusTreeSearchType;
  typedef ANNPriorityTreeSearch< ListSampleType >    ANNPriorityTreeSearchType;

  typedef typename BinaryKNNTreeSearchType::IndexArrayType    IndexArrayType;
  typedef typename BinaryKNNTreeSearchType::DistanceArrayType DistanceArrayType;

  typedef typename DerivativeType::ValueType        DerivativeValueType;
  typedef typename TransformJacobianType::ValueType TransformJacobianValueType;

  /** Typedefs for multi-threading. */
  typedef itk::MultiThreader                             KNNThreaderType;
  typedef typename KNNThreaderType::ThreadInfoStruct     KNNThreadInfoType;

  /**
   * *** Set trees: ***
   * Currently kd, bd, and brute force trees are supported.
   */

  /** Set ANNkDTree. */
  void SetANNkDTree( unsigned int bucketSize, std::string splittingRule );

  /** Set ANNkDTree. */
  void SetANNkDTree( unsigned int bucketSize, std::string splittingRuleFixed,
    std::string splittingRuleMoving, std::string splittingRuleJoint );

  /** Set ANNbdTree. */
  void SetANNbdTree( unsigned int bucketSize, std::string splittingRule,
    std::string shrinkingRule );

  /** Set ANNbdTree. */
  void SetANNbdTree( unsigned int bucketSize, std::string splittingRuleFixed,
    std::string splittingRuleMoving, std::string splittingRuleJoint,
    std::string shrinkingRuleFixed, std::string shrinkingRuleMoving,
    std::string shrinkingRuleJoint );

  /** Set ANNBruteForceTree. */
  void SetANNBruteForceTree( void );

  /**
   * *** Set tree searchers: ***
   * Currently standard, fixed radius, and priority tree searchers are supported.
   */

  /** Set ANNStandardTreeSearch. */
  void SetANNStandardTreeSearch( unsigned int kNearestNeighbors,
    double errorBound );

  /** Set ANNFixedRadiusTreeSearch. */
  void SetANNFixedRadiusTreeSearch( unsigned int kNearestNeighbors,
    double errorBound, double squaredRadius );

  /** Set ANNPriorityTreeSearch. */
  void SetANNPriorityTreeSearch( unsigned int kNearestNeighbors,
    double errorBound );

  /**
   * *** Standard metric stuff: ***
   */

  /** Initialize the metric. */
  virtual void Initialize( void ) throw ( ExceptionObject );

  /** Get the derivatives of the match measure. */
  void GetDerivative( const TransformParametersType & parameters,
    DerivativeType & Derivative ) const;

  /** Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType & parameters ) const;

  /** Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivativeSingleThreaded( const TransformParametersType & parameters,
	MeasureType & Value, DerivativeType & Derivative ) const;

  void GetValueAndDerivative( const TransformParametersType & parameters,
    MeasureType & Value, DerivativeType & Derivative ) const;

  /** Set alpha from alpha - mutual information. */
  itkSetClampMacro( Alpha, double, 0.0, 1.0 );

  /** Get alpha from alpha - mutual information. */
  itkGetConstReferenceMacro( Alpha, double );

  /** Avoid division by a small number. */
  itkSetClampMacro( AvoidDivisionBy, double, 0.0, 1.0 );

  /** Avoid division by a small number. */
  itkGetConstReferenceMacro( AvoidDivisionBy, double );

  /** Set the use Multi-thread. */
  itkSetMacro( UseKNNMultiThread, bool );

  /** Set the number of threads if Multi-thread is used. */
  itkSetMacro( KNNNumberOfThreads, unsigned int );

protected:

  /** Constructor. */
  KNNGraphAlphaMutualInformationImageToImageMetric();

  /** Destructor. */
  virtual ~KNNGraphAlphaMutualInformationImageToImageMetric();

  /** PrintSelf. */
  virtual void PrintSelf( std::ostream & os, Indent indent ) const;

  /** Member variables. */
  BinaryKNNTreePointer m_BinaryKNNTreeFixed;
  BinaryKNNTreePointer m_BinaryKNNTreeMoving;
  BinaryKNNTreePointer m_BinaryKNNTreeJoint;

  BinaryKNNTreeSearchPointer m_BinaryKNNTreeSearcherFixed;
  BinaryKNNTreeSearchPointer m_BinaryKNNTreeSearcherMoving;
  BinaryKNNTreeSearchPointer m_BinaryKNNTreeSearcherJoint;

  double m_Alpha;
  double m_AvoidDivisionBy;

  bool                     m_UseKNNMultiThread;
  unsigned int             m_KNNNumberOfThreads;
  KNNThreaderType::Pointer m_KNNThreader;

private:

  KNNGraphAlphaMutualInformationImageToImageMetric( const Self & ); // purposely not implemented
  void operator=( const Self & );                                   // purposely not implemented

  /** Typedef's for the computation of the derivative. */
  typedef typename Superclass::FixedImagePointType            FixedImagePointType;
  typedef typename Superclass::MovingImagePointType           MovingImagePointType;
  typedef typename Superclass::MovingImageDerivativeType      MovingImageDerivativeType;
  typedef typename Superclass::MovingImageContinuousIndexType MovingImageContinuousIndexType;
  typedef std::vector< TransformJacobianType >                TransformJacobianContainerType;
  //typedef std::vector<ParameterIndexArrayType>           TransformJacobianIndicesContainerType;
  typedef std::vector< NonZeroJacobianIndicesType > TransformJacobianIndicesContainerType;
  typedef Array2D< double >                         SpatialDerivativeType;
  typedef std::vector< SpatialDerivativeType >      SpatialDerivativeContainerType;
  typedef std::vector< IndexArrayType >             IndexArrayContainerType;
  typedef std::vector< DistanceArrayType >          DistanceArrayContainerType;

  /** This function takes the fixed image samples from the ImageSampler
   * and puts them in the listSampleFixed, together with the fixed feature
   * image samples. Also the corresponding moving image values and moving
   * feature values are computed and put into listSampleMoving. The
   * concatenation is put into listSampleJoint.
   * If desired, i.e. if doDerivative is true, then also things needed to
   * compute the derivative of the cost function to the transform parameters
   * are computed:
   * - The sparse Jacobian of the transformation (dT/dmu).
   * - The spatial derivatives of the moving (feature) images (dm/dx).
   */
  virtual void ComputeListSampleValuesAndDerivativePlusJacobian(
    const ListSamplePointer & listSampleFixed,
    const ListSamplePointer & listSampleMoving,
    const ListSamplePointer & listSampleJoint,
    const bool & doDerivative,
    TransformJacobianContainerType & jacobians,
    TransformJacobianIndicesContainerType & jacobiansIndices,
    SpatialDerivativeContainerType & spatialDerivatives ) const;

  /** This function calculates the spatial derivative of the
   * featureNr feature image at the point mappedPoint.
   * \todo move this to base class.
   */
  virtual void EvaluateMovingFeatureImageDerivatives(
    const MovingImagePointType & mappedPoint,
    SpatialDerivativeType & featureGradients ) const;

  /** This function essentially computes D1 - D2, but also takes
   * care of going from a sparse matrix (hence the indices) to a
   * full sized matrix.
   */
  virtual void UpdateDerivativeOfGammas(
    const SpatialDerivativeType & D1sparse,
    const SpatialDerivativeType & D2sparse_M,
    const SpatialDerivativeType & D2sparse_J,
    //const ParameterIndexArrayType & D1indices,
    //const ParameterIndexArrayType & D2indices_M,
    //const ParameterIndexArrayType & D2indices_J,
    const NonZeroJacobianIndicesType & D1indices,
    const NonZeroJacobianIndicesType & D2indices_M,
    const NonZeroJacobianIndicesType & D2indices_J,
    const MeasurementVectorType & diff_M,
    const MeasurementVectorType & diff_J,
    const MeasureType & distance_M,
    const MeasureType & distance_J,
    DerivativeType & dGamma_M,
    DerivativeType & dGamma_J ) const;

  /** Helper structs that multi-threads the computation of
  * the metric derivative using ITK threads.
  */
  struct KNNMultiThreaderParameterType
  {
	// To give the threads access to all members.
	KNNGraphAlphaMutualInformationImageToImageMetric * st_Metric;
	
	// Used for some parameters of accumulating derivatives
	ListSamplePointer st_FixedListSample;
	ListSamplePointer st_MovingListSample;
	ListSamplePointer st_JointListSample;
	TransformJacobianContainerType st_JacobianContainer;
	TransformJacobianIndicesContainerType st_JacobianIndicesContainer;
	SpatialDerivativeContainerType st_SpatialDerivativesContainer;

	IndexArrayContainerType st_MovingIndexArrayContainer;
	IndexArrayContainerType st_JointIndexArrayContainer;
	DistanceArrayContainerType st_FixedDistArrayContainer;
	DistanceArrayContainerType st_MovingDistArrayContainer;
	DistanceArrayContainerType st_JointDistArrayContainer;

	unsigned int st_KNearestNeighbours;
	DerivativeType st_KNNContrib;
  };
  mutable KNNMultiThreaderParameterType m_KNNThreaderMetricParameters;

  /** Most metrics will perform multi-threading by letting
  * each thread compute a part of the value and derivative.
  *
  * These parameters are initialized at every call of GetValueAndDerivative
  * in the function InitializeThreadingParameters(). Since GetValueAndDerivative
  * is const, also InitializeThreadingParameters should be const, and therefore
  * these member variables are mutable.
  */

  // test per thread struct with padding and alignment
  struct KNNGetValueAndDerivativePerThreadStruct
  {
	unsigned long st_numberOfPixelsCounted;

	ListSamplePointer st_fixedListSample;
	ListSamplePointer st_movingListSample;
	ListSamplePointer st_jointListSample;	
	TransformJacobianContainerType st_jacobianContainer;
	TransformJacobianIndicesContainerType st_jacobianIndicesContainer;
	SpatialDerivativeContainerType st_spatialDerivativesContainer;
	
	MeasureType st_sumG;
	DerivativeType st_contribution;
  };
  itkPadStruct( ITK_CACHE_LINE_ALIGNMENT, KNNGetValueAndDerivativePerThreadStruct,
	PaddedKNNGetValueAndDerivativePerThreadStruct );
  itkAlignedTypedef( ITK_CACHE_LINE_ALIGNMENT, PaddedKNNGetValueAndDerivativePerThreadStruct,
	AlignedKNNGetValueAndDerivativePerThreadStruct );
  mutable AlignedKNNGetValueAndDerivativePerThreadStruct * m_KNNGetValueAndDerivativePerThreadVariables;
  mutable ThreadIdType                                     m_KNNGetValueAndDerivativePerThreadVariablesSize;

  /** Initialize some multi-threading related parameters. */
  void InitializeKNNThreadingParameters( void ) const;

  /** KNNListSamplesAndDerivativePlusJacobian threader callback function. */
  static ITK_THREAD_RETURN_TYPE KNNListSamplesAndDerivativePlusJacobianThreaderCallback( void * arg );

  /** KNNGetValueAndDerivative threader callback function. */
  static ITK_THREAD_RETURN_TYPE KNNGetValueAndDerivativeThreaderCallback( void * arg );

  /** KNNMergeDerivative threader callback function. */
  static ITK_THREAD_RETURN_TYPE KNNMergeDerivativeThreaderCallback( void * arg );

  /** Get listsamples and derivatives plus jacobian for each thread. */
  void KNNThreadedListSamplesAndDerivativePlusJacobian( ThreadIdType threadId );

  /** Get value and derivatives for each thread. */
  void KNNThreadedGetValueAndDerivative( ThreadIdType threadId );

  /** Gather derivatives from all threads. */
  void KNNThreadedMergeDerivative( ThreadIdType threadId );

  /** Gather listsamples and derivatives plus jacobian from all threads. */
  void AfterKNNThreadedListSamplesAndDerivativePlusJacobian( void ) const;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkKNNGraphAlphaMutualInformationImageToImageMetric.hxx"
#endif

#endif // end #ifndef __itkKNNGraphAlphaMutualInformationImageToImageMetric_h
