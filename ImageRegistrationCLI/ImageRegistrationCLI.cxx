#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkCastImageFilter.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNormalizeImageFilter.h"

#include "itkCenteredTransformInitializer.h"
#include "itkEuler2DTransform.h"
#include "itkExhaustiveOptimizer.h"

#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"

#include "itkTimeProbesCollectorBase.h"

#include "itkPluginUtilities.h"
#include "ImageRegistrationCLICLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
typedef double FixedPixelType;
typedef double MovingPixelType;

typedef itk::Image<FixedPixelType, 2>  FixedInputImageType;
typedef itk::Image<MovingPixelType, 2> MovingInputImageType;

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef itk::ImageFileReader<FixedInputImageType>  FixedReaderType;
  typedef itk::ImageFileReader<MovingInputImageType> MovingReaderType;

   FixedReaderType::Pointer fixedReader   = FixedReaderType::New();
   MovingReaderType::Pointer movingReader = MovingReaderType::New();

  fixedReader->SetFileName( fixedImage.c_str() );
  movingReader->SetFileName( movingImage.c_str() );

  try
   {
      fixedReader->Update();
      movingReader->Update();
   }
  catch ( itk::ExceptionObject & e )
   {
      std::cerr << "exception in file reader " << std::endl;
      std::cerr << e.GetDescription() << std::endl;
      std::cerr << e.GetLocation() << std::endl;
      return EXIT_FAILURE;
   }

  FixedInputImageType::Pointer  fixImage = fixedReader->GetOutput();
  MovingInputImageType::Pointer movImage = movingReader->GetOutput();

  itk::TimeProbesCollectorBase firstTimer;
  firstTimer.Start("Image Resampling");
  // Normalize the images
  typedef itk::NormalizeImageFilter<FixedInputImageType, FixedInputImageType> NormalizeFilterType;
 
  NormalizeFilterType::Pointer fixedNormalizer = NormalizeFilterType::New();
  NormalizeFilterType::Pointer movingNormalizer = NormalizeFilterType::New();
 
  fixedNormalizer->SetInput(  fixImage );
  movingNormalizer->SetInput( movImage );

  fixedNormalizer->Update();
  movingNormalizer->Update();

  // Resample moving image to have the same size as the fixed image

  typedef itk::IdentityTransform<double, 2> IdentityTransformType;
   IdentityTransformType::Pointer identityTransform = IdentityTransformType::New(); 

  typedef itk::ResampleImageFilter< MovingInputImageType, MovingInputImageType > ResampleFilterType;
   ResampleFilterType::Pointer resampler = ResampleFilterType::New();
 
  identityTransform->SetIdentity();

  unsigned int fixedImageSize = fixedNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()[0] * fixedNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
  unsigned int movingImageSize = movingNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()[0] * movingNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()[1];

  std::cout<<"Fixed image size: "<<fixedNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()<<std::endl;
  std::cout<<"Moving image size: "<<movingNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()<<std::endl;

  double OutputSpacing[2];

  if(movingImageSize < fixedImageSize)
  {
	  double oldWidth = movingNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
	  double newWidth = fixedNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()[0];  

	  double oldHeight = movingNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
      double newHeight = fixedNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()[1]; 

	  OutputSpacing[0] = fixedNormalizer->GetOutput()->GetSpacing()[0] * (double) oldWidth / (double) newWidth;
	  OutputSpacing[1] = fixedNormalizer->GetOutput()->GetSpacing()[1] * (double) oldHeight / (double) newHeight;

      resampler->SetInput( movingNormalizer->GetOutput() );
	  resampler->SetTransform( identityTransform );
	  resampler->SetSize( fixedNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize() );
	  resampler->SetOutputOrigin(  fixedNormalizer->GetOutput()->GetOrigin() );
	  resampler->SetOutputSpacing( OutputSpacing );
	  resampler->SetOutputDirection( fixImage->GetDirection() );
	  resampler->UpdateLargestPossibleRegion();
	  resampler->SetDefaultPixelValue(100);  
  }
  else
  {
	  double oldWidth  = fixedNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()[0];  
	  double newWidth  = movingNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
    
	  double oldHeight = fixedNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()[1]; 
	  double newHeight = movingNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize()[1];

	  OutputSpacing[0] = movingNormalizer->GetOutput()->GetSpacing()[0] * (double) oldWidth / (double) newWidth;
	  OutputSpacing[1] = movingNormalizer->GetOutput()->GetSpacing()[1] * (double) oldHeight / (double) newHeight;
 
	  resampler->SetInput( fixedNormalizer->GetOutput() );
	  resampler->SetTransform( identityTransform );
	  resampler->SetSize( movingNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize() );
	  resampler->SetOutputOrigin(  movingNormalizer->GetOutput()->GetOrigin() );
	  resampler->SetOutputSpacing( OutputSpacing );
	  resampler->SetOutputDirection( movImage->GetDirection() );
	  resampler->UpdateLargestPossibleRegion();
	  resampler->SetDefaultPixelValue(100);  
  }

  resampler->Update();
  firstTimer.Stop("Image Resampling");
  // Register the two images

  typedef itk::Euler2DTransform<double>                                                       	 TransformType;
  typedef itk::ExhaustiveOptimizer                                                   	 		 OptimizerType; 
  typedef itk::MattesMutualInformationImageToImageMetric< FixedInputImageType, MovingInputImageType >    MetricType;
  typedef itk::LinearInterpolateImageFunction< MovingInputImageType, double >	         		 InterpolatorType;
  typedef itk::CenteredTransformInitializer< TransformType, FixedInputImageType,  MovingInputImageType > TransformInitializerType;

  typedef itk::MultiResolutionImageRegistrationMethod< FixedInputImageType, MovingInputImageType >	 RegistrationType;
  typedef itk::MultiResolutionPyramidImageFilter< FixedInputImageType, FixedInputImageType >   		 FixedImagePyramidType;
  typedef itk::MultiResolutionPyramidImageFilter< MovingInputImageType, MovingInputImageType >   	 MovingImagePyramidType;

   TransformType::Pointer            transform    = TransformType::New();
   MetricType::Pointer               metric       = MetricType::New();
   OptimizerType::Pointer            optimizer    = OptimizerType::New();
   InterpolatorType::Pointer         interpolator = InterpolatorType::New();
   RegistrationType::Pointer         registration = RegistrationType::New();
   TransformInitializerType::Pointer initializer  = TransformInitializerType::New();

  itk::TimeProbesCollectorBase timer;
  itk::TimeProbesCollectorBase singleTimer;
  timer.Start("Full Process");

  // Set the number of samples (radius) in each dimension, with a default step size of 1.0

  unsigned int angles = anglesNumber;

  OptimizerType::StepsType steps( transform->GetNumberOfParameters() );
  steps[0] = int(angles);
  steps[1] = 0;
  steps[2] = 0;

  optimizer->SetNumberOfSteps( steps );

  // Utilize the scale to set the step size for each dimension
  OptimizerType::ScalesType scales( transform->GetNumberOfParameters() );
  scales[0] = 2.0 * vnl_math::pi / angles;
  scales[1] = 1.0;
  scales[2] = 1.0; 
  optimizer->SetScales( scales );
  optimizer->SetStepLength(1.0);
  std::cout<<optimizer->GetStepLength()<<std::endl;
  //
  //  Initialize a rigid transform by using Image Intensity Moments
  //
  singleTimer.Start("Initializer");
  initializer->SetTransform(   transform );
  if(movingImageSize < fixedImageSize)
  {
	  initializer->SetFixedImage(  fixedNormalizer->GetOutput()  );
	  initializer->SetMovingImage( resampler->GetOutput() );
  }
  else
  {
	  initializer->SetFixedImage(  resampler->GetOutput()  );
	  initializer->SetMovingImage( movingNormalizer ->GetOutput() );
  }
  initializer->GeometryOn();
  initializer->InitializeTransform();
  singleTimer.Stop("Initializer");

  metric->SetNumberOfHistogramBins( 128 );
  const unsigned int numberOfSamples = int (fixedNormalizer->GetOutput()->GetBufferedRegion().GetNumberOfPixels() * 60.0 / 100.0 );
  metric->SetNumberOfSpatialSamples( numberOfSamples );
  metric->SetUseExplicitPDFDerivatives(0);
  metric->ReinitializeSeed( 76926294 );

  // Initialize registration

  singleTimer.Start("Registration");
  FixedImagePyramidType::Pointer fixedImagePyramid   = FixedImagePyramidType::New();
  MovingImagePyramidType::Pointer movingImagePyramid = MovingImagePyramidType::New();

  registration->SetMetric(      metric    );
  registration->SetOptimizer(   optimizer );
  registration->SetInterpolator(  interpolator  );
  if(movingImageSize < fixedImageSize)
  {
	  registration->SetFixedImage(  fixedNormalizer->GetOutput()  );
 	  registration->SetFixedImageRegion( fixedNormalizer->GetOutput()->GetBufferedRegion() );
	  registration->SetMovingImage( resampler->GetOutput()  );
  }
  else
  {
	  registration->SetFixedImage(  resampler->GetOutput()  );
 	  registration->SetFixedImageRegion( resampler->GetOutput()->GetBufferedRegion() );
	  registration->SetMovingImage( movingNormalizer->GetOutput()  );
  }

  registration->SetFixedImagePyramid(  fixedImagePyramid  );
  registration->SetMovingImagePyramid( movingImagePyramid  );
  registration->SetInitialTransformParameters( transform->GetParameters() );
  registration->SetTransform( transform );
  registration->SetNumberOfLevels(2);

  std::cout << "resampler dim  = " << resampler->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;
  //std::cout << "Initial Metric value  = " << metric->GetValue( registration->GetInitialTransformParameters() ) << std::endl;

  try
  {
    registration->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }


  singleTimer.Stop("Registration");
  //  The value of the image metric corresponding to the best set of parameters
  //  can be obtained with the GetMinimumMetricValuePosition method of the optimizer.

  typedef  TransformType::ParametersType ParametersType;
  ParametersType bestParameters = optimizer->GetMinimumMetricValuePosition();
  std::cout << "Final parameters: " << bestParameters << std::endl;

  // The best value of the image metric can be obtained with the  GetMinimumMetricValue
  // method of the optimizer
  const float bestValue = optimizer->GetMinimumMetricValue();

  // The best rotation angle is converted from rad into degrees 
  const float bestAngleInDegrees = bestParameters[0] * 180.0 / vnl_math::pi;

  // Print out results
  //
  optimizer->Print(std::cout);
  std::cout << "Result = " << std::endl;
  std::cout << " Metric value  = "    << bestValue << std::endl;
  std::cout << " Angle (rad) =  "     << bestParameters[0]  << std::endl;
  std::cout << " Angle (degrees) =  " << bestAngleInDegrees  << std::endl;
  std::cout << " Rotation Center = "  << transform->GetCenter() << std::endl;
  std::cout << registration->GetOptimizer()->GetStopConditionDescription() << std::endl;

  registrationAngle = (int)bestAngleInDegrees;
  std::ofstream rts;
  rts.open(returnParameterFile.c_str() );
  rts << "registrationAngle = " << registrationAngle << std::endl;
  rts.close();

  timer.Stop("Full Process");
  timer.Report( std::cout );
  singleTimer.Report( std::cout );
  firstTimer.Report( std::cout );
  return EXIT_SUCCESS;
}
