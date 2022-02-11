#include "imageanalysis.h"

#include <opencv2/core/cvstd.hpp>

ImageAnalysis::ImageAnalysis()
{

}

/*
 * function: MaskImage
 *  purpose: Uses the grayscale image (mask) to clear the source (image). If a pixel
 *           in mask is < threshold, then the corresponding pixel(channels) in the source
 *           image are cleared.  Note, the size of the mask and source image must be
 *           the same.
 */
cv::Mat ImageAnalysis::MaskImage(cv::Mat image, int maskThreshold) {
    cv::Mat maskedImage;
    int maskChannel = 0; // Use the IR image in the blue channel to mask

    maskedImage = image.clone();

    for (int y = 0; y < image.rows; y++ ) {
        for (int x = 0; x < image.cols; x++ ) {
            if (image.at<cv::Vec3b>(y,x)[maskChannel] <= maskThreshold) {
                maskedImage.at<cv::Vec3b>(y,x)[0] = 0;
                maskedImage.at<cv::Vec3b>(y,x)[1] = 0;
                maskedImage.at<cv::Vec3b>(y,x)[2] = 0;
            }
        } //x
    } // y

    return maskedImage;
} //MaskImage


/*
 * function: imageHistogram
 *  purpose: generates and prints to stdout a histogram of pixel values of the
 *           specified image.
 */
void ImageAnalysis::imageHistogram(cv::Mat image, long int histogram[GRAY_LEVELS],
                    long int *minPixel, long int *maxPixel, bool print)
{
    int x, y, c, i, j, pixel;
    float scale, maxCnt=0;

    *maxPixel = -1;
    *minPixel = 1000;

    printf("GRAY_LEVELS = %d\n",GRAY_LEVELS);

    for (i=0;i<GRAY_LEVELS;i++) // initialize the histogram counts
        histogram[i] = 0;

    for (y = 0; y < image.rows; y++)
        for (x = 0; x < image.cols; x++)
            for (c = 0; c < image.channels(); c++) { // only one channel in grayscale...
                pixel = image.at<cv::Vec3b>(y,x)[c];
                histogram[pixel]++;
                if (histogram[pixel]> maxCnt)
                    maxCnt = histogram[pixel];
                if (pixel < *minPixel) *minPixel = pixel;
                if (pixel > *maxPixel) *maxPixel = pixel;
            }

    if (print) {
        /* calculate scale factor for 60 column histogram */
        scale =  60.0/ maxCnt;
        printf("Image Histogram: (%f)\n", scale);
        for (i=0; i<GRAY_LEVELS; i++){
            if (histogram[i] > 0) {
                printf("[%3d]:%5ld ",i,histogram[i]);
                for (j = 0; j < histogram[i]*scale + 1; j++)
                    printf("*");
                printf("\n");
            } // if
        } // i
        puts("\n");
    }
} // imageHistogram


/*
 *
 * References:
 *  Robert Laganiere Learning OpenCV 3:: computer vision
 *  http://datahacker.rs/opencv-sobel-operator-and-image-gradient/
 */
void ImageAnalysis::EdgeDetection(cv::Mat image) {
    cv::Mat image_X, image_Y;
    cv::Mat sobelImage, image_Sobel_thresholded;
    cv::Mat image_Laplacian, image_Laplacian_thresholded;
    double max_value, min_value, max_value1, min_value1;
    double sobmin, sobmax;

    Sobel(image, image_X, CV_8UC1, 1, 0); // horizontal edge detector
    cv::imshow("Sobel image: Horizontal Edges", image_X);

    Sobel(image, image_Y, CV_8UC1, 0, 1);  // vertical edge detector.
    cv::imshow("Sobel image: Vertical Edges", image_Y);

    cv::Mat sobel = image_X + image_Y;
    cv::imshow("Sobel - L1 norm", sobel);

    /*
     * Convert the non-edges to white values and the edges to
     * dark values. This is accomplished by:
     *         sobelImage = - alpha * sobel +  255;
     */

    cv::minMaxLoc(sobel, &sobmin, &sobmax);
    sobel.convertTo(sobelImage, CV_8UC1, -255./sobmax, 255);
    cv::imshow("Edges with a sobel detector", sobelImage);

    cv::minMaxLoc(sobelImage, &min_value, &max_value);
    //threshold(sobelImage, image_Sobel_thresholded, 20, 255, cv::THRESH_BINARY);
    adaptiveThreshold(sobelImage, image_Sobel_thresholded, 255, cv::ADAPTIVE_THRESH_GAUSSIAN_C,
                      cv::THRESH_BINARY,5, 20);
    cv::imshow("Thresholded Sobel", image_Sobel_thresholded);
    cv::waitKey();


    //image_Laplacian = image_Laplacian / max_value * 255;

    // Apply low pass filtering in order to better detect edges
    // uncomment this line and the result will be much poorer.
    GaussianBlur(image, image, cv::Size(5,5), 1);
    Laplacian(image, image_Laplacian, CV_8UC1);
    cv::imshow("The Laplacian", image_Laplacian);

    minMaxLoc(image_Laplacian, &min_value1, &max_value1);
    image_Laplacian = image_Laplacian / max_value * 255;


    threshold(sobel, image_Laplacian_thresholded, 70, 220, cv::THRESH_BINARY);
    cv::imshow("Thresholded Laplacian", image_Laplacian_thresholded);
    cv::waitKey();
} // EdgeDetection



/*
 * function: CalculateNDVI
 * modified: 25-June-2020
 *
 * channel 0 = NIR channel stored in the blue channel
 * channel 1 = NIR channel duplicated in the green channel
 * channel 2 = the red channel of the Aerial image
 * NDVI = (R(NIR) - R(RED)) /  (R(NIR) + R(RED)))
 */
cv::Mat ImageAnalysis::CalculateNDVI(cv::Mat image){
    cv::Mat NDVI_Image;
    int x, y;
    int red, nir, ndvi;

    NDVI_Image = image.clone();

    for (y = 0; y < image.rows; y++)
        for (x = 0; x < image.cols; x++) {
            red = image.at<cv::Vec3b>(y,x)[2];
            nir = image.at<cv::Vec3b>(y,x)[0];
            ndvi = (nir-red)/(nir+red);

            printf("%d,%d=%d ", x, y, ndvi);

            NDVI_Image.at<cv::Vec3b>(y,x)[0] = ndvi;
            NDVI_Image.at<cv::Vec3b>(y,x)[1] = ndvi;
            NDVI_Image.at<cv::Vec3b>(y,x)[2] = ndvi;
        } // x

    return NDVI_Image;
} //CalculateNDVI


/**  @function Erosion  */
cv::Mat ImageAnalysis::Erosion(int erosion_elem, int erosion_size, cv::Mat src)
{
  int erosion_type;

  if (erosion_elem == 0 ){ erosion_type = cv::MORPH_RECT; }
  else if( erosion_elem == 1 ){ erosion_type = cv::MORPH_CROSS; }
  else if( erosion_elem == 2) { erosion_type = cv::MORPH_ELLIPSE; }

  cv::Mat element = cv::getStructuringElement( erosion_type,
                                       cv::Size( 2*erosion_size + 1, 2*erosion_size+1 ),
                                       cv::Point( erosion_size, erosion_size ) );
  cv::Mat erodedImage = src.clone();
  erode(src, erodedImage, element );

  return erodedImage;
}


/** @function Dilation */
cv::Mat ImageAnalysis::Dilation(int dilation_elem, int dilation_size, cv::Mat src)
{
  int dilation_type;

  if( dilation_elem == 0 ){ dilation_type = cv::MORPH_RECT; }
  else if( dilation_elem == 1 ){ dilation_type = cv::MORPH_CROSS; }
  else if( dilation_elem == 2) { dilation_type = cv::MORPH_ELLIPSE; }

  cv::Mat element = cv::getStructuringElement( dilation_type,
                                       cv::Size( 2*dilation_size + 1, 2*dilation_size+1 ),
                                       cv::Point( dilation_size, dilation_size ) );
  cv::Mat dilatedImage = src.clone();
  dilate( src, dilatedImage, element );

  return dilatedImage;
}
