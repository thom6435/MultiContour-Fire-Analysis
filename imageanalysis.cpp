#include "imageanalysis.h"
#include <opencv2/core/cvstd.hpp>

#include <math.h>
#include <stdio.h>
#include <vector>
#include <malloc.h>
#include <algorithm>


ImageAnalysis::ImageAnalysis()
{

}

ImageAnalysis::ImageAnalysis(int width, int height) {
    this->imageWidth = width;
    this->imageHeight = height;
}

void ImageAnalysis::setImageWidth(int width) {
    this->imageWidth = width;
}

void ImageAnalysis::setImageHeight(int height) {
    this->imageHeight = height;
}

/**
 * @function MaskImage
 * @Note Uses the grayscale image (mask) to clear the source (image). If a pixel
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


/**
 * @function MaskImage
 * @date 4-January-2021
 * @param source image the masking operation is applied to
 * @param mask the thermal image that is used as the mask.
 * @param maskThreshold integer specifying the base value used in the mask operation
 * @return masked copy of the source image.
 */
cv::Mat ImageAnalysis::MaskImage(cv::Mat source, cv::Mat mask, float maskThreshold) {
    cv::Mat maskedImage;
    maskedImage = source.clone();

    printf("maskThreshold = %f\n", maskThreshold);

    if (maskThreshold < 20)
        maskThreshold = 20;

    for (int y = 0; y < source.rows; y++ ) {
        for (int x = 0; x < source.cols; x++ ) {
            if (mask.at<cv::Vec3b>(y,x)[0] <= maskThreshold) {
                maskedImage.at<cv::Vec3b>(y,x)[0] = 0;
                maskedImage.at<cv::Vec3b>(y,x)[1] = 0;
                maskedImage.at<cv::Vec3b>(y,x)[2] = 0;
            }
        } //x
    } // y

    return maskedImage;
} // maskedImage


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


/** @function setContrastBalance */
cv::Mat ImageAnalysis::setContrastBalance(double contrast_alpha, int brightness_beta,cv::Mat src)
{
    cv::Mat adjusted = src.clone();
    for (int y = 0; y < src.rows; y++ ) {
        for (int x = 0; x < src.cols; x++ ) {
            for (int c = 0; c < src.channels(); c++ ) {
                adjusted.at<cv::Vec3b>(y,x)[c] =
                  cv::saturate_cast<uchar>(contrast_alpha*src.at<cv::Vec3b>(y,x)[c] + brightness_beta ); //Vec3b
            }
        }
    }

    return adjusted;
} // setContrastBalance


/** @function FireLine_LineFit
 *  @date 28 December 2020
 *  @note This method fits a line to the fire line progression vector origination
 *        points.  This establishes a line through and along the flame front
 *        from which the flame front depth can be measured.
 */
void ImageAnalysis::FireLine_LineFit(std::vector<point3F> p, cv::Point *A, cv::Point *B) {
    size_t i, j;
    float sumX=0, sumX2=0, sumY=0, sumXY=0;
    float meanX=0, meanY=0, xSD, ySD, a=0;
    float filteredCount=0;
    bool done;


    float n = p.size();
    for (j=0; j < p.size(); j++) {
        meanX += p[j].x;
        meanY += p[j].y;
    }
    meanX = meanX / n;
    meanY = meanY / n;

    for (j=0; j < p.size(); j++) {
        float xdif = p[j].x - meanX;
        float ydif = p[j].y - meanY;

        sumX += xdif * xdif;
        sumY += ydif * ydif;
    }
    xSD = sqrt(sumX / (n-1));
    ySD = sqrt(sumY / (n-1));

    //printf("X std dev: %f  Y std dev: %f\n", xSD, ySD);

    sumY = 0; sumX = 0;
    float minX = meanX - 2.0 * xSD;
    float maxX = meanX + 2.0 * xSD;
    float minY = meanY - 2.0 * ySD;
    float maxY = meanY + 2.0 * ySD;

    for (j=0; j < p.size(); j++) {
        if (p[j].x >= minX && p[j].x <= maxX &&
            p[j].y >= minY && p[j].y <= maxY) {
            sumX += p[j].x;
            sumY += p[j].y;
            filteredCount++;
        }
    } // j

    meanX = sumX / filteredCount;
    meanY = sumY / filteredCount;

    float meanDif=0, xMeanDif2=0, m, b;
    for (j=0; j < p.size(); j++) {
        if (p[j].x >= minX && p[j].x <= maxX &&
            p[j].y >= minY && p[j].y <= maxY) {

            meanDif += (p[j].x-meanX) * (p[j].y-meanY);
            xMeanDif2 += pow((p[j].x-meanX),2);
        }
    } // j

    m = meanDif / xMeanDif2;
    b = meanY - m * meanX;

    // y = mx + b

    /*
     * Pick two reasonable, but likely less than optimal, starting points on the
     * fitted line for the flame front depth analysis
     */
    done = false;
    j = 2;
    while (!done) {
        if (p[j].x >= minX && p[j].x <= maxX &&
            p[j].y >= minY && p[j].y <= maxY) {
            A->x = p[j].x;
            A->y = m * A->x + b;
            done = true;
        }

        j++;
        if (j > p.size()) {
            A->x = 0;
            A->y = 0;
            //puts("Done without match, incr loop");
            done = true;
        }
    } // done

    done = false;
    j = p.size() - 2;
    while (!done) {
        if (p[j].x >= minX && p[j].x <= maxX &&
            p[j].y >= minY && p[j].y <= maxY) {
            B->x = p[j].x;
            B->y = m * B->x + b;
            done = true;
        }

        j--;
        if (j < 0) {
            B->x = 0;
            B->y = 0;
            //puts("Done without match, decr loop");
            done = true;
        }
    } // done

    //printf("FireLine fit Initial A:  %d, %d\n", A->x, A->y);
    //printf("FireLine fit Initial B:  %d, %d\n", B->x, B->y);

    /*
     * Using the sub-optimal starting points, Optimize the starting points such
     * that they are the further apart than any other point pairs.
     */

    float bestDistance = sqrt(pow(A->x - B->x,2.0)+pow(A->y-B->y,2.0));
    float thisDistance;
    point3F C, D;
    for (i=0; i < p.size()-1; i++) {
        C.x = p[i].x;
        C.y = m * C.x + b;
        for (j=i+1; j<p.size(); j++) {
            D.x = p[j].x;
            D.y = m * D.x + b;
            if (D.y >= 0 && D.y <=this->imageHeight &&
                C.y >=0 && C.y <= this->imageHeight) {
                thisDistance = sqrt(pow(C.x - D.x, 2.0) + pow(C.y - D.y, 2.0));
                if (thisDistance > bestDistance) {
                    bestDistance = thisDistance;
                    A->x = C.x;
                    A->y = C.y;
                    B->x = D.x;
                    B->y = D.y;
                }
            }
        } // j
    } // i

    /*
    printf("n= %f  filteredCount=%f\n", n, filteredCount);
    printf("FireLine Line fit results m=%f  b=%f\n", m, b);
    printf("FireLine fit A:  %d, %d\n", A->x, A->y);
    printf("FireLine fit B:  %d, %d\n", B->x, B->y);
    */
}// FireLine_LineFit


/** @function AnalyzeProgressionVector
 *  @date 9 December 2020
 *  @note Processes the fire line progression normals(vectors) and determines a mean
 *        progression direction and speed. Removes observations > 2 x std deviation.
 *        May want to explore a more rigorous outlier removal method, ie X^2, if this
 *        overall processing approach is valid.
 */
void ImageAnalysis::AnalyzeProgressionVector(std::vector<point3F> p,
                                             std::vector<Vector3F> d,
                                             std::vector<float> magnitude)
{
    int j;
    float filteredCount;
    float difX, difY, difDX, difDY, difMagnitude;
    float sumX=0, sumY=0, sumDX=0, sumDY=0, sumMagnitude=0;
    float sdX, sdY, sdDX, sdDY, sdMagnitude;
    float minX, minY, minDX, minDY, minMagnitude;
    float maxX, maxY, maxDX, maxDY, maxMagnitude;
    int n = p.size();

    this->flameFrontCenter.x = 0.0;
    this->flameFrontCenter.y = 0.0;
    this->flameFrontDirection.x = 0.0;
    this->flameFrontDirection.y = 0.0;
    this->flameFrontSpreadMagnitude = 0.0;

    for (j=0; j < n; j++) {
        this->flameFrontCenter.x += p[j].x;
        this->flameFrontCenter.y += p[j].y;
        this->flameFrontDirection.x += d[j].x;
        this->flameFrontDirection.y += d[j].y;
        this->flameFrontSpreadMagnitude += magnitude[j];
    }

    this->flameFrontCenter.x = this->flameFrontCenter.x / n;
    this->flameFrontCenter.y = this->flameFrontCenter.y / n;
    this->flameFrontDirection.x = this->flameFrontDirection.x / n;
    this->flameFrontDirection.y = this->flameFrontDirection.y / n;
    this->flameFrontSpreadMagnitude = this->flameFrontSpreadMagnitude / n;

    printf("Original Vectors:");
    printf("  Direction Summary: (%f,%f) --> (%f,%f)   <%f>\n",
           this->flameFrontCenter.x, this->flameFrontCenter.y,
           this->flameFrontDirection.x, this->flameFrontDirection.y,
           this->flameFrontSpreadMagnitude);

    for (j=0; j < n; j++) {
        difX = p[j].x - this->flameFrontCenter.x;
        difY = p[j].y - this->flameFrontCenter.y;
        difDX = d[j].x - this->flameFrontDirection.x;
        difDY = d[j].y - this->flameFrontDirection.y;
        difMagnitude = magnitude[j] - this->flameFrontSpreadMagnitude;

        sumX += difX * difX;
        sumY += difY * difY;
        sumDX += difDX * difDX;
        sumDY += difDY * difDY;
        sumMagnitude += difMagnitude * difMagnitude;
    } //j

    sdX = sqrt(sumX / (n-1));
    sdY = sqrt(sumY / (n-1));
    sdDX = sqrt(sumDX / (n-1));
    sdDY = sqrt(sumDY / (n-1));
    sdMagnitude = sqrt(sumMagnitude / (n-1));

    minX = this->flameFrontCenter.x - 2 * sdX; maxX = this->flameFrontCenter.x + 2 * sdX;
    minY = this->flameFrontCenter.y - 2 * sdY; maxY = this->flameFrontCenter.y + 2 * sdY;
    minDX = this->flameFrontDirection.x - 2 * sdDX; maxDX = this->flameFrontDirection.x + 2 * sdDX;
    minDY = this->flameFrontDirection.y - 2 * sdDY; maxDY = this->flameFrontDirection.y + 2 * sdDY;
    minMagnitude = this->flameFrontSpreadMagnitude - 2 * sdMagnitude;
    maxMagnitude = this->flameFrontSpreadMagnitude + 2 * sdMagnitude;

    this->flameFrontCenter.x = 0; this->flameFrontCenter.y = 0;
    this->flameFrontDirection.x = 0; this->flameFrontDirection.y = 0;
    this->flameFrontSpreadMagnitude = 0;
    filteredCount = 0;
    for (j=0; j < n; j++) {
        if (p[j].x >= minX && p[j].x <= maxX &&
            p[j].y >= minY && p[j].y <= maxY &&
            d[j].x >= minDX && d[j].x <= maxDX &&
            d[j].y >= minDY && d[j].y <= maxDY &&
            magnitude[j] >= minMagnitude &&
            magnitude[j] <= maxMagnitude)
        {
            this->flameFrontCenter.x += p[j].x;
            this->flameFrontCenter.y += p[j].y;
            this->flameFrontDirection.x += d[j].x;
            this->flameFrontDirection.y += d[j].y;
            this->flameFrontSpreadMagnitude += magnitude[j];
            filteredCount++;
        }
    } // j

    this->flameFrontCenter.x = this->flameFrontCenter.x / filteredCount;
    this->flameFrontCenter.y = this->flameFrontCenter.y / filteredCount;
    this->flameFrontDirection.x = this->flameFrontDirection.x / filteredCount;
    this->flameFrontDirection.y = this->flameFrontDirection.y / filteredCount;
    this->flameFrontSpreadMagnitude = this->flameFrontSpreadMagnitude / filteredCount;
}// AnalyzeProgressionVector


/**
 * @function getFlameFrontDirection
 * @date 10-December-20
 * @return
 */
Vector2F ImageAnalysis::getFlameFrontDirection() {
    return this->flameFrontDirection;
}


/**
 * @function getFlameFrontCenter
 * @date 10-December-20
 * @return
 */
 point2F ImageAnalysis::getFlameFrontCenter() {
    return this->flameFrontCenter;
}


/**
 * @function getFlameFrontSpreadMagnitude
 * @date 10-December-20
 * @return
 */
float ImageAnalysis::getFlameFrontSpreadMagnitude() {
    return this->flameFrontSpreadMagnitude;
}

/**
 * @function printProgressionSummary
 * @date 29-December-2020
 * @note Prints the fire line progression information in a tabular format.  This method
 *       displays the progression data in translated metric units and compass direction.
 */
void ImageAnalysis::printProgressionSummary(std::vector<point2F> fireCenter,
                                            std::vector<Vector2F> fireSpreadDirection,
                                            std::vector<float> fireSpreadMagnitude,
                                            std::vector<float> flameFrontDepth,
                                            std::vector<float> meanThermal,
                                            std::vector<float> medianThermal,
                                            std::vector<int> minThermal,
                                            std::vector<int> maxThermal,
                                            std::vector<float> radiometerFlameWidth,
                                            std::vector<float> radiometerFlameDepth,
                                            std::vector<float> radiometerA,
                                            std::vector<float> radiometerB,
                                            std::vector<float> fractionalFireArea,
                                            MCA_Config config, int imageWidth, int imageHeight)
{
    int i;
    float exposureTime=0;
    float FOV = config.getFOV();
    float pixelsPerMeter= FOV / (float) imageWidth;
    float compassHeading;

    // print the header
    printf(" Image Width (meters): %5.1f\n", (float) imageWidth * pixelsPerMeter);
    printf("Image Height (meters): %5.1f\n", (float) imageHeight * pixelsPerMeter);
    printf("     Meters Per pixel: %7.4f\n", pixelsPerMeter);
    printf("                                                                      Flame Front   ----------- Thermal -----------    ------------ Radiometer --------------\n");
    printf("           Time        Fireline Center*     Progression     Speed        Depth      -------------------------------     Flame  Flame               Fractional\n");
    printf(" Frame   (seconds)      x pos.   y pos       Direction     meter/sec    (meters)     Mean     Median     Min    Max     Width  Depth     A     B      Area\n");

    for (i=0; i<fireCenter.size(); i++) {
/*
        printf("radiometerFlameWidth %5.2f \n", radiometerFlameWidth[i]);
        printf("radiometerFlameDepth %5.2f\n", radiometerFlameDepth[i]);
        printf("radiometerA  %6.3f\n", radiometerA[i]);
        printf("radiometerB  %6.3f\n", radiometerB[i]);
        printf("fractionalFireArea  %6.3f\n", fractionalFireArea[i]);
*/
        if (fireCenter[i].x == -1) {
            // fire line not present in image / or insufficient data
            printf(" %4d     %6.2f        I-D       I-D           I_D           I_D         I_D       %5.1f   %5.0f      %4d    %4d       I_D    I_D     I_D    I_D     I_D\n",
                   i+1, exposureTime, meanThermal[i], medianThermal[i], minThermal.at(i), maxThermal.at(i));
        } else {
            compassHeading = getCompassHeading(fireSpreadDirection[i]);
            printf(" %4d     %6.2f      %5.1f     %5.1f         %5.1f         %5.2f", i+1,
                   exposureTime, fireCenter[i].x * pixelsPerMeter, fireCenter[i].y * pixelsPerMeter,
                   compassHeading, fireSpreadMagnitude[i]*pixelsPerMeter);
            if (flameFrontDepth[i] > 0)
                printf("        %5.2f    ",flameFrontDepth[i] * pixelsPerMeter);
            else
                printf("         I_D     ");
            printf("  %5.1f   %5.0f      %4d    %4d   ",
                   meanThermal[i], medianThermal[i], minThermal.at(i), maxThermal.at(i));
            printf("  %5.2f  %5.2f  %6.3f  %6.3f  %6.3f\n",
                   radiometerFlameWidth[i], radiometerFlameDepth[i],
                   radiometerA[i], radiometerB[i], fractionalFireArea[i]);
        } // Sufficient Data?
        exposureTime += config.getElapsedTime(); // elapsed time in seconds between photos
    } // i

    printf("-------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("* Fire line center is in meters from top, left corner of image.\n");
    printf(" Progression Direction: is compass degree heading of dominant progression.\n");
    printf("I-D:: Insufficient data from image processing to determine fire position and progression.\n\n");
}

/**
 * @function getCompassHeading
 * @date 14-December-2020
 * @note  determines the direction of flame front based on x,y vector direction.
 * @todo add code to adjust compass heading based on photo orientation specified in
 *       config file, ie: NORTH is left edge of photo.
 * @param fireSpreadDirection
 * @return float compass heading direction of flame front.
 */
int ImageAnalysis::getCompassHeading(Vector2F fireSpreadDirection)
{
    float x = fireSpreadDirection.x;
    float y = fireSpreadDirection.y* -1.0; // invert sign on Y to correspond to image Y. Where
                                           // image y at top of frame is 0, bottom of frame is
                                           // 480 (typically with current imagery dimensions).

    int angle = atan2(y,x) * 57.29578;
    if (angle < 0)
        angle += 360;

    angle = (450 - angle) % 360;
    return angle;
}

/**
 * @function getMeanBackground
 * @note Commented out code to loop through # of channels attached to the image. The thermal
 *       imagery is stored as a 3 channel gray scale.  Thus, we need to examine only one channel.
 *       In this case, processing is limited to the first channel (channel[0]).
 * @date 14-December-2020
 * @param image
 * @return
 */
float ImageAnalysis::getMeanBackground(cv::Mat image, thermalArea tAOI) {
    int x, y;
    float sum = 0, count = 0, mean;

    for (y = 0; y < image.rows; y++)
        for (x = 0; x < image.cols; x++) {
            //for (c = 0; c < image.channels(); c++) { // only one channel in grayscale...
            if (tAOI.isInAOI(x,y)) {
                sum += image.at<cv::Vec3b>(y, x)[0];
                count++;
            }
        }
    mean = sum / count;

    return mean;
}

/**
 * @function getMedianBackground
 * @date 14-December-2020
 * @note Commented out code to loop through # of channels attached to the image. The thermal
 *       imagery is stored as a 3 channel gray scale.  Thus, we need to examine only one channel.
 *       In this case, processing is limited to the first channel (channel[0]).
 * @param image
 * @return
 */
float ImageAnalysis::getMedianBackground(cv::Mat image, thermalArea tAOI) {
    std::vector<float> pixel;
    int x, y, midpt;

    for (y = 0; y < image.rows; y++)
        for (x = 0; x < image.cols; x++) {
            //for (c = 0; c < image.channels(); c++)  // only one channel in grayscale...
            if (tAOI.isInAOI(x,y))
                pixel.push_back(image.at<cv::Vec3b>(y, x)[0]);
        }
    std::sort(pixel.begin(), pixel.end());

    if (pixel.size() <= 1)
        return 0;

    midpt = pixel.size() / 2;

    return pixel.at(midpt);
}

int ImageAnalysis::getMinBackground(cv::Mat image, thermalArea tAOI) {
    int min = 1000;

    for (int y = 0; y < image.rows; y++)
        for (int x = 0; x < image.cols; x++) {
            //for (c = 0; c < image.channels(); c++) { // only one channel in grayscale...
            if (tAOI.isInAOI(x,y)) {
                 if (image.at<cv::Vec3b>(y, x)[0] < min)
                     min = image.at<cv::Vec3b>(y, x)[0];
             }
        }

    return min;
}

int ImageAnalysis::getMaxBackground(cv::Mat image, thermalArea tAOI) {
    int max = -1;

    for (int y = 0; y < image.rows; y++)
        for (int x = 0; x < image.cols; x++) {
            //for (c = 0; c < image.channels(); c++) { // only one channel in grayscale...
            if (tAOI.isInAOI(x,y)) {
                if (image.at<cv::Vec3b>(y, x)[0] > max)
                    max = image.at<cv::Vec3b>(y, x)[0];
            }
        }
    return max;
} //getMaxBackground


/**
 * @function optimizeFlameFrontLine
 * @date 4-January-2021
 * @note Determines opposite normal vectors for the line AB. The method then
 *       moves the line along the normals to a location that locally maximizes
 *       the actively burning area described by AB.
 * @param A starting point of a line fitted to the flame front
 * @param B ending point of a line fitted to the flame front
 * @param fireRegions vector of multiple contours that define the burning area
 */
void ImageAnalysis::optimizeFlameFrontLine(cv::Point *A, cv::Point *B,
                                           std::vector<std::vector<cv::Point>> fireRegions,
                                           int Image_X_Max, int Image_Y_Max)
{
    Vector3F normalVector1, normalVector2;
    point3F intercept1, intercept2;
    point3F P;
    cv::Point A_prime, B_prime;
    cv::Point A_optimal, B_optimal;
    float optimalArea = 0.0, thisArea;

    Normal N = Normal(*A, *B, Image_X_Max, Image_Y_Max);
    N.CalculateNormals(&normalVector1, &normalVector2);

    Vector3F AB;
    AB.x = B->x - A->x;
    AB.y = B->y - A->y;
    float magnitudeAB = sqrt(pow(AB.x,2)+pow(AB.y,2));
    // Normalize AB
    AB.x = AB.x/magnitudeAB;
    AB.y = AB.y/magnitudeAB;
    AB.z = 0.0;

    A_optimal.x = A->x; A_optimal.y = A->y;
    B_optimal.x = B->x; B_optimal.y = B->y;

    for (int j=-150; j<=150; j+=15) {
        thisArea = 0.0;
        A_prime.x = A->x + j * normalVector1.x; A_prime.y = A->y + j * normalVector1.y;
        B_prime.x = B->x + j * normalVector1.x; B_prime.y = B->y + j * normalVector1.y;

        for (int i=0; i<fireRegions.size(); i++) {
            if (fireRegions[i].size() >= FLAMEFRONT_MIN_PT_SIZE) {
                for (float d = 0; d < magnitudeAB; d++) {
                    P.x = A_prime.x + d * AB.x;  //AB is direction vector
                    P.y = A_prime.y + d * AB.y;
                    P.z = 0.0;

                    if (N.intercept(P, normalVector1, fireRegions[i], &intercept1) &&
                        N.intercept(P, normalVector2, fireRegions[i], &intercept2))
                        thisArea += sqrt(pow(intercept2.x - intercept1.x, 2.0) +
                                         pow(intercept2.y - intercept1.y, 2.0));
                } // d
            } // flame front contour area is large enough to process
        } //i

        if (thisArea > optimalArea) {
            A_optimal.x = A_prime.x;
            A_optimal.y = A_prime.y;
            B_optimal.x = B_prime.x;
            B_optimal.y = B_prime.y;
            optimalArea = thisArea;
        }

        if (j > 0) {
            thisArea = 0.0;
            A_prime.x = A->x + j * normalVector2.x; A_prime.y = A->y + j * normalVector2.y;
            B_prime.x = B->x + j * normalVector2.x; B_prime.y = B->y + j * normalVector2.y;

            for (int i=0; i<fireRegions.size(); i++) {
                if (fireRegions[i].size() >= FLAMEFRONT_MIN_PT_SIZE) {
                    for (float d = 0; d < magnitudeAB; d++) {
                        P.x = A_prime.x + d * AB.x;  //AB is direction vector
                        P.y = A_prime.y + d * AB.y;
                        P.z = 0.0;

                        if (N.intercept(P, normalVector1, fireRegions[i], &intercept1) &&
                            N.intercept(P, normalVector2, fireRegions[i], &intercept2))
                            thisArea += sqrt(pow(intercept2.x - intercept1.x, 2.0) +
                                             pow(intercept2.y - intercept1.y, 2.0));
                    } // d
                } // flame front contour area is large enough to process
            } //i

            if (thisArea > optimalArea) {
                A_optimal.x = A_prime.x;
                A_optimal.y = A_prime.y;
                B_optimal.x = B_prime.x;
                B_optimal.y = B_prime.y;
                optimalArea = thisArea;
            }
        } // examine opposing normal
    } // i

    A->x = A_optimal.x;
    A->y = A_optimal.y;
    B->x = B_optimal.x;
    B->y = B_optimal.y;
} // optimizeFlameFrontLine


/**
 * @function getFlameFrontDepth
 * @date 9-February-2021
 * @note Determines the average depth of the flame front
 * @param A starting point of a line fitted to the flame front
 * @param B ending point of a line fitted to the flame front
 * @param fireRegions vector of multiple contours that define the burning area
 * @param InterceptPt1 first intercept point of normal from AB to fireRegion
 * @param InterceptPt2 second intercept point of normal from AB to fireRegion
 * @param Image_X_Max pixel width of imagery
 * @param Image_Y_Max pixel height of imagery
 * @return float, average cross-sectional (perpendicular to the fitted line) width of
 *                the flame front.
 */
float ImageAnalysis::getFlameFrontDepth(cv::Point A, cv::Point B,
                                        std::vector<std::vector<cv::Point>> fireRegions,
                                        std::vector<cv::Point> *InterceptPt1,
                                        std::vector<cv::Point> *InterceptPt2,
                                        int Image_X_Max, int Image_Y_Max)
{
    float AvgFlameDepth=0;
    Vector3F normalVector1, normalVector2;
    point3F intercept1, intercept2;

    Normal N = Normal(A,B, Image_X_Max, Image_Y_Max);
    N.CalculateNormals(&normalVector1, &normalVector2);

    Vector3F AB;
    AB.x = B.x - A.x;
    AB.y = B.y - A.y;
    float magnitudeAB = sqrt(pow(AB.x,2)+pow(AB.y,2));
    // Normalize AB
    AB.x = AB.x/magnitudeAB;
    AB.y = AB.y/magnitudeAB;
    AB.z = 0.0;
/*
    printf("getFlameFrontDepth:: A: %d,%d  B: %d,%d  Magnitude:%f\n",
           A.x, A.y, B.x, B.y, magnitudeAB);
    printf("getFlameFrontDepth:: NormalVector1: P --> (%f,%f,%f)\n",
            normalVector1.x, normalVector1.y, normalVector1.z);
    printf("getFlameFrontDepth:: NormalVector2: P --> (%f,%f,%f)\n",
            normalVector2.x, normalVector2.y, normalVector2.z);
    printf("getFlameFrontDepth::  Number of Regions %ld\n",fireRegions.size());
*/

    for (int i=0; i<fireRegions.size(); i++) {
        if (fireRegions[i].size() >= FLAMEFRONT_MIN_PT_SIZE) {
            /*
             * step along the AB line in approx 1 pixel increments.  For each step
             * determine the intersection point of the normals 1 & 2 with the fire line
             * contour.  the sum of the magnitudes of the intersecting vectors is the
             * flame front depth at that point.
             */
            point3F P;
            for (float d = 0; d < magnitudeAB; d++) {
                P.x = A.x + d * AB.x;  //AB is direction vector
                P.y = A.y + d * AB.y;
                P.z = 0.0;

                if (N.intercept(P, normalVector1, fireRegions[i], &intercept1) &&
                    N.intercept(P, normalVector2, fireRegions[i], &intercept2)) {

                    // distance between intercepts 1 & 2 is flame front depth.
                    float flameDepth = sqrt(pow(intercept2.x - intercept1.x, 2.0) +
                                            pow(intercept2.y - intercept1.y, 2.0));
                    AvgFlameDepth += flameDepth;

                    //printf("getFlameFrontDepth:: (%f) d:%f Point:%f,%f  Intercept1:%f,%f  Intercept2:%f,%f  Depth:%f\n",
                    //      interceptCount, d, P.x, P.y, intercept1.x, intercept1.y, intercept2.x, intercept2.y, flameDepth);

                    InterceptPt1->push_back(cv::Point(intercept1.x, intercept1.y));
                    InterceptPt2->push_back(cv::Point(intercept2.x, intercept2.y));
                } // Both Intercepts found
            } // d
        }
    } //i

    float interceptCount = InterceptPt1->size();

    if (interceptCount > 0)
        AvgFlameDepth = AvgFlameDepth/interceptCount;
    else
        AvgFlameDepth = -1;

    return AvgFlameDepth;
}


/**
 * @function printProgressionSummary
 * @date 8-February-2021
 * @note Prints the fire line progression information in a comma-separated (CSV) format
 *       to the file path specified in the reportPath parameter.  This method
 *       displays the progression data in translated metric units and compass direction.
 */
int ImageAnalysis::printProgressionSummary(char reportPath[PATH_SIZE], std::vector<point2F> fireCenter,
                                           std::vector<Vector2F> fireSpreadDirection,
                                           std::vector<float> fireSpreadMagnitude, std::vector<float> flameFrontDepth,
                                           std::vector<float> meanThermal, std::vector<float> medianThermal,
                                           std::vector<int> minThermal,
                                           std::vector<int> maxThermal,
                                           std::vector<float> radiometerFlameWidth,
                                           std::vector<float> radiometerFlameDepth,
                                           std::vector<float> radiometerA,
                                           std::vector<float> radiometerB,
                                           std::vector<float> fractionalFireArea,
                                           MCA_Config config,
                                           int imageWidth, int imageHeight)
{
    int i;
    float exposureTime=0;
    float FOV = config.getFOV();
    float pixelsPerMeter= FOV / imageWidth;
    float compassHeading;
    FILE *reportStream;

    if (!(reportStream = fopen(reportPath,"w"))) {
        printf("Error Unable to open report file %s\n", reportPath);
        return 0;
    }

    // print the header
    fprintf(reportStream,"\"Image Width (meters)\", %5.1f\n", imageWidth * pixelsPerMeter);
    fprintf(reportStream,"\"Image Height (meters)\", %5.1f\n", imageHeight * pixelsPerMeter);
    fprintf(reportStream,"\"Meters Per pixel\", %7.4f\n", pixelsPerMeter);
    fprintf(reportStream,"\" \", \" \", \" \", \" \", \" \", \" \", \" \", \"Flame\", \" \", \" \", \" \", \" \", \"Radiometer\",\"Radiometer\",\"Radiometer\",\"Radiometer\",\"Radiometer\"\n");
    fprintf(reportStream,"\" \", \" \", \" \", \"Center(1)\", \"Center(1)\", \"Progression\", \"Progression\", \"Front\", \"Thermal\", \"Thermal\", \"Thermal\", \"Thermal\",\"Flame\",\"Flame\",\"AOI\",\"AOI\",\"Fractional\"\n");
    fprintf(reportStream," \"Frame\", \"Time\", \"x pos\", \"y pos \", \"Direction(2)\", \"Speed(7)\", \"Depth(5)\", \"Mean\", \"Median\", \"Min\", \"Max\", \"Width(6)\", \"Depth\",\"A(3)\",\"B(4)\",\"Area\"\n");
    fprintf(reportStream," \"\", \"(seconds)\", \"(meters)\", \"(meters)\",   \"(degrees)\", \"(meter/sec)\", \"(meters)\",\"-\", \"-\",\"-\",\"-\",\"(meters)\", \"(meters)\",\"(sq. meters)\",\"(sq. meters)\",\"-\"\n");
    for (i=0; i<fireCenter.size(); i++) {
        compassHeading = getCompassHeading(fireSpreadDirection[i]);
        fprintf(reportStream,"%4d,%6.2f,%5.1f,%5.1f,%5.1f,%6.3f,", i+1,
            exposureTime, fireCenter[i].x * pixelsPerMeter, fireCenter[i].y * pixelsPerMeter,
            compassHeading, fireSpreadMagnitude[i]*pixelsPerMeter);

        fprintf(reportStream,"%5.2f,%5.1f,%5.0f,%4d,%4d,",
                flameFrontDepth[i] * pixelsPerMeter, meanThermal[i],
                medianThermal[i], minThermal.at(i), maxThermal.at(i));

        fprintf(reportStream,"%5.2f,%5.2f,%6.3f,%6.3f,%6.3f\n",
               radiometerFlameWidth[i], radiometerFlameDepth[i],
               radiometerA[i], radiometerB[i], fractionalFireArea[i]);

        exposureTime += config.getElapsedTime(); // elapsed time in seconds between photos
    } // i

    fprintf(reportStream,"\"(1) Fire line center is in meters from top, left corner of image.\"\n");
    fprintf(reportStream,"\"(2) Progression Direction: is compass degree heading of dominant\"\n");
    fprintf(reportStream,"\"    progression direction.\"\n");
    fprintf(reportStream,"\"(3) Active burning area inside Radiometer AOI in square meters.\"\n");
    fprintf(reportStream,"\"(4) Non-burning area inside Radiometer AOI in square meters.\"\n");
    fprintf(reportStream,"\"(5) Average depth in meters of flame front across the active burning area.\"\n");
    fprintf(reportStream,"\"    Determined by calculating normals to the fitted progression line,\"\n");
    fprintf(reportStream,"\"    then finding the intersection points of the normals with the burning\"\n");
    fprintf(reportStream,"\"    area contour. Depth is average of summed opposing normal vector\"\n");
    fprintf(reportStream,"\"    magnitudes\"\n");
    fprintf(reportStream,"\"(6) Flame front width in meters. Determined by fitted line to the\"\n");
    fprintf(reportStream,"\"    flame front progression normals and the lines intersection\"\n");
    fprintf(reportStream,"\"    points (if any) with the radiometer AOI boundary.\"\n");
    fprintf(reportStream,"\"(7) Progression speed is the mean length of the normal vectors.\"\n");
    fprintf(reportStream,"\"    across the flame front.  These vectors describe the distance the\"\n");
    fprintf(reportStream,"\"    flame spread between frames at the specified time interval between images\"\n");
    fprintf(reportStream,"\"Important Note: Negative value means insufficient data from image\"\n");
    fprintf(reportStream,"\"    processing to make a determination.\"\n");

    return 1;
}

/**
 * @function printProgressionSummary
 * @date 4-January-2021
 * @note This method "blanks out" pixels of the source image that are not inside the
 *       radiometer pixel AOI.
 * @param radiometerMask defines boundaries/areas of the radiometer pixel
 * @param image to which the mask operation is applied.
 * @returns cv:Mat masked image
 */
cv::Mat ImageAnalysis::logicalAndImage(radiometer radiometerMask, cv::Mat image) {
    cv::Mat result = image.clone();

    for (int y = 0; y < result.rows; y++ ) {
        for (int x = 0; x < result.cols; x++) {
            if (result.channels() == 1) {
                if (!radiometerMask.isInAOI(x, y))
                    result.at<uchar>(y, x) = 0;
            } else {
                for (int c = 0; c < result.channels(); c++) {
                    if (!radiometerMask.isInAOI(x, y))
                        result.at<cv::Vec3b>(y, x) = 0;
                } // channels
            } // single channel image?
        } //x
    } // y

    return result;
}// logicalAndImage


/**
 * @function ImageAnalysis
 * @date 10-February-2021
 * @note This method takes the optimized flame front line and clips it to the flame front
 *       contour contained withing the radiometer AOI.
 * @param A
 * @param B
 * @param radContours
 */
void ImageAnalysis::clipFrameFrontLine(cv::Point *A, cv::Point *B,
                                       std::vector<std::vector<cv::Point>> radContours)
{
    Vector3F normal_PA, normal_PB;
    point3F thisIntercept, intercept1, intercept2;
    float thisMagnitude;
    float bestMagnitudePA = -10000000.0;
    float bestMagnitudePB = -10000000.0;
    Normal N = Normal(*A, *B, this->imageWidth, this->imageHeight);

    //printf("IA:clipFrameFrontLine: initial A: %d,%d\n", A->x, A->y);
    //printf("IA:clipFrameFrontLine: initial B: %d,%d\n", B->x, B->y);

    intercept1.x = A->x;    intercept1.y = A->y;
    intercept2.x = B->x;    intercept2.y = B->y;

    point3F P;
    P.x = (A->x + B->x)/2.0;
    P.y = (A->y + B->y)/2.0;
    P.z = 0.0;

    //printf("IA:clipFrameFrontLine: MidPoint %f,%f\n", P.x, P.y);

    // First intercept is from the mid point P toward point A.
    float PA_magnitude = sqrt(pow((A->x-P.x),2.0)+pow((A->y-P.y),2.0));
    normal_PA.x = (A->x - P.x)/PA_magnitude;
    normal_PA.y = (A->y - P.y)/PA_magnitude;
    normal_PA.z = 0.0;

    // second intercept is from the mid point P toward point B.
    float PB_magnitude = sqrt(pow((B->x-P.x),2.0)+pow((B->y-P.y),2.0));
    normal_PB.x = (B->x - P.x)/PB_magnitude;
    normal_PB.y = (B->y - P.y)/PB_magnitude;
    normal_PB.z = 0.0;

    for (int i=0; i<radContours.size(); i++) {
        if (radContours[i].size() >= FLAMEFRONT_MIN_PT_SIZE) {
            if (N.interceptRadiometerAOI(P, normal_PA, radContours[i], &thisIntercept)) {

                    thisMagnitude = sqrt(pow((thisIntercept.x - P.x), 2.0) +
                                         pow((thisIntercept.y - P.y), 2.0));
                    if (thisMagnitude > bestMagnitudePA) {
                        bestMagnitudePA = thisMagnitude;
                        intercept1.x = thisIntercept.x;
                        intercept1.y = thisIntercept.y;
                        intercept1.z = thisIntercept.z;
                    } // more distant solution point?

            }

            if (N.interceptRadiometerAOI(P, normal_PB, radContours[i], &thisIntercept)) {

                    thisMagnitude = sqrt(pow((thisIntercept.x - P.x), 2.0) +
                                         pow((thisIntercept.y - P.y), 2.0));
                    if (thisMagnitude > bestMagnitudePB) {
                        bestMagnitudePB = thisMagnitude;
                        intercept2.x = thisIntercept.x;
                        intercept2.y = thisIntercept.y;
                        intercept2.z = thisIntercept.z;
                    } // more distant solution point?

            }
        } // sufficient size contour for analysis
    } // i

    /**
     * At this point we should have the two most distant points in the contour set
     * contained within the radiometer AOI.  So need to update the initial A & B
     * intercept points.
     */
    A->x = intercept1.x;    A->y = intercept1.y;
    B->x = intercept2.x;    B->y = intercept2.y;

    //printf("IA:clipFrameFrontLine: clipped A: %d,%d\n", A->x, A->y);
    //printf("IA:clipFrameFrontLine: clipped B: %d,%d\n", B->x, B->y);
}// clipFrameFrontLine

line2F ImageAnalysis::getFlameFrontEdgeDepth(std::vector<cv::Point> *InterceptPt1,
                                             std::vector<cv::Point> *InterceptPt2,
                                             int startEdge)
{
    line2F l;
    int start, stop, incr;

    if (startEdge == FROM_UPPER_EXTREME_EDGE) {
        start = InterceptPt1->size() - 1;
        incr = -1;
        stop = 0;
    } else {
        start = 0;
        incr = 1;
        stop = InterceptPt1->size() - 1;
    }

    for (int i = start; i <= stop; i+= incr) {

    } // i

    return l;
}

