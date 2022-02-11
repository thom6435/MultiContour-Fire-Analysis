#ifndef IMAGEANALYSIS_H
#define IMAGEANALYSIS_H

#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/xphoto.hpp"
#include "opencv2/highgui.hpp"
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>

#define GRAY_LEVELS 256

class ImageAnalysis
{

public:
    ImageAnalysis();
    cv::Mat CalculateNDVI(cv::Mat image);
    cv::Mat MaskImage(cv::Mat image, int maskThreshold);
    void imageHistogram(cv::Mat image, long int histogram[GRAY_LEVELS],
                        long int *minPixel, long int *maxPixel, bool print);
    void EdgeDetection(cv::Mat image);
    cv::Mat Erosion(int erosion_elem, int erosion_size, cv::Mat src);
    cv::Mat Dilation(int dilation_elem, int dilation_size, cv::Mat src);
};

#endif // IMAGEANALYSIS_H
