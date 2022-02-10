#ifndef IMAGEANALYSIS_H
#define IMAGEANALYSIS_H

#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
//#include "opencv2/xphoto.hpp"
#include "opencv2/highgui.hpp"
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>

#include "normal.h"
#include "mca_config.h"
#include "fireline.h"
#include "MultiContour.hpp"
#include "thermalArea.h"
#include "radiometer.h"

#define GRAY_LEVELS 256
#define FLAMEFRONT_MIN_PT_SIZE 75 // was 100
#define FROM_UPPER_EXTREME_EDGE -1
#define FROM_LOWER_EXTREME_EDGE 1

typedef struct line2F_Struct {
    float x1, y1; // point A
    float x2, y2; // point B
} line2F;


class ImageAnalysis
{
private:
    int imageWidth;
    int imageHeight;
    point2F flameFrontCenter;
    Vector2F flameFrontDirection;
    float flameFrontSpreadMagnitude;

public:
    ImageAnalysis();
    ImageAnalysis(int width, int height);
    void setImageWidth(int width);
    void setImageHeight(int height);
    cv::Mat CalculateNDVI(cv::Mat image);
    cv::Mat MaskImage(cv::Mat image, int maskThreshold);
    cv::Mat MaskImage(cv::Mat source, cv::Mat mask, float maskThreshold);
    void imageHistogram(cv::Mat image, long int histogram[GRAY_LEVELS],
                        long int *minPixel, long int *maxPixel, bool print);
    void EdgeDetection(cv::Mat image);
    cv::Mat Erosion(int erosion_elem, int erosion_size, cv::Mat src);
    cv::Mat Dilation(int dilation_elem, int dilation_size, cv::Mat src);
    cv::Mat setContrastBalance(double contrast_alpha, int brightness_beta, cv::Mat src);
    float getMeanBackground(cv::Mat image, thermalArea tAOI);
    float getMedianBackground(cv::Mat image, thermalArea tAOI);
    int getMinBackground(cv::Mat image, thermalArea tAOI);
    int getMaxBackground(cv::Mat image, thermalArea tAOI);
    cv::Mat logicalAndImage(radiometer mask, cv::Mat image);
    point2F getFlameFrontCenter();
    Vector2F getFlameFrontDirection();
    float getFlameFrontSpreadMagnitude();
    line2F getFlameFrontEdgeDepth(std::vector<cv::Point> *InterceptPt1,
                                  std::vector<cv::Point> *InterceptPt2,
                                  int startEdge);
    float getFlameFrontDepth(cv::Point A, cv::Point B,
                             std::vector<std::vector<cv::Point>> fireRegions,
                             std::vector<cv::Point> *InterceptPt1,
                             std::vector<cv::Point> *InterceptPt2,
                             int Image_X_Max, int Image_Y_Max);
    void FireLine_LineFit(std::vector<point3F> p, cv::Point *A, cv::Point *B);
    void optimizeFlameFrontLine(cv::Point *A, cv::Point *B,
                                std::vector<std::vector<cv::Point>> fireRegions,
                                int Image_X_Max, int Image_Y_Max);
    void AnalyzeProgressionVector(std::vector<point3F> p, std::vector<Vector3F> d,
                                  std::vector<float> magnitude);
    void printProgressionSummary(std::vector<point2F> fireCenter,
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
                                 MCA_Config config,
                                 int imageWidth, int imageHeight);
    int printProgressionSummary(char reportPath[PATH_SIZE],
                                std::vector<point2F> fireCenter,
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
                                 MCA_Config config,
                                 int imageWidth, int imageHeight);
    void clipFrameFrontLine(cv::Point *A, cv::Point *B,
                            std::vector<std::vector<cv::Point>> radContours);
    int getCompassHeading(Vector2F fireSpreadDirection);
};

#endif // IMAGEANALYSIS_H
