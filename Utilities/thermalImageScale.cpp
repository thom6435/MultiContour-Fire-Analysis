/*
 *  program: thermalImageScale
 * modified: 21-January-2021
 */

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <stdio.h>

#include <vector>

#define TOP_TRIM    56
#define BOTTOM_TRIM  4
#define SIDE_TRIM   20

int main(int argc, char* argv[])
{
  char destination[256];
  strcpy(destination,"Scaled/");
  strcat(destination,argv[1]);

  printf("Thermal Imagery: %s\n", argv[1]);
  printf("Destination: %s\n",destination);
 
  cv::Mat thermalImage = cv::imread(argv[1], cv::IMREAD_UNCHANGED); // reads 2 byte, 1 channel png.
  int thermalImageWidth = thermalImage.cols;
  int thermalImageHeight = thermalImage.rows;
  
  printf("Dimensions: %d x %d\n", thermalImageWidth, thermalImageHeight);
  printf("Depth: %d\n", thermalImage.depth());
  printf("Channels: %d\n", thermalImage.channels());
  
  if (thermalImageWidth == 80) {
    cv::resize(thermalImage, thermalImage, cv::Size(720, 540), cv::INTER_LANCZOS4);
    thermalImageWidth = thermalImage.cols;
    thermalImageHeight = thermalImage.rows;
  }
  
  cv::Rect thermalROI;
  thermalROI.x = 0;
  thermalROI.y = TOP_TRIM;
  thermalROI.width = thermalImageWidth - SIDE_TRIM;
  thermalROI.height = thermalImageHeight - (TOP_TRIM + BOTTOM_TRIM);
  
  cv::Mat thermalCropped = thermalImage(thermalROI);
  cv::resize(thermalCropped, thermalCropped, cv::Size(720, 480), cv::INTER_CUBIC);
 
  cv::imwrite(destination, thermalCropped);
  
  return 0;
}
