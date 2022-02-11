/*
 * Note: based on source from 
 *     https://stackoverflow.com/questions/10315551/opencv-2-3-c-how-to-isolate-object-inside-image/10317919#10317919
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
  printf("RGB Imagery: %s\n", argv[1]);
  printf("Thermal Imagery: %s\n", argv[2]);
  
  cv::Mat rgbImage = cv::imread(argv[1]);
  cv::Mat thermalImage = cv::imread(argv[2]);

  int rgbImageWidth = rgbImage.cols;
  int rgbImageHeight = rgbImage.rows;
  int thermalImageWidth = thermalImage.cols;
  int thermalImageHeight = thermalImage.rows;
    
  cv::Rect rgbROI;
  rgbROI.x = SIDE_TRIM;
  rgbROI.y = 0;
  rgbROI.width = rgbImageWidth - SIDE_TRIM;
  rgbROI.height = rgbImageHeight;
  cv::Mat rgbCropped = rgbImage(rgbROI);
  cv::resize(rgbCropped, rgbCropped, cv::Size(720, 480), cv::INTER_CUBIC);


  //if (rgbImage.cols != thermalImage.cols ||
  //    rgbImage.rows != thermalImage.rows) 

 
  cv::Rect thermalROI;
  thermalROI.x = 0;
  thermalROI.y = TOP_TRIM;
  thermalROI.width = thermalImageWidth - SIDE_TRIM;
  thermalROI.height = thermalImageHeight - (TOP_TRIM + BOTTOM_TRIM);

  //printf("thermalROI:: %d,%d (%d,%d)\n", thermalROI.x, thermalROI.y, thermalROI.width, thermalROI.height);
  
  cv::Mat thermalCropped = thermalImage(thermalROI);
  cv::resize(thermalCropped, thermalCropped, cv::Size(720, 480), cv::INTER_CUBIC);
  
  cv::imshow("thermalCropped", thermalCropped);
  cv::imshow("rgbCropped", rgbCropped);
  cv::waitKey(0);
  
  std::vector<int> compression_params;
  compression_params.push_back(cv::IMWRITE_JPEG_QUALITY); // compression technique
  compression_params.push_back(98); // quality level
  
  cv::imwrite("rgbCropped.jpg", rgbCropped, compression_params);
  cv::imwrite("thermalCropped.jpg", thermalCropped, compression_params);

  return 0;
}
