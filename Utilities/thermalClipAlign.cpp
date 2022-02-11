/*
 * Note: based on source from 
 *     https://stackoverflow.com/questions/10315551/opencv-2-3-c-how-to-isolate-object-inside-image/10317919#10317919
 */

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <stdio.h>

#include <vector>

#define RGB_IMG_HEIGHT 480

// Define the offset of the RGB imagery relative to the thermal imagery.
#define TOP_TRIM    85 //56
#define SIDE_TRIM   31 //20

int main(int argc, char* argv[])
{
  char destination[256];
  strcpy(destination,"Processed/");
  strcat(destination,argv[1]);

  printf("Thermal Imagery: %s\n", argv[1]);
  printf("Destination: %s\n",destination);
 
  cv::Mat thermalImage = cv::imread(argv[1]);
  int thermalImageWidth = thermalImage.cols;
  int thermalImageHeight = thermalImage.rows;


  if (thermalImageWidth == 80) {
    cv::resize(thermalImage, thermalImage, cv::Size(720, 540), cv::INTER_LANCZOS4);
    thermalImageWidth = thermalImage.cols;
    thermalImageHeight = thermalImage.rows;
  }

  cv::Rect thermalROI;
  thermalROI.y = TOP_TRIM;

  if (SIDE_TRIM < 0) {
    thermalROI.x = 0;
     thermalROI.width = thermalImageWidth + SIDE_TRIM;
  } else {
    thermalROI.x = SIDE_TRIM;
    thermalROI.width = thermalImageWidth - SIDE_TRIM;
  }

 
  if (thermalImageHeight - TOP_TRIM < RGB_IMG_HEIGHT)
    thermalROI.height = thermalImageHeight - TOP_TRIM;
  else
    thermalROI.height = RGB_IMG_HEIGHT;
 
      //thermalROI.height = thermalImageHeight - (TOP_TRIM + BOTTOM_TRIM);
  
  cv::Mat thermalCropped = thermalImage(thermalROI);
  cv::resize(thermalCropped, thermalCropped, cv::Size(720, 480), cv::INTER_CUBIC);
  
  std::vector<int> compression_params;
  compression_params.push_back(cv::IMWRITE_JPEG_QUALITY); // compression technique
  compression_params.push_back(98); // quality level
  
  cv::imwrite(destination, thermalCropped, compression_params);

  return 0;
}
