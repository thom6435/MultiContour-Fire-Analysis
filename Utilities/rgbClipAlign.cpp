/*
 * Note: based on source from 
 *     https://stackoverflow.com/questions/10315551/opencv-2-3-c-how-to-isolate-object-inside-image/10317919#10317919
 */

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <stdio.h>
#include <string.h>
#include <vector>

// Define the offset of the RGB imagery relative to the thermal imagery.
#define TOP_TRIM    85 //56
#define SIDE_TRIM   31 //20

int main(int argc, char* argv[])
{
  char destination[256];
  strcpy(destination,"Processed/");
  strcat(destination,argv[1]);
  printf("RGB Imagery: %s\n", argv[1]);
  printf("Destination: %s\n",destination);
 
  
  cv::Mat rgbImage = cv::imread(argv[1]);

  int rgbImageWidth = rgbImage.cols;
  int rgbImageHeight = rgbImage.rows;
    
  cv::Rect rgbROI;
  if (SIDE_TRIM > 0) {
    rgbROI.width = rgbImageWidth - SIDE_TRIM;
    rgbROI.x = 0;
  } else {
    rgbROI.width = rgbImageWidth + SIDE_TRIM;
    rgbROI.x = SIDE_TRIM * -1;
  }
 
  rgbROI.y = 0;
  rgbROI.height = rgbImageHeight;


  printf("RGB ROI:  %d, %d  W=%d H=%d\n", rgbROI.x, rgbROI.y,  rgbROI.width,  rgbROI.height);
  
  cv::Mat rgbCropped = rgbImage(rgbROI);
  cv::resize(rgbCropped, rgbCropped, cv::Size(720, 480), cv::INTER_CUBIC);
 
  std::vector<int> compression_params;
  compression_params.push_back(cv::IMWRITE_JPEG_QUALITY); // compression technique
  compression_params.push_back(98); // quality level
  
  cv::imwrite(destination, rgbCropped, compression_params);

  return 0;
}
