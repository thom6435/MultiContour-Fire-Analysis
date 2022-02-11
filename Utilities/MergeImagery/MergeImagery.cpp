/*
 *  program: MergeImagery
 * modified: 9-September-2020
 *   author: Ed Thomas
 *  purpose: Experiment to examine potential of using OpenCV multi-contour image
 *           processing methods to analyze fire imagery.
 *  Library: OpenCV 4.2
 *
 */

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/xphoto.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core/cvstd.hpp>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "MergeImagery.hpp"
#include "fireline.h"
#include "contourregion.h"
#include "normal.h"
#include "mca_config.h"
#include "imageanalysis.h"


int main (int argc, char *argv[]) {
  char aerialFilePath[MAX_AVG_SET][PATH_SIZE], infraRedFilePath[MAX_AVG_SET][PATH_SIZE];
  char contourFilePath[PATH_SIZE], grayscaleFilePath[PATH_SIZE];
  char contourProgressionPath[2][PATH_SIZE];
  cv::Mat contourProgression[2];
  cv::Mat channel[3];
  cv::Mat IRchannel[3];
  int i, id, SeriesCount;
  double threshold, area, TotalArea;
  double alpha = 1.25; /* Simple contrast control */
  int beta = 4;        /* Simple brightness control */
  int regionID;
  int Imagery_X_Res, Imagery_Y_Res;
  cv::Mat FireLineProgression;
  
  cv::Scalar color[MAXCONTOURS];
  color[0] = BLUEBLACK;
  color[1] = VIOLETBLUE4;

    //if (argc < 2) {
    //  printf("Usage MultiContour <ConfigurationFile.MCA>\n");
    //  return 0;
    //}

    MCA_Config config;
    //config.Read_MCA_Config(argv[1]);
    char configPath[256];
    strcpy(configPath, "/home/thomas/Desktop/Fire-Imagery/config.MCA");
    config.Read_MCA_Config(configPath);
    config.printConfig();

    char progressionMapPath[PATH_SIZE], progressionVectorPath[PATH_SIZE];

    strcpy(progressionMapPath,  config.getOutputDirectory());
    strcat(progressionMapPath, "/");
    strcat(progressionMapPath, config.getAnalysisName());
    strcat(progressionMapPath, ".jpg");

    strcpy(progressionVectorPath,  config.getOutputDirectory());
    strcat(progressionVectorPath, "/");
    strcat(progressionVectorPath, config.getAnalysisName());
    strcat(progressionVectorPath, ".jpg");

    FILE *VectorTxtStream = fopen(progressionVectorPath, "w");
    if (!VectorTxtStream) {
        printf("Error unable to open output summary/analysis file:\n\t%s\n", progressionVectorPath);
        exit(0);
    }

    ImageAnalysis ia;
    FireLine prevFireLine, currentFireLine;
    ContourRegion prevRegion[MAXREGIONS];

    SeriesCount = config.getFileCount();

    /*
     * Loop through the image series
     */
    for (id=0; id < SeriesCount-(MAX_AVG_SET-1); id++) {
        char idBuffer[5];
        sprintf(idBuffer,"-%d-",id);

        /*
         * Build the Input / Ouptput File Names
         */
	for (i=0; i<MAX_AVG_SET; i++) {	
	  strcpy(aerialFilePath[i],   config.getInputDirectory());
	  strcat(aerialFilePath[i],   "/");
	  strcat(aerialFilePath[i],   config.getRGB_Input(id+i));
	  
	  strcpy(infraRedFilePath[i], config.getInputDirectory());
	  strcat(infraRedFilePath[i], "/");
	  strcat(infraRedFilePath[i], config.getIR_Input(id+i));
	}
	
/*
        strcpy(contourFilePath,  config.getOutputDirectory());
        strcat(contourFilePath,  "/");
        strcat(contourFilePath,  config.getAnalysisName());
        strcat(contourFilePath,  idBuffer);
        strcat(contourFilePath,  "contour.jpg");
*/
        strcpy(grayscaleFilePath, config.getOutputDirectory());
        strcat(grayscaleFilePath, "/");
        strcat(grayscaleFilePath, config.getAnalysisName());
        strcat(grayscaleFilePath, idBuffer);
        strcat(grayscaleFilePath, "grayscale.jpg");

        strcpy(contourProgressionPath[id%2], config.getOutputDirectory());
        strcat(contourProgressionPath[id%2], "/");
        strcat(contourProgressionPath[id%2], config.getAnalysisName());
        strcat(contourProgressionPath[id%2], idBuffer);
        strcat(contourProgressionPath[id%2], "progression.jpg");

        printf("  Aerial Imagery:\n");
	for (i=0; i<MAX_AVG_SET; i++)
	  printf("\t\t%s\n",aerialFilePath[i]);
	printf("infraRed Imagery:\n");
	for (i=0; i<MAX_AVG_SET; i++)
	  printf("\t\t%s\n", infraRedFilePath[i]);
	
	cv::Mat infraRed[MAX_AVG_SET];
	cv::Mat aerial[MAX_AVG_SET];
	
	for (i=0; i<MAX_AVG_SET; i++) {        
	  infraRed[i] = cv::imread(infraRedFilePath[i], 1);
	  if (!infraRed[i].data) {
	    puts("Could not open or find the infra-red imagery!\n");
	    return -1;
	  }
	  
	  aerial[i] = cv::imread(aerialFilePath[i], 1);
	  if (!aerial[i].data) {
	    puts("Could not open or find the aerial imagery!\n");
	    return -1;
	  }
	} // i

	// Resize all Imagery to InfraRed imagery specifications
	for (i=0; i<MAX_AVG_SET; i++)
	  cv::resize(aerial[i], aerial[i], infraRed[0].size());
	
	Imagery_X_Res = aerial[0].cols;
        Imagery_Y_Res = aerial[1].rows;
	
	printf("Aerial image resized to match infra-red image dimensions (%d x %d)\n",
	       Imagery_X_Res,Imagery_Y_Res);
	
        if (id ==0)
            FireLineProgression = aerial[0].clone();

        contourProgression[id%2] = aerial[0].clone();

	cv::Mat aerialImage = aerial[0].clone();
	cv::Mat infraRedImage = infraRed[0].clone();

	//cv::add(aerial[0], aerial[1], aerialImage);
	//cv::add(infraRed[0], infraRed[1], infraRedImage);

	float gamma = 1.0 / MAX_AVG_SET;
	
	cv::addWeighted(aerial[0], gamma, aerial[1], gamma, 0.0, aerialImage, aerial[0].depth());
	cv::addWeighted(infraRed[0], gamma, infraRed[1], gamma, 0.0, infraRedImage, infraRed[0].depth());

	
	
	for (i=2; i<MAX_AVG_SET; i++) {
	  cv::addWeighted(aerial[i], gamma,  aerialImage, i*gamma, 0.0, aerialImage, aerial[i].depth());
	  cv::addWeighted(infraRed[i], gamma, infraRedImage, i*gamma, 0.0, infraRedImage, infraRed[i].depth());
	}
	
	cv::namedWindow("Aerial Merged", cv::WINDOW_AUTOSIZE );
	cv::imshow("Aerial Merged", aerialImage);
	cv::namedWindow("infraRed Merged", cv::WINDOW_AUTOSIZE );
	cv::imshow("infraRed Merged", infraRedImage);
	cv::waitKey(0);
	
        cv::Mat redImage = aerialImage.clone();
        printf("Extracting Red channel from Aerial Image...\n");
        split(aerialImage, channel);
        channel[0]=cv::Mat::zeros(aerialImage.rows, aerialImage.cols, CV_8UC1);//Set blue channel to 0
        channel[1]=cv::Mat::zeros(aerialImage.rows, aerialImage.cols, CV_8UC1);//Set green channel to 0
        merge(channel,3,redImage);// Create the image from the red channel and the nullified blue & green channels

        cv::Mat infraRedBGR = redImage.clone();
        split(infraRedImage, IRchannel);
        insertChannel(IRchannel[2], infraRedBGR, 0); //Insert IR RedChannel as blue channel
        insertChannel(IRchannel[2], infraRedBGR, 1); //Insert IR RedChannel as green channel

        //Mat NDVI = CalculateNDVI(infraRedBGR);
        //namedWindow("NDVI_Image", WINDOW_AUTOSIZE );
        //imshow("NDVI_Image", NDVI);
        //waitKey(0);

        for (int y = 0; y < infraRedBGR.rows; y++ )
            for (int x = 0; x < infraRedBGR.cols; x++ )
                for (int c = 0; c < infraRedBGR.channels(); c++ )
                    infraRedBGR.at<cv::Vec3b>(y,x)[c] =
                        cv::saturate_cast<uchar>( alpha*infraRedBGR.at<cv::Vec3b>(y,x)[c] + beta );

        infraRedBGR = ia.MaskImage(infraRedBGR, 125); // was 100
        //GaussianBlur(infraRedBGR,infraRedBGR, Size(5,5), 1);
        cv::blur(infraRedBGR, infraRedBGR, cv::Size(3,3));
        //cv::imshow("blurred", infraRedBGR);

        //infraRedBGR= Dilation(2, 5, infraRedBGR);
        //imshow("Dilation 2, 5", infraRedBGR);

        cv::Mat grayScaleImg = infraRedBGR.clone();
        cvtColor(infraRedBGR, grayScaleImg, cv::COLOR_BGR2GRAY);
        //Mat contourImage = infraRedBGR.clone();
        cv::Mat contourImage = aerialImage.clone();

        //printHistogram(grayScaleImg);
        equalizeHist(grayScaleImg, grayScaleImg);
       
        threshold = 20;
        int erosion_size = 7, dilation_size;

        cv::Mat erodeElement1 = getStructuringElement(cv::MORPH_ELLIPSE,
                             cv::Size( 2*erosion_size + 1, 2*erosion_size+1 ),
                             cv::Point( erosion_size, erosion_size ) );
        cv::erode(grayScaleImg, grayScaleImg, erodeElement1);

        dilation_size = 5;
        cv::Mat dilationElement1 = getStructuringElement(cv::MORPH_RECT,
                             cv::Size( 2*dilation_size + 1, 2*dilation_size+1 ),
                             cv::Point( dilation_size, dilation_size ) );
        cv::dilate(grayScaleImg, grayScaleImg, dilationElement1);

        dilation_size = 7;
        cv::Mat dilationElement2 = getStructuringElement(cv::MORPH_ELLIPSE,
                             cv::Size( 2*dilation_size + 1, 2*dilation_size+1 ),
                             cv::Point( dilation_size, dilation_size ) );
        cv::dilate(grayScaleImg, grayScaleImg, dilationElement2);
       
        cv::blur(grayScaleImg, grayScaleImg, cv::Size(10, 10));
        cv::blur(grayScaleImg, grayScaleImg, cv::Size(5, 5));

        cv::Mat thresholdImage = grayScaleImg.clone();
        //for (i=0; i<MAXCONTOURS; i++) {
        i=0;
        //cv::threshold(grayScaleImg, thresholdImage, threshold, threshold+CONTOUR_INTERVAL, 0);
        cv::threshold(grayScaleImg, thresholdImage, 0, 100, 0); // 0, 100, 0

        std::vector<std::vector<cv::Point> > contours;
        cv::findContours( thresholdImage, contours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE );  // alt method RETR_TREE, CHAIN_APPROX_SIMPLE
        TotalArea = 0.0;
        printf("contour level: %d  Number of Regions %ld\n", i, contours.size());


        int regionCount = contours.size();
        ContourRegion cr[regionCount];
        int areaPointCount;
        currentFireLine = FireLine(regionCount, Imagery_X_Res, Imagery_Y_Res);

        for (regionID = 0; regionID < regionCount; regionID++) { // loop through contour regions
            cv::drawContours(contourImage, contours, regionID, color[i]);
            cv::drawContours(FireLineProgression, contours, regionID, color[i]);
            cv::drawContours(contourProgression[id%2], contours, regionID, GREEN1);
            if (id > 0)
               cv::drawContours(contourProgression[(id-1)%2], contours, regionID, BLUE1);

            const std::vector<cv::Point>& c = contours[regionID];

            areaPointCount = (int) c.size();
            cr[regionID] = ContourRegion(areaPointCount);
            area = cv::contourArea(c, false);

            for (size_t j=0; j < c.size(); j++)
               cr[regionID].addPoint(c[j].x, c[j].y);

            cr[regionID].setArea(area);
            cr[regionID].CalculateCenterPoint();

            TotalArea += area;
        } // idx < contours.size

        currentFireLine.addContourRegions(cr);
        currentFireLine.setRegionCount(regionCount);
        currentFireLine.filterRegions();

        puts("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv");
        puts("-----------------CURRENT FIRE LINE-----------------------");
        printf("Total Area in contour level: %f\n\n", TotalArea);
        currentFireLine.print();
        puts("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");

        if (id>0) {
            puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
            puts("-----------------PREVIOUS FIRE LINE----------------------");
            prevFireLine.print();
            puts("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM");
        }

        if (id > 0)  { // Normal creation and analysis can be done starting with the second image
            prevFireLine.progressionAnalysis(currentFireLine);

            int normalCount = prevFireLine.getNormalCount();
            for (int i = 0; i<normalCount; i++) {
                Ray3F n = prevFireLine.getNormal(i);
                if (n.Pt.x < Imagery_X_Res-1 && n.Pt.x > 0 &&  // need to push this filtering down in sub classes
                    n.Pt.y < Imagery_Y_Res-1 && n.Pt.y > 0) {
                    cv::Point A, B;
                    A.x = n.Pt.x;
                    A.y = n.Pt.y;
                    B.x = A.x + n.D.x;
                    B.y = A.y + n.D.y;
                    if (sqrt(pow(B.x-A.x,2) + pow(B.y-A.y,2)) < 250) {  //Filter out large rays
                        cv::arrowedLine(contourImage, A, B, RED2, 1, 8, 0, 0.1);
                        cv::arrowedLine(FireLineProgression, A, B, RED2, 1, 8, 0, 0.1);
                        cv::arrowedLine(contourProgression[(id-1)%2], A, B, RED2, 1, 8, 0, 0.1);
                    }
                }
            } // i
        } // at least one prevous image has been processed.

        /*
         * Copy the current fireline to the previous fireline
         */
        prevFireLine = FireLine();
        prevFireLine.setRegionCount(regionCount);
        prevFireLine.setImageX(Imagery_X_Res);
        prevFireLine.setImageY(Imagery_Y_Res);

        for (int k=0; k<regionCount; k++) {
            prevRegion[k] = ContourRegion(cr[k].getPointCount());
            prevRegion[k].setArea(cr[k].getArea());
            prevRegion[k].setEncapsulationFlag(cr[k].isEncapsulated());
            prevRegion[k].setAnalysisExcluded(cr[k].isAnalysisExcluded());
            point2D center = cr[k].getCenter();
            prevRegion[k].setCenter(center);
            prevRegion[k].copy(cr[k].getAllPoints(), cr[k].getArea());
        }
        prevFireLine.addContourRegions(prevRegion);

        threshold += CONTOUR_INTERVAL;
    //}

        //cv::imwrite(contourFilePath, contourImage);
        cv::imwrite(grayscaleFilePath, grayScaleImg);
        cv::imwrite(progressionMapPath, FireLineProgression);
        if (id > 0)
            cv::imwrite(contourProgressionPath[(id-1)%2], contourProgression[(id-1)%2]);
    } // i loop through all images in series


    fclose(VectorTxtStream);

    return 0;
} //contourAnalysis main














