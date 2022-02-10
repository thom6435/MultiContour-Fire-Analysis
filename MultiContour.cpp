/*
 *  program: MultiContour
 * modified: 10-February-2021
 *   author: Ed Thomas
 *  purpose: Experiment to examine potential of using OpenCV multi-contour image
 *           processing methods to analyze fire imagery.
 *  Library: OpenCV 4.2
 *
 */

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <string.h>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "MultiContour.hpp"
#include "fireline.h"
#include "contourregion.h"
#include "normal.h"
#include "mca_config.h"
#include "imageanalysis.h"
#include "radiometer.h"
#include "thermalArea.h"


int main (int argc, char *argv[]) {
    char aerialFilePath[MAX_AVG_SET][PATH_SIZE], infraRedFilePath[MAX_AVG_SET][PATH_SIZE];
    char grayscaleFilePath[PATH_SIZE], flameFrontPath[PATH_SIZE];
    char radiometerAOIPath[PATH_SIZE];
    char contourProgressionPath[2][PATH_SIZE];
    cv::Mat contourProgression[2], FireLineProgression;
    cv::Mat channel[3], IRchannel[3];
    int i, id, SeriesCount;
    int regionCount, areaPointCount;
    double threshold, area, TotalArea, radiometerFireArea;
    double alpha = 1.25; /* Simple contrast control */
    int beta = 4;        /* Simple brightness control */
    int regionID;
    int Imagery_X_Res, Imagery_Y_Res;
    float medianThermalBackground, meanThermalBackground;
    float minThermalBackground, maxThermalBackground;
    float globalMinimalMedian = 1000, tempMedian;
    double pixels2Meters;
    cv::Mat tripleRed;

    // Vectors containing summarized flame front data.
    std::vector<point2F> fireCenter;
    std::vector<Vector2F> fireSpreadDirection;
    std::vector<float> fireSpreadMagnitude;
    std::vector<float> flameFrontDepth;
    std::vector<float> medianThermal, meanThermal;
    std::vector<int> minThermal, maxThermal;

    // Vectors containing Radiometer data
    std::vector<float> radiometerFlameWidth;
    std::vector<float> radiometerFlameDepth;
    std::vector<float> radiometerA; // Burning area in radiometer AOI
    std::vector<float> radiometerB; // Non-burning area in radiometer AOI
    std::vector<float> fractionalFireArea;

    cv::Scalar color[MAXCONTOURS];
    color[0] = BLUEBLACK;
    color[1] = VIOLETBLUE4;

    if (argc < 2) {
        printf("Usage MultiContour <ConfigurationFile.MCA>\n");
        return 0;
    }

    MCA_Config config;
    config.Read_MCA_Config(argv[1]);
    //config.Read_MCA_Config("/home/thomas/Desktop/Fire-Imagery/Config-4_19_2017Cam6.MCA");
    config.printConfig();

    radiometer RMeterAOI(config.getRadiometerFOV());
    double radiometerArea = RMeterAOI.getArea();
    std::vector<std::vector<cv::Point>> radiometerAOIContour;
    thermalArea thermalAOI;
    cv::Mat radiometerMask;
    cv::Mat thermalMask;
    cv::Mat radiometerAOI;
    cv::Mat radiometerAOIThreshold;

    char reportPath[PATH_SIZE];
    strcpy(reportPath,  config.getOutputDirectory());
    strcat(reportPath, "/");
    strcat(reportPath, config.getAnalysisName());
    strcat(reportPath, "-report.csv");

    char progressionMapPath[PATH_SIZE];
    strcpy(progressionMapPath,  config.getOutputDirectory());
    strcat(progressionMapPath, "/");
    strcat(progressionMapPath, config.getAnalysisName());
    strcat(progressionMapPath, ".jpg");

    char progressionVectorPath[PATH_SIZE];
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
     */ //for (id=0; id < 70; id++) {
    for (id=0; id < SeriesCount-(MAX_AVG_SET-1); id++) {
        char idBuffer[12];
        sprintf(idBuffer,"-%d-",id);
        printf("SEQUENCE %d -----------------------------------------------RFD---\n", id);

        /*
         * Build the Input / Output File Names
         */
        for (i=0; i<MAX_AVG_SET; i++) {
           strcpy(aerialFilePath[i],   config.getInputDirectory());
           strcat(aerialFilePath[i],   "/");
           strcat(aerialFilePath[i],   config.getRGB_Input(id+i));
	  
           strcpy(infraRedFilePath[i], config.getInputDirectory());
           strcat(infraRedFilePath[i], "/");
           strcat(infraRedFilePath[i], config.getIR_Input(id+i));
        } // i

        strcpy(grayscaleFilePath, config.getOutputDirectory());
        strcat(grayscaleFilePath, "/");
        strcat(grayscaleFilePath, config.getAnalysisName());
        strcat(grayscaleFilePath, idBuffer);
        strcat(grayscaleFilePath, "grayscale.jpg");

        strcpy(flameFrontPath, config.getOutputDirectory());
        strcat(flameFrontPath, "/");
        strcat(flameFrontPath, config.getAnalysisName());
        strcat(flameFrontPath, idBuffer);
        strcat(flameFrontPath, "flame-front.jpg");

        strcpy(radiometerAOIPath, config.getOutputDirectory());
        strcat(radiometerAOIPath, "/");
        strcat(radiometerAOIPath, config.getAnalysisName());
        strcat(radiometerAOIPath, idBuffer);
        strcat(radiometerAOIPath, "radAOI.jpg");

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

        /*
         * Load the aerial imagery set
         */
        for (i=0; i<MAX_AVG_SET; i++) {
            aerial[i] = cv::imread(aerialFilePath[i], 1);
            if (!aerial[i].data) {
                puts("Could not open or find the aerial imagery!\n");
                return -1;
            }
        }
        Imagery_X_Res = aerial[0].cols;
        Imagery_Y_Res = aerial[0].rows;
        ia.setImageWidth(Imagery_X_Res);
        ia.setImageHeight(Imagery_Y_Res);






        /*
         * Setup the Radiometer area of Interest
         */
        if (!RMeterAOI.isInitialized()) {
            RMeterAOI.setupAOI(Imagery_X_Res, Imagery_Y_Res, config.getFOV());

            // setup an image mask for example display and potentially use for area calculation.
            radiometerMask = aerial[0].clone();

            for (int y = 0; y < aerial[0].rows; y++ )
                for (int x = 0; x < aerial[0].cols; x++ ) {
                    if (RMeterAOI.isInAOI(x,y)) {
                        radiometerMask.at<cv::Vec3b>(y, x)[0] = 255;
                        radiometerMask.at<cv::Vec3b>(y, x)[1] = 255;
                        radiometerMask.at<cv::Vec3b>(y, x)[2] = 255;
                    } else {
                        radiometerMask.at<cv::Vec3b>(y, x)[0] = 0;
                        radiometerMask.at<cv::Vec3b>(y, x)[1] = 0;
                        radiometerMask.at<cv::Vec3b>(y, x)[2] = 0;
                    }
                }

            //radiometerAOI = radiometerMask.clone();
            cv::Mat radiometerAOIThreshold;

            cv::cvtColor(radiometerMask, radiometerAOIThreshold,cv::COLOR_BGR2GRAY);
            cv::threshold(radiometerAOIThreshold, radiometerAOIThreshold, 0, 255, 0);
            cv::findContours(radiometerAOIThreshold, radiometerAOIContour, cv::RETR_LIST, cv::CHAIN_APPROX_NONE );

            printf("Number of contours in radiometerAOIContour: %d\n", radiometerAOIContour.size());

            //cv::drawContours(aerial[0]/*radiometerAOIThreshold*/, radiometerAOIContour, 0, RED2);
            //cv::imshow("radiometer AOI Contour", aerial[0]/*radiometerAOIThreshold*/);
            //cv::waitKey();
        } // radiometer AOI initialized?

        if (!thermalAOI.isInitialized())
            thermalAOI.setupAOI(Imagery_X_Res, Imagery_Y_Res, config.getFOV());

        /*
         * Load the thermal imagery set.
         */
        for (i=0; i<MAX_AVG_SET; i++) {
            infraRed[i] = cv::imread(infraRedFilePath[i], 1);
            if (id == 0 && i == 0) {
                medianThermalBackground = ia.getMedianBackground(infraRed[0], thermalAOI);
                meanThermalBackground = ia.getMeanBackground(infraRed[0], thermalAOI);
                minThermalBackground = ia.getMinBackground(infraRed[0], thermalAOI);
                maxThermalBackground = ia.getMaxBackground(infraRed[0], thermalAOI);
            }

            if (!infraRed[i].data) {
                puts("Could not open or find the infra-red imagery!\n");
                return -1;
            }

            float medianThermal = ia.getMedianBackground(infraRed[i], thermalAOI);
            float meanThermal = ia.getMeanBackground(infraRed[i], thermalAOI);
            int minThermal = ia.getMinBackground(infraRed[i], thermalAOI);
            int maxThermal = ia.getMaxBackground(infraRed[i], thermalAOI);

            //printf ("::: F:%d I:%d  %f  %f  %d  %d\n", id, i, medianThermal, meanThermal, minThermal, maxThermal);
        } // i


        if (id ==0)
          FireLineProgression = aerial[0].clone();

        contourProgression[id%2] = aerial[0].clone();
	    cv::Mat aerialImage = aerial[0].clone();
	    cv::Mat infraRedImage = infraRed[0].clone();

	    //cv::add(aerial[0], aerial[1], aerialImage);
	    //cv::add(infraRed[0], infraRed[1], infraRedImage);

	    //float gamma = 1.0 / MAX_AVG_SET;
	
	    cv::addWeighted(aerial[0], 0.5, aerial[1], 0.5, 0.0, aerialImage, aerial[0].depth());
	    cv::addWeighted(aerial[2], 0.33,  aerialImage, 0.67, 0.0, aerialImage, aerial[2].depth());

	    cv::addWeighted(infraRed[0], 0.5, infraRed[1], 0.5, 0.0, infraRedImage, infraRed[0].depth());
	    cv::addWeighted(infraRed[2], 0.33, infraRedImage, 0.67, 0.0, infraRedImage, infraRed[2].depth());

	    /*
	    for (i=2; i<MAX_AVG_SET; i++) { // was i=2 and above code enabled
	       cv::addWeighted(aerial[i], gamma,  aerialImage, i*gamma, 0.0, aerialImage, aerial[i].depth());
	       cv::addWeighted(infraRed[i], gamma, infraRedImage, i*gamma, 0.0, infraRedImage, infraRed[i].depth());
	    }
        */

	    cv::Mat middleAerial = aerialImage.clone();
	    int adjustedRed;
        for (int y = 0; y < middleAerial.rows; y++ )
            for (int x = 0; x < middleAerial.cols; x++ ) {
                middleAerial.at<cv::Vec3b>(y, x)[0] *= 0.400;  // scale back blue channel
                middleAerial.at<cv::Vec3b>(y, x)[1] *= 0.700;  // green
                adjustedRed = middleAerial.at<cv::Vec3b>(y, x)[2];
                adjustedRed *= 1.10;
                if (adjustedRed > 255) {
                    //printf("%d,%d (%d)\n", x, y, adjustedRed);
                    adjustedRed = 255;
                }
                middleAerial.at<cv::Vec3b>(y, x)[2] = adjustedRed;  // enhance red
            }

	    meanThermal.push_back(ia.getMeanBackground(infraRedImage, thermalAOI));
	    minThermal.push_back(ia.getMinBackground(infraRedImage, thermalAOI));
	    maxThermal.push_back(ia.getMaxBackground(infraRedImage, thermalAOI));
	    tempMedian = ia.getMedianBackground(infraRedImage, thermalAOI);
	    medianThermal.push_back(tempMedian);

	    if (tempMedian < globalMinimalMedian)
             globalMinimalMedian = tempMedian;

	    //cv::namedWindow("Aerial Merged", cv::WINDOW_AUTOSIZE );
	    //cv::imshow("Aerial Merged", aerialImage);
	    //cv::namedWindow("infraRed Merged", cv::WINDOW_AUTOSIZE );
	    //cv::imshow("infraRed Merged", infraRedImage);

        cv::Mat redImage = aerialImage.clone();
        printf("Extracting Red channel from Aerial Image...\n");
        split(aerialImage, channel);
        channel[0]=cv::Mat::zeros(aerialImage.rows, aerialImage.cols, CV_8UC1);//Set blue channel to 0
        channel[1]=cv::Mat::zeros(aerialImage.rows, aerialImage.cols, CV_8UC1);//Set green channel to 0
        merge(channel,3,redImage);// Create the image from the red channel and the nullified blue & green channels

        /*
         * TripleRed Image Experiment
         */
        //cv::Mat
        tripleRed = redImage.clone();
        insertChannel(channel[2], tripleRed, 0);
        insertChannel(channel[2], tripleRed, 1);

        cv::erode(tripleRed, tripleRed, getStructuringElement(cv::MORPH_ELLIPSE,
                                                              cv::Size(9, 7),
                                                              cv::Point(3, 3)));
        cv::dilate(tripleRed, tripleRed, getStructuringElement(cv::MORPH_ELLIPSE,
                                                            cv::Size(13, 11),
                                                            cv::Point(5, 5)));
        GaussianBlur(tripleRed, tripleRed, cv::Size(11,11), 1);
        GaussianBlur(tripleRed, tripleRed, cv::Size(5,5), 1);
        GaussianBlur(tripleRed, tripleRed, cv::Size(3,3), 1);

        cv::erode(tripleRed, tripleRed, getStructuringElement(cv::MORPH_ELLIPSE,
                                                              cv::Size(7, 5),
                                                              cv::Point(2, 2)));
        cv::dilate(tripleRed, tripleRed, getStructuringElement(cv::MORPH_ELLIPSE,
                                                            cv::Size(9, 7),
                                                            cv::Point(3, 3)));
        // next two operations added 5-Jan-2021
        cv::erode(tripleRed, tripleRed, getStructuringElement(cv::MORPH_ELLIPSE,
                                                              cv::Size(11, 9),
                                                              cv::Point(2, 2)));
        cv::dilate(tripleRed, tripleRed, getStructuringElement(cv::MORPH_ELLIPSE,
                                                            cv::Size(9, 11),
                                                            cv::Point(3, 3)));

        GaussianBlur(tripleRed, tripleRed, cv::Size(11,11), 1);
        GaussianBlur(tripleRed, tripleRed, cv::Size(9,9), 1); // 1-5-21
        GaussianBlur(tripleRed, tripleRed, cv::Size(5,5), 1); // 1-5-21

        tempMedian *= .50;
        tripleRed = ia.MaskImage(tripleRed, infraRedImage, tempMedian);

        tripleRed = ia.setContrastBalance((double)1.7, (int) -20, tripleRed);

        //cv::imshow("tripleRed Image, processed", tripleRed);

        cv::Mat TR_thresholdImage = tripleRed.clone();
        cv::Mat TR_grayScale = tripleRed.clone();
        cvtColor(tripleRed, TR_grayScale, cv::COLOR_BGR2GRAY);
        cv::threshold(TR_grayScale, TR_thresholdImage, FIRELINE_MIN_THRESHOLD, FIRELINE_MAX_THRESHOLD, 0); // 0, 100, 0

        std::vector<std::vector<cv::Point>> TR_contours;
        cv::findContours( TR_thresholdImage, TR_contours, cv::RETR_LIST, cv::CHAIN_APPROX_SIMPLE);  // alt method RETR_TREE, CHAIN_APPROX_NONE, CHAIN_APPROX_SIMPLE

        for (regionID = 0; regionID < TR_contours.size(); regionID++) {// loop through contour regions
            cv::drawContours(tripleRed, TR_contours, regionID, RED1);
            cv::drawContours(middleAerial, TR_contours, regionID, RED1);
        }
        //cv::namedWindow("tripleRed Image, processed", cv::WINDOW_AUTOSIZE );
        //cv::imshow("tripleRed Image, processed", middleAerial);

        cv::Mat infraRedBGR = redImage.clone();
        split(infraRedImage, IRchannel);
        insertChannel(IRchannel[2], infraRedBGR, 0); //Insert IR RedChannel as blue channel
        insertChannel(IRchannel[2], infraRedBGR, 1); //Insert IR RedChannel as green channel

        cv::Mat radiometerIRBGR = infraRedBGR.clone();

        //cv::imshow("radiometerIRBGR A", radiometerIRBGR);
        //cv::imshow("infraRedBGR Image", infraRedBGR);

        for (int y = 0; y < infraRedBGR.rows; y++ )
            for (int x = 0; x < infraRedBGR.cols; x++ )
                for (int c = 0; c < infraRedBGR.channels(); c++ )
                    infraRedBGR.at<cv::Vec3b>(y,x)[c] =
                        cv::saturate_cast<uchar>( alpha*infraRedBGR.at<cv::Vec3b>(y,x)[c] + beta );

        infraRedBGR = ia.MaskImage(infraRedBGR, 125); // was 100
        cv::blur(infraRedBGR, infraRedBGR, cv::Size(3,3));

        cv::Mat grayScaleImg = infraRedBGR.clone();
        cvtColor(infraRedBGR, grayScaleImg, cv::COLOR_BGR2GRAY);
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


        /**************************************************************************
         ********************* Process the radiometer pixel AOI *******************
         *************************************************************************/
        int radRegionCount;
        radiometerIRBGR = ia.logicalAndImage(RMeterAOI, radiometerIRBGR);

        cv::Mat radiometerClipped = TR_grayScale.clone();
        cv::Mat radiometerThreshold = TR_grayScale.clone();
        radiometerClipped = ia.logicalAndImage(RMeterAOI, radiometerClipped);

        cv::threshold(radiometerClipped, radiometerThreshold, 200, 255, 0); // 150, 255, 0
        std::vector<std::vector<cv::Point> > radiometerContours;
        cv::findContours(radiometerThreshold, radiometerContours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE );
        radRegionCount = radiometerContours.size();
        ContourRegion radCR[radRegionCount];

        radiometerFireArea = 0.0;
        for (regionID = 0; regionID < radRegionCount; regionID++) { // loop through contour regions
            const std::vector<cv::Point> &c = radiometerContours[regionID];

            areaPointCount = (int) c.size();
            if (areaPointCount > 50 &&
                areaPointCount < 1048576) {
                cv::drawContours(radiometerIRBGR, radiometerContours, regionID, BLUE2);
                radCR[regionID] = ContourRegion(areaPointCount);
                area = cv::contourArea(c, false);

                for (size_t j = 0; j < c.size(); j++)
                    radCR[regionID].addPoint(c[j].x, c[j].y);

                radCR[regionID].setArea(area);
                radCR[regionID].CalculateCenterPoint();
                radiometerFireArea += area;
            }
        } // idx < contours.size

        pixels2Meters = Imagery_X_Res / config.getFOV();
        radiometerFireArea = radiometerFireArea / pow(pixels2Meters,2.0);
        printf("Total Radiometer Fire Area: %f of %f (sq. pixels)\n", radiometerFireArea, radiometerArea);

        radiometerA.push_back(radiometerFireArea);
        radiometerB.push_back(radiometerArea - radiometerFireArea);
        fractionalFireArea.push_back(radiometerFireArea/radiometerArea);

        //cv::imshow("radiometerIRBGR", radiometerIRBGR);
        //cv::imshow("radiometer mask", radiometerMask);
        //cv::imshow("radiometer thresholded", radiometerThreshold);


        /**************************************************************************
         ************* Perform fire-line analysis using entire scene **************
         **************************************************************************/

        cv::Mat thresholdImage = grayScaleImg.clone();
        //for (i=0; i<MAXCONTOURS; i++) {
        i=0;
        //cv::threshold(grayScaleImg, thresholdImage, threshold, threshold+CONTOUR_INTERVAL, 0);
        cv::threshold(grayScaleImg, thresholdImage, 0, 100, 0); // 0, 100, 0
        std::vector<std::vector<cv::Point> > contours;
        cv::findContours( thresholdImage, contours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE );  // alt method RETR_TREE, CHAIN_APPROX_SIMPLE
        TotalArea = 0.0;
        printf("contour level: %d  Number of Regions %ld\n", i, contours.size());

        regionCount = contours.size();
        ContourRegion cr[regionCount];
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

        if (id > 0) { // Normal creation and analysis can be done starting with the second image
            prevFireLine.progressionAnalysis(currentFireLine);
            std::vector<point3F> p;
            std::vector<Vector3F> d;
            std::vector<float> magnitude;

            int count = prevFireLine.getNormalCount();

            //printf("LEVEL: %d  normal Count: %d\n", id, count);

            for (int i = 0; i < count; i++) {
                Ray3F n = prevFireLine.getNormal(i);

                if (n.Pt.x < Imagery_X_Res - 1 && n.Pt.x > 0 &&  // need to push this filtering down in sub classes
                    n.Pt.y < Imagery_Y_Res - 1 && n.Pt.y > 0)    // filter out normals that originate on image edge.
                {
                    cv::Point A, B;
                    A.x = n.Pt.x;
                    A.y = n.Pt.y;
                    B.x = A.x + n.D.x;
                    B.y = A.y + n.D.y;

                    float m = sqrt(pow(B.x - A.x, 2) + pow(B.y - A.y, 2));
                    if (m < ABNORMAL_RAY_THRESHOLD) {  //Filter out abnormally long rays
                        // printf("xxx>I::%d mag: %f\n", i, m);
                        p.push_back(n.Pt);
                        d.push_back(n.D);
                        magnitude.push_back(m);

                        cv::arrowedLine(contourImage, A, B, RED2, 1, 8, 0, 0.1);
                        cv::arrowedLine(FireLineProgression, A, B, RED2, 1, 8, 0, 0.1);
                        cv::arrowedLine(contourProgression[(id - 1) % 2], A, B, RED2, 1, 8, 0, 0.1);
                    } // magnitude within normal range
                } // normal origination not on image edge
            } // i

            //cv::imshow("ContourImage", contourImage);

            /*
             * Process and summarize the fire line progression vectors
             */
            if (p.size() > 10) {
                /**
                 * Use the normal vector origin points to fit a line for use in determining
                 * flame front depth.
                 */
                cv::Point A, B;
                ia.FireLine_LineFit(p, &A, &B);
                cv::line(tripleRed, A, B, VIOLET1, 2);
                cv::line(middleAerial, A, B, VIOLET1, 2);

                /**
                 * Process the radiometer AOI to determine flame front width and depth
                 * Use A & B as starting points and find the intersection points of the line
                 * with the contour defined by radiometerContours
                 */
                cv::Point radA=A, radB=B;

                std::vector<cv::Point> radInterceptPt1, radInterceptPt2;
                ia.optimizeFlameFrontLine(&radA, &radB, radiometerContours,
                                          Imagery_X_Res, Imagery_Y_Res);
                ia.clipFrameFrontLine(&radA, &radB, radiometerAOIContour);
                ia.clipFrameFrontLine(&radA, &radB, radiometerContours);

                double RFD = ia.getFlameFrontDepth(radA, radB, radiometerContours,
                                                   &radInterceptPt1, &radInterceptPt2,
                                                   Imagery_X_Res, Imagery_Y_Res);
                double RFW = sqrt(pow(radA.x - radB.x, 2.0) + pow(radA.y - radB.y, 2.0));

                RFD = RFD / pixels2Meters;
                RFW = RFW / pixels2Meters;

                if (RFD < 0) {
                    RFD = -1.0;
                    RFW = -1.0;
                } else {
                    cv::line(radiometerIRBGR, radA, radB, BLUE2, 2);
                }
                radiometerFlameDepth.push_back(RFD);
                radiometerFlameWidth.push_back(RFW);

                // plot flame depth analysis back to radiometer AOI
                for (int k=0; k<radInterceptPt1.size(); k++)
                    cv::line(radiometerIRBGR, radInterceptPt1[k], radInterceptPt2[k], GREEN2, 1);

                //cv::imshow("radiometerIRBGR", radiometerIRBGR);

                /**
                 * Process the full frame imagery to  analyze fire line width and depth.
                 */

                std::vector<cv::Point> InterceptPt1, InterceptPt2;
                ia.optimizeFlameFrontLine(&A, &B, TR_contours, Imagery_X_Res, Imagery_Y_Res);

                cv::line(middleAerial, A, B, BLUE2, 2);
                flameFrontDepth.push_back(ia.getFlameFrontDepth(A, B, TR_contours,
                                                       &InterceptPt1, &InterceptPt2,
                                                       Imagery_X_Res, Imagery_Y_Res));

                for (int k=0; k<InterceptPt1.size(); k++) {
                    cv::line(tripleRed, InterceptPt1[k], InterceptPt2[k], ORANGE2, 1);
                    cv::line(middleAerial, InterceptPt1[k], InterceptPt2[k], ORANGE2, 1);
                }

                //cv::namedWindow("tripleRed Image, processed", cv::WINDOW_AUTOSIZE );
                //cv::imshow("tripleRed Image, processed", middleAerial);

                ia.AnalyzeProgressionVector(p, d, magnitude);
                fireCenter.push_back(ia.getFlameFrontCenter());
                fireSpreadDirection.push_back(ia.getFlameFrontDirection());
                fireSpreadMagnitude.push_back(ia.getFlameFrontSpreadMagnitude());
            } else {
                point2F nullCenter;
                nullCenter.x = -1;
                nullCenter.y = -1;
                fireCenter.push_back(nullCenter);
                Vector2F nullDirection;
                nullDirection.x = -1;
                nullDirection.y = -1;
                flameFrontDepth.push_back(-1.0);
                fireSpreadDirection.push_back(nullDirection);
                fireSpreadMagnitude.push_back(-1.0);
                radiometerFlameDepth.push_back(-1.0);
                radiometerFlameWidth.push_back(-1.0);
            }
        } // at least one previous image has been processed.

        //cv::waitKey(0);

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

        //cv::imwrite(grayscaleFilePath, grayScaleImg);
        cv::imwrite(flameFrontPath, middleAerial);
        //cv::imwrite(radiometerAOIPath, radiometerIRBGR);
        cv::imwrite(radiometerAOIPath, tripleRed);
        cv::imwrite(progressionMapPath, FireLineProgression);

        if (id > 0)
            cv::imwrite(contourProgressionPath[(id-1)%2], contourProgression[(id-1)%2]);

    } // i loop through all images in series

    ia.printProgressionSummary(fireCenter, fireSpreadDirection, fireSpreadMagnitude,
                               flameFrontDepth, meanThermal, medianThermal,
                               minThermal, maxThermal, radiometerFlameWidth, radiometerFlameDepth,
                               radiometerA, radiometerB, fractionalFireArea,
                               config, Imagery_X_Res, Imagery_Y_Res);

    ia.printProgressionSummary(reportPath, fireCenter, fireSpreadDirection, fireSpreadMagnitude,
                               flameFrontDepth, meanThermal, medianThermal,
                               minThermal, maxThermal,  radiometerFlameWidth, radiometerFlameDepth,
                               radiometerA, radiometerB, fractionalFireArea,
                               config, Imagery_X_Res, Imagery_Y_Res);


    fclose(VectorTxtStream);
    return 0;
} //contourAnalysis main