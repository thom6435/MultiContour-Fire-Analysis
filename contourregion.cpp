/*
 *    class: ContourRegion
 * modified: 25-August-2020
 *  purpose:
 *
 */

#include "contourregion.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


ContourRegion::ContourRegion() {
    this->area = 0;
    this->pointCount = 0;
    this->index = 0;
    this->analysisExcluded = false;
}

ContourRegion::ContourRegion(int pointCount)
{
    this->fireline = (point2D *) malloc(sizeof (point2D) * pointCount);
    this->area = 0;
    this->pointCount = pointCount;
    this->index = 0;
    this->analysisExcluded = false;
}

void ContourRegion::copy(point2D *p, float area) {
    this->area = area;

    for (int i=0; i< pointCount; i++) {
       this->fireline[i].x = p[i].x;
       this->fireline[i].y = p[i].y;
    }
} // copy


void ContourRegion::setArea(float area){
    this->area = area;
}


void ContourRegion::setEncapsulationFlag(bool flag){
    this->encapsulated = flag;
}


int ContourRegion::getPointCount(){
    return this->pointCount;
}


float ContourRegion::getArea(){
    return this->area;
}


void ContourRegion::addPoint(int x, int y){
    fireline[index].x = x;
    fireline[index].y = y;
    index++;
}

/*
 *  method: getPoint
 * purpose: returns points in a continuous series, independent of the
 *          ends of the point vector.
 */
point2D ContourRegion::getPoint(int index){
    if (index < 0) {
        index = pointCount + index;
    } else if (index>=pointCount) {
        index = index-pointCount;
    }

    return fireline[index];
}

point2D* ContourRegion::getAllPoints() {
    return fireline;
}

point2D ContourRegion::getCenter(){
    return this->center;
}

void ContourRegion::setCenter(point2D cntr) {
    this->center = cntr;
}

point2D ContourRegion::CalculateCenterPoint(){
    int i;
    float X=0.0, Y=0.0;

    for (i=0; i< this->pointCount; i++){
        X += fireline[i].x;
        Y += fireline[i].y;
    } // i

    //puts("CENTERPOINT CALCULATION");
    //printf("X sum = %f  Y Sum = %f\n", X, Y);
    //printf("Total Points: %d\n", this->pointCount);

    this->center.x = (int) X/this->pointCount;
    this->center.y = (int) Y/this->pointCount;

    //printf("center: X: %d   Y: %d\n", this->center.x, this->center.y);

    return this->center;
} // CalculateCenterPoint


/*
 *  method: getNearestPointIndex
 * purpose: Finds the point in the contourRegion that is the nearest to the point
 *          specified in the argument list. Returns -1 if outer region contains no
 *          no points or a nearest point could not be determined (unlikely).
 */
int ContourRegion::getNearestPointIndex(point2D p){
   int index=-1, i;
   double distance, bestDistance = 1000000;


   for (i=0; i<this->pointCount; i++) {
       point2D firelinePt = this->getPoint(i);
       distance = sqrt(pow((double) (p.y - firelinePt.y),2) +
                       pow((double) (p.x - firelinePt.x), 2));
       if (distance < bestDistance) {
           index = i;
           bestDistance = distance;
       }
   } // i

   return index;
}// getNearestPointIndex


/*
 *  method: EncapsulationTest
 * purpuse: Determines if the ContourRegion in the argument is contained by this
 *          ContourRegion.
 *    note: method simply determines if the center point of the ContourRegion argument
 *          is within the min/max x/y bounds of this contourRegion.
 * Warning: This is a simplistic approach that can fail in some conditions. Esp. with
 *          convoluted outer or container contours. A more robust method might be needed
 *          here.
 */
bool ContourRegion::EncapsulationTest(ContourRegion cr){
    int minX=1000000, minY=1000000, maxX=-1, maxY=-1;
    int i;
    point2D crCenter = cr.getCenter();

    for (i=0; i<this->pointCount; i++) {
        if (fireline[i].x < minX)
            minX = fireline[i].x;
        if (fireline[i].y < minY)
            minY = fireline[i].y;
        if (fireline[i].x > maxX)
            maxX = fireline[i].x;
        if (fireline[i].y > maxY)
            maxY = fireline[i].y;
    } // i

    //printf("EncapTest:: minX %d  minY %d  maxX %d  maxY %d\n", minX, minY, maxX, maxY);
    //printf("center: %d,%d\n", crCenter.x, crCenter.y);


    if (crCenter.x >= minX && crCenter.x <= maxX &&
        crCenter.y >= minY && crCenter.y <= maxY)
        return true;

    //printf("not encapsulated\n");
    return false;
}

bool ContourRegion::isAnalysisExcluded() {
    return this->analysisExcluded;
}


void ContourRegion::setAnalysisExcluded(bool b) {
    this->analysisExcluded = b;
}


bool ContourRegion::isEncapsulated(){
    return this->encapsulated;
}

void ContourRegion::PrintCoordinates()
{
    int i;

    printf("Number of Points in Region: %d\n", pointCount);
    for (i=0;i<pointCount; i++)
        printf("%3d: %4d,%4d\n", i, fireline[i].x, fireline[i].y);

}

void ContourRegion::PrintSummary(char *message) {
    printf("Contour Region Summary <<%s>>:\n", message);
    printf("Area: %f\n", area);
    printf("pointCount: %d\n", pointCount);
    printf("Loading Index: %d\n", index);
    if (encapsulated)
        puts("Region is encapsulated");
    else
        puts("Region is NOT Encapsulated");
    if (this->analysisExcluded)
        puts("Region EXCLUDED from processing");

    printf("Region CenterPoint: %d, %d\n", center.x, center.y);
    puts("First region coordinates");

    int reportLen = 10;
    if (reportLen > pointCount)
        reportLen = pointCount;
    for (int i = 0; i<reportLen; i++)
        printf("%d:: %d,%d\n", i, fireline[i].x, fireline[i].y);
} // printSummary
