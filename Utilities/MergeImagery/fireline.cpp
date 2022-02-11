#include "fireline.h"
#include <stdio.h>
#include <stdlib.h>


FireLine::FireLine(int RegionCount, int maxX, int maxY)
{
   this->RegionCount = RegionCount;
   this->max_Image_X = maxX;
   this->max_Image_Y = maxY;
   this->region = (ContourRegion *) malloc(sizeof (ContourRegion) * RegionCount);
   this->nvCount = 0;
}


FireLine::FireLine() {
    this->RegionCount = 0;
    this->nvCount = 0;
}


ContourRegion FireLine::initializeContourRegion(int index, int pointCount){
    this->region[index] = ContourRegion(pointCount);

    return this->region[index];
}


ContourRegion FireLine::getContourRegion(int index){
    return this->region[index];
}


void FireLine::addContourRegions(ContourRegion cr[]){
    this->region = cr;
}

point2D FireLine::getPoint2D(int ContourID, int PointID){
    return this->region[ContourID].getPoint(PointID);
}


void FireLine::addNextPoint(int ContourID, float x, float y){
    this->region[ContourID].addPoint(x, y);
}


int FireLine::MinContourPointThreshold() {
    return MINCOORDS;
}

void FireLine::setImageY(int y){
    this->max_Image_Y = y;
}

void FireLine::setImageX(int x){
    this->max_Image_X = x;
}

int FireLine::getImageY(){
    return this->max_Image_Y;
}

int FireLine::getImageX(){
    return this->max_Image_X;
}

/*
 *  method: setRegionCount
 * purpose: Used to update the regionCount once contours containing less than the min.
 *          threshold number of points are encountered.
 */
void FireLine::setRegionCount(int regionCount){
    this->RegionCount = regionCount;
}


int FireLine::getRegionCount(){
    return this->RegionCount;
}


void FireLine::print(){
    printf("FireLine Contains %d  Regions\n", this->RegionCount);

    for (int i=0; i<this->RegionCount; i++) {
        printf("Contour Region: %d\n", i);
        this->region[i].PrintSummary(" ");
    }
} // print


/*
 *   method: progressionAnalysis
 * modified: 18-August-2020
 *  purpose:
 */
void FireLine::progressionAnalysis(FireLine nextFireLine) {
    int i, j, k, index;

    this->nvCount = 0;

    for (i=0; i<nextFireLine.RegionCount; i++) {
        if (!nextFireLine.region[i].isAnalysisExcluded()) {           
            ContourRegion outerRegion = nextFireLine.region[i];
            for (j=0; j<this->RegionCount; j++) {
                if (!this->region[j].isAnalysisExcluded()) {
                if (outerRegion.EncapsulationTest(this->region[j])) {
                    printf("OuterRegion %d encapsulates InnerRegion %d\n", i, j);
                    ContourRegion innerRegion = this->region[j];
                    for (k=0; k<outerRegion.getPointCount(); k+=3) {
                        point2D outerPoint = outerRegion.getPoint(k);
                        if (outerPoint.x != 0 && outerPoint.x != this->max_Image_X &&
                            outerPoint.y != 0 && outerPoint.y != this->max_Image_Y) {

                            index = innerRegion.getNearestPointIndex(outerPoint);

                            point2D innerPoint = innerRegion.getPoint(index);
                            if (innerPoint.x != 0 && innerPoint.x != this->max_Image_X &&
                                innerPoint.y != 0 && innerPoint.y != this->max_Image_Y) {

                                point2D fittingPoints[FIT_POINT_MAX];
                                int VecIndex = 0;
                                for (int p = -3; p < 4; p++) {
                                    fittingPoints[VecIndex] = innerRegion.getPoint(index+p);
                                    VecIndex++;
                                }
                                Normal flNormal = Normal(fittingPoints,
                                                         this->max_Image_X, this->max_Image_Y);
                                flNormal.CalculateNormal(outerPoint);
                                // Determine if the normal intersects the innerregion fire contour.
                                // if so, it is not a valid solution.  If the normal in the opposite direction
                                // still intercepts the inner fire contour. Then no solution exists for that
                                // point and is likely the product of a complext, convoluted fire line.
                                if (flNormal.innerContourIntercept(innerRegion)) {
                                    puts("--------> ERROR Normal intersected originating fire line contour <-----------");
                                } else {
                                if (flNormal.intercept(outerRegion)) {
                                    /* As the normals are saved from the temporary class used to generate them,
                                     * may want to check for multiple intersections from a common point.  Ie, this
                                     * could happen in the case of a very convoluted fire line.  In this case, we want
                                     * to keep the shortest, and discard the longer (ones).
                                     */
//                                   if (nvCount > 10000)
//                                       puts("oops");

                                    this->normals[this->nvCount] = flNormal.getNormal();
                                    this->nvCount++;                                    
                                    flNormal.print();
                                }
                                }
                            } // closest inner point is not on image edge
                        } // OuterRegion point is not on image edge

                        //flNormal.FindNormalRegionIntersect(outerRegion);
                    } // K
                } // InnerRegion is encapsulated within outerRegion
                }  else { // inner region ! analysis excluded
                    puts("INNER REGION ANALYSIS EXCLUDED-------------------------------------");
                }
            } // j
        } // contour in second fireline is not encapsulated
    } // i
}// progressionAnalysis


/*
 *   method: filterRegions
 * modified: 28-August-2020
 *  purpose: Filters out the regions that compose the fireline and flags regions that:
 *             1) are encapsulated or enclosed inside a larger region
 *             2) are small in area.  How small needs to be defined by Matt.
 *           The main reason for doing this is to avoid meaningless analysis and progression
 *           direction graphs.  The small, encapsulated areas come and go and confuse the
 *           overall approach.
 *     note: Other alternatives tried included different thresholding limits for the contours
 *           and different blur settings.
 */
void FireLine::filterRegions() {
    int i, j;

    for (i=0; i<this->RegionCount; i++)
        if (region[i].getArea() <= MINAREA) {
            region[i].setAnalysisExcluded(true);
            printf("MIN AREA region excluded: %d :: %f\n", i, region[i].getArea());
        } else if (region[i].getPointCount() < MINCOORDS) {
            region[i].setAnalysisExcluded(true);
            printf("MIN COORD region excluded: %d :: %d\n", i, region[i].getPointCount());
        }

    for (i=0; i<this->RegionCount-1; i++)
        if (!region[i].isAnalysisExcluded()) {
            for (j=i+1; j<this->RegionCount; j++) {
                if (!region[j].isAnalysisExcluded()) { // if excluded from analysis, don't need to test again.
                    if (region[i].EncapsulationTest(region[j]))
                        region[i].setAnalysisExcluded(true);
                }
            } // j
        } // !analysisExcluded
} // filterRegions


Ray3F FireLine::getNormal(int index){
    return this->normals[index];
}


int FireLine::getNormalCount(){
    return this->nvCount;
}


void FireLine::DetermineCenterPoints() {
    for (int i=0; i<this->RegionCount-1; i++)
         region[i].CalculateCenterPoint();

} // DetermineCenterPoints
