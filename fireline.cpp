#include "fireline.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


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
                                // point and is likely the product of a complex, convoluted fire line.
                                if (flNormal.innerContourIntercept(innerRegion)) {
                                    //puts("--------> ERROR Normal intersected originating fire line contour <-----------");
                                } else {
                                    if (flNormal.intercept(outerRegion)) {
                                        /* As the normals are saved from the temporary class used to generate them,
                                         * may want to check for multiple intersections from a common point.  Ie, this
                                         * could happen in the case of a very convoluted fire line.  In this case, we want
                                         * to keep the shortest, and discard the longer (ones).
                                         */

                                        this->normals[this->nvCount] = flNormal.getNormal();
                                        this->exclude[this->nvCount] = false;
                                        this->intersections[this->nvCount] = 0;
                                        this->nvCount++;
                                        //flNormal.print();
                                   }
                                }
                            } // closest inner point is not on image edge
                        } // OuterRegion point is not on image edge
                    } // K
                } // InnerRegion is encapsulated within outerRegion
                } /* else { // inner region ! analysis excluded
                    puts("INNER REGION ANALYSIS EXCLUDED-------------------------------------");
                } */
            } // j
        } // contour in second fire line is not encapsulated
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
            //printf("MIN AREA region excluded: %d :: %f\n", i, region[i].getArea());
        } else if (region[i].getPointCount() < MINCOORDS) {
            region[i].setAnalysisExcluded(true);
            //printf("MIN COORD region excluded: %d :: %d\n", i, region[i].getPointCount());
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


/**
 * @brief Normal::FilterNormals
 * @date 10-September-2020
 * Description: This method flags and removes normals that intersect other normals. Normals that
 *     intersect others are the product of a convoluted fireline and provide no meaningful
 *     direction or spread magnitude information.
 *     2) Flags normals that originate on the edge of the image.
 */
void FireLine::FilterNormals() {
    int i, j;

    // Test for Normal originating on/near edge of image;
    for (i=0; i<this->nvCount; i++) {
       if (this->normals[i].Pt.x > max_Image_X-2 ||
           this->normals[i].Pt.x < 1 ||
           this->normals[i].Pt.y > max_Image_Y-2 ||
           this->normals[i].Pt.y < 1)
        this->exclude[i] = true;
    } // i

    // Test for Normal to Normal Intersection;
    for (i=0; i<this->nvCount; i++) {
       for (j=0; j<this->nvCount; j++) {
           if (j!= i && !this->exclude[i]) { // avoid comparing normal to self;
           Ray3F A = getNormal(i);
           Ray3F B = getNormal(j);
           if (NormalToNormalIntersection(A, B))
               this->intersections[i]++;
           } // not self or excluded
       }// j
    } // i

    for (i=0; i<this->nvCount; i++)
        if (this->intersections[i] >=2)
            this->exclude[i] = true;

    // condense the array and remove the excluded normals
    int cleanCount = 0;
    for (i=0; i< this->nvCount; i++) {
        if (!this->exclude[i]) {
            this->normals[cleanCount].D.x = this->normals[i].D.x; // Not efficient or smart approach...
            this->normals[cleanCount].D.y = this->normals[i].D.y;
            this->normals[cleanCount].D.z = this->normals[i].D.z;
            this->normals[cleanCount].Pt.x = this->normals[i].Pt.x;
            this->normals[cleanCount].Pt.y = this->normals[i].Pt.y;
            this->normals[cleanCount].Pt.z = this->normals[i].Pt.z;
            this->exclude[cleanCount] = false;
            cleanCount++;
        } // excluded?
    } // i
} // FilterNormals


/**
 * @brief FireLine::NormalToNormalIntersection
 * @date 10-September-2020
 * this code is a variation of the normal::innerContourIntercept method
 * @param A Ray3F defining first normal vectorRay3F
 * @param B Ray3F defining second normal vector
 * @return true if Ray A intersects Ray B
 */
bool FireLine::NormalToNormalIntersection(Ray3F A, Ray3F B) {
    //int i=0;
    float r, s, t, a, b, D;
    float uu, uv, vv, wu, wv;
    Vector3F u, v, n, w0, w;
    point3F B1, B2, B3, Intercept;

    // Project the normal vector A outwards
    A.D.x *= 1000.0;
    A.D.y *= 1000.0;

    /*
     * Clockwise facet winding
     */

    // Assume right angle facet which extends into Z dimension 1 unit from the first point of the contour.
    B1.x = B.Pt.x;  B1.y = B.Pt.y;  B1.z = -1;
    B2.x = B.Pt.x + B.D.x * 1000.0; // project vector outwards to permit intersection detection.
    B2.y = B.Pt.y + B.D.y * 1000.0;
    B2.z = 0;
    B3.x = B1.x; B3.y = B1.y; B3.z = 1;

    Intercept.x = 0.0;
    Intercept.y = 0.0;
    Intercept.z = 0.0;

    u.x = B2.x - B1.x;
    u.y = B2.y - B1.y;
    u.z = B2.z - B1.z;
    v.x = B3.x - B1.x;
    v.y = B3.y - B1.y;
    v.z = B3.z - B1.z;

    n.x = u.y*v.z - u.z*v.y;
    n.y = u.z*v.x - u.x*v.z;
    n.z = u.x*v.y - u.y*v.x;

    // see if facet is degenerate -- if so, don't handle this case
    if (!(n.x==0 && n.y==0 && n.z==0)) {
        w0.x = A.Pt.x - B1.x;
        w0.y = A.Pt.y - B1.y;
        w0.z = A.Pt.z - B1.z;

        a = -(n.x*w0.x + n.y*w0.y + n.z*w0.z); // a = -( n dot w0)
        b = n.x*A.D.x + n.y*A.D.y +n.z*A.D.z;  // b = n dot dir

        if (fabs(b) > SMALL_NUM) {
                // determine intersect point of ray with facet
            r = a/b;
            if (r > 0.0) {
                Intercept.x = A.Pt.x + r*A.D.x;
                Intercept.y = A.Pt.y + r*A.D.y;
                Intercept.z = A.Pt.z + r*A.D.z;

                uu = u.x*u.x + u.y*u.y + u.z*u.z;  // u dot u
                uv = u.x*v.x + u.y*v.y + u.z*v.z;  // u dot v
                vv = v.x*v.x + v.y*v.y + v.z*v.z;  // v dot v

                w.x = Intercept.x - B1.x;
                w.y = Intercept.y - B1.y;
                w.z = Intercept.z - B1.z;

                wu = w.x*u.x + w.y*u.y + w.z*u.z;  // w dot u
                wv = w.x*v.x + w.y*v.y + w.z*v.z;  // w dot v

                D = uv*uv - uu*vv;
                s = (uv*wv - vv*wu) / D;
                t = (uv*wu - uu*wv) / D;

                if (!(s<0.0 || s>1.0) && !(t<0.0 || s+t>1.0) &&
                    !isnan(Intercept.x) && !isnan(Intercept.y) && !isnan(Intercept.z))
                    return true; // Intercept point with inner fireline contour
            }  // r > 0, ray points toward contour line facet
        } // normal not parallel to outercontour plane.
    } // degenerate facet test


    /*
     * Counter-Clockwise facet winding
     */

    // Assume right angle facet which extends into Z dimension 1 unit from the first point of the contour.
    B1.x = B.Pt.x;  B1.y = B.Pt.y;  B1.z = 1;
    B3.x = B.Pt.x + B.D.x * 1000.0; // project vector outwards to permit intersection detection.
    B3.y = B.Pt.y + B.D.y * 1000.0;
    B3.z = 0;
    B2.x = B1.x; B2.y = B1.y; B2.z = -1;

    Intercept.x = 0.0;
    Intercept.y = 0.0;
    Intercept.z = 0.0;

    u.x = B2.x - B1.x;
    u.y = B2.y - B1.y;
    u.z = B2.z - B1.z;

    v.x = B3.x - B1.x;
    v.y = B3.y - B1.y;
    v.z = B3.z - B1.z;

    n.x = u.y*v.z - u.z*v.y;
    n.y = u.z*v.x - u.x*v.z;
    n.z = u.x*v.y - u.y*v.x;

    // see if facet is degenerate -- if so, don't handle this case
    if (!(n.x==0 && n.y==0 && n.z==0)) {
        w0.x = A.Pt.x - B1.x;
        w0.y = A.Pt.y - B1.y;
        w0.z = A.Pt.z - B1.z;

        a = -(n.x*w0.x + n.y*w0.y + n.z*w0.z); // a = -( n dot w0)
        b = n.x*A.D.x + n.y*A.D.y +n.z*A.D.z;  // b = n dot dir

        if (fabs(b) > SMALL_NUM) {
                // determine intersect point of ray with facet
            r = a/b;
            if (r > 0.0) {
                Intercept.x = A.Pt.x + r*A.D.x;
                Intercept.y = A.Pt.y + r*A.D.y;
                Intercept.z = A.Pt.z + r*A.D.z;

                uu = u.x*u.x + u.y*u.y + u.z*u.z;  // u dot u
                uv = u.x*v.x + u.y*v.y + u.z*v.z;  // u dot v
                vv = v.x*v.x + v.y*v.y + v.z*v.z;  // v dot v

                w.x = Intercept.x - B1.x;
                w.y = Intercept.y - B1.y;
                w.z = Intercept.z - B1.z;

                wu = w.x*u.x + w.y*u.y + w.z*u.z;  // w dot u
                wv = w.x*v.x + w.y*v.y + w.z*v.z;  // w dot v

                D = uv*uv - uu*vv;
                s = (uv*wv - vv*wu) / D;
                t = (uv*wu - uu*wv) / D;

                if (!(s<0.0 || s>1.0) && !(t<0.0 || s+t>1.0) &&
                    !isnan(Intercept.x) && !isnan(Intercept.y) && !isnan(Intercept.z))
                    return true; // Intercept point with inner fireline contour
            }  // r > 0, ray points toward contour line facet
        } // normal not parallel to outercontour plane.
    } // degenerate facet test


    return false; // No Intercept point found!
}
