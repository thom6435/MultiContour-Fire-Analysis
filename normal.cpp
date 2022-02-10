/*
 *       Class: Normal
 *    Modified: 10-September-2020
 * Description:
 *
 */

#include "normal.h"
#include "contourregion.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


Normal::Normal()
{
}

Normal::Normal(point2D fittingPoints[FIT_POINT_MAX], int MaxX, int MaxY){
   this->imageMaxX = MaxX;
   this->imageMaxY = MaxY;

    for (int i=0; i<FIT_POINT_MAX; i++) {
        ptVector[i].x = (float) fittingPoints[i].x;
        ptVector[i].y = (float) fittingPoints[i].y;
    }
    this->pointCount = FIT_POINT_MAX;

    FitBestLine();
}


Normal::Normal(point2F fittingPoints[FIT_POINT_MAX], int MaxX, int MaxY){
    this->imageMaxX = MaxX;
    this->imageMaxY = MaxY;

    for (int i=0; i<FIT_POINT_MAX; i++) {
        ptVector[i].x = fittingPoints[i].x;
        ptVector[i].y = fittingPoints[i].y;
    }
    this->pointCount = FIT_POINT_MAX;

    FitBestLine();
}

Normal::Normal(cv::Point A, cv::Point B, int MaxX, int MaxY) {
    this->imageMaxX = MaxX;
    this->imageMaxY = MaxY;
    this->A.x = A.x;
    this->A.y = A.y;
    this->B.x = B.x;
    this->B.y = B.y;

    this->pointCount = 2;

    this->FittedLineVector.x = this->B.x - this->A.x;
    this->FittedLineVector.y = this->B.y - this->A.y;

    // Calculate the midpoint and store it as the origin of the normal.
    N.Pt.x = (A.x + B.x)/2.0;
    N.Pt.y = (A.y + B.y)/2.0;
    N.Pt.z = 0;
}

void Normal::print(){
    printf("NORMAL:: (%8.4f,%8.4f) ---> <%8.4f,%8.4f>\n", this->N.Pt.x, this->N.Pt.y, this->N.D.x, this->N.D.y);
}

Ray3F Normal::getNormal(){
    return this->N;
}


/*
 *   method: FitBestLine
 * modified: 18-August-2020
 *  purpose: Uses a simple least-squares method to fit a line to as many as five contiguous points
 *           on the inner fireline contour.  The more accurate the line fit, the more meaningful
 *           the fire line progression direction and magnitude estimate of the normal vector.
 */
void Normal::FitBestLine() {
    double meanX=0.0, meanY=0.0;
    double xy_dif=0.0, x_dif=0.0;
    double slope, yIntercept;
    //point2F A, B;
    int i;

    for (i=0; i<FIT_POINT_MAX; i++) {
        meanX += this->ptVector[i].x;
        meanY += this->ptVector[i].y;
    } // i

    meanX = meanX / FIT_POINT_MAX;
    meanY = meanY / FIT_POINT_MAX;

    for (i=0; i<FIT_POINT_MAX; i++) {
        xy_dif += (this->ptVector[i].x - meanX) * (this->ptVector[i].y - meanY);
        x_dif += pow(this->ptVector[i].x - meanX, 2);
    } // i

    slope = xy_dif / x_dif;

    if (isnan(slope)) {
        slope = 0;
        A.x = ptVector[0].x;
        A.y = ptVector[0].y;
        B.x = ptVector[FIT_POINT_MAX-1].x;
        B.y = ptVector[FIT_POINT_MAX-1].y;
    } else {
        if ((xy_dif < 0.000001 && xy_dif > -0.000001) || (slope < 0.000001 && slope > -0.000001)) {
            slope = 0;
            A.x = ptVector[0].x;
            A.y = ptVector[0].y;
            B.x = ptVector[FIT_POINT_MAX-1].x;
            B.y = ptVector[FIT_POINT_MAX-1].y;
        } else {
            yIntercept = meanY - slope * meanX;

            A.x = ptVector[0].x;
            A.y = slope*A.x + yIntercept;
            B.x = ptVector[6].x;
            B.y = slope*B.x + yIntercept;
        }
    }

    this->FittedLineVector.x = B.x - A.x;
    this->FittedLineVector.y = B.y - A.y;

    // Calculate the midpoint and store it as the origin of the normal.
    N.Pt.x = (A.x + B.x)/2.0;
    N.Pt.y = (A.y + B.y)/2.0;
    N.Pt.z = 0;
/*
    printf("FitBestLine:  A: %f, %f\n", A.x, A.y);
    printf("              B: %f, %f\n", B.x, B.y);
    printf("              N: %f, %f\n",N.Pt.x, N.Pt.y);
    printf("     FittedVctr: %f, %f\n", this->FittedLineVector.x, this->FittedLineVector.y);
*/
} //FitBestLine


/*
 *   method: FindNormalRegionIntersect
 * modified: 21-August-2020
 *  purpose:
 *  returns: true  If normal was found to intersect the fire contour line
 *           false If the normal was not found to intersect the contour (likely
 *                 a logic/math error to be flagged and examined)
 */
bool Normal::FindNormalRegionIntersect(ContourRegion OuterFireContour) {
    int i=0, j;
    Ray3F F;
    float t, s, div;
    point2F Intercept;

    int pointCount = OuterFireContour.getPointCount();

    while (i < pointCount) {
        F.Pt.x = OuterFireContour.getPoint(i).x;
        F.Pt.y = OuterFireContour.getPoint(i).y;

        j = i + 1;
        if (j >= pointCount)
            j = 0;
        F.D.x = OuterFireContour.getPoint(j).x - F.Pt.x;
        F.D.y = OuterFireContour.getPoint(j).y - F.Pt.y;

        div = N.D.y * F.D.x - N.D.x * F.D.y;

        //printf("%d ----------------------------------------------\n", i);
        //this->print();
        //printf("F:: (%f,%f)---->%f,%f\n", F.Pt.x, F.Pt.y, F.D.x, F.D.y);
        //printf("Div = %f\n", div);

        if (div != 0) { // Rays intersect, determine where.
            t = (F.D.x * (F.Pt.y - N.Pt.y)) + (F.D.y *(N.Pt.x - F.Pt.x)) / div;
            s = (N.D.x * (F.Pt.y - N.Pt.y)) + (N.D.y *(N.Pt.x - F.Pt.x)) / div;

            /*
             * S is the distance along the fireContour ray F to the intersection point.
             * T is the distance along the Normal N to the intersection point.
             * If S is less than or equal to the length of the line segment described by the
             * F Ray, then the Normal intersects F.
             */
            float magnitudeF = sqrt(pow(F.D.x,2) + pow(F.D.y,2));
            if (s <= magnitudeF) {
                Intercept.x = N.Pt.x + t*N.D.x;
                Intercept.y = N.Pt.y + t*N.D.y;


                // See if the intersection point is inside the bounds of the image
                if (Intercept.x >=0  && Intercept.x < this->imageMaxX &&
                    Intercept.y >=0  && Intercept.y < this->imageMaxY) {

                    // update the Normal Vector to reflect the intersection point
                    N.D.x = N.D.x * t;
                    N.D.y = N.D.y * t;
                    N.D.z = 0;

                    printf("NormalVector Intercept @ %f,%f  NV=(%f,%f)---> {%f,%f}\n",
                            Intercept.x,Intercept.y,N.Pt.x,N.Pt.y,N.D.x,N.D.y);

                    return true;
                } // intercept in bounds
            }
        }
        i++;
    } // while

    return false;
} // FindNormalRegionIntersect


/**
 @function intercept
 @date 29-December 2020
 @note This version of the intercept method uses the point (p) and normal (nv)
       instead of the N class variable for determining the intercept point.
 @param p origination point of vector for intercept test
 @param nv perpendicular to line fitted to leading edge of the fireline
 @param fireline vector of points that describe the active burning region
 @param InterceptPt point to discovered intercept, if found
 @returns bool true if intercept found, false otherwise
 */
bool Normal::intercept(point3F p, Vector3F nv, std::vector<cv::Point> fireline,
                       point3F *InterceptPt)
{
    int i=0, j;
    float r, s, t, a, b, D;
    float uu, uv, vv, wu, wv;
    Vector3F u, v, n, w0, w;
    point3F A, B, C, I;
    bool interceptFound = false;
    float bestMagnitude = -1000000000; //absurdly small initial value...

    while (i < fireline.size()) {
        // Assume right angle facet which extends into Z dimension 1 unit from the first point of the contour.
        I.x = 0.0;
        I.y = 0.0;
        I.z = 0.0;

        A.x = fireline[i].x;
        A.y = fireline[i].y;
        A.z = 0.0;

        j = i + 1;
        if (j >= fireline.size()) // if at end of list, connect back to start of list
            j = 0;
        B.x = fireline[j].x;
        B.y = fireline[j].y;
        B.z = 0.0;

        C.x = A.x;
        C.y = A.y;
        C.z = 1.0;

        u.x = B.x - A.x;
        u.y = B.y - A.y;
        u.z = B.z - A.z;

        v.x = C.x - A.x;
        v.y = C.y - A.y;
        v.z = C.z - A.z;

        n.x = u.y*v.z - u.z*v.y;
        n.y = u.z*v.x - u.x*v.z;
        n.z = u.x*v.y - u.y*v.x;

        // test that facet is not degenerate -- if so, don't handle this case
        if (!(n.x==0 && n.y==0 && n.z==0)) {
            w0.x = p.x - A.x;
            w0.y = p.y - A.y;
            w0.z = p.z - A.z;

            a = -(n.x*w0.x + n.y*w0.y + n.z*w0.z); // a = -( n dot w0)
            b = n.x*nv.x + n.y*nv.y +n.z*nv.z;  // b = n dot dir

            if (fabs(b) > SMALL_NUM) {
                // determine intersect point of ray with facet
                r = a/b;
                if (r > 0.0) { // ray points toward contour line facet.
                    I.x = p.x + r*nv.x;
                    I.y = p.y + r*nv.y;
                    I.z = p.z + r*nv.z;

                    uu = u.x*u.x + u.y*u.y + u.z*u.z;  // u dot u
                    uv = u.x*v.x + u.y*v.y + u.z*v.z;  // u dot v
                    vv = v.x*v.x + v.y*v.y + v.z*v.z;  // v dot v

                    w.x = I.x - A.x;
                    w.y = I.y - A.y;
                    w.z = I.z - A.z;

                    wu = w.x*u.x + w.y*u.y + w.z*u.z;  // w dot u
                    wv = w.x*v.x + w.y*v.y + w.z*v.z;  // w dot v

                    D = uv*uv - uu*vv;
                    s = (uv*wv - vv*wu) / D;
                    t = (uv*wu - uu*wv) / D;

                    if (!(s<0.0 || s>1.0) && !(t<0.0 || s+t>1.0) &&
                        !isnan(I.x) && !isnan(I.y) && !isnan(I.z)) {

                        if (I.x > .5 && I.x < this->imageMaxX-1.5 &&  // exclude intercepts with the image edge
                            I.y > .5 && I.y < this->imageMaxY-1.5) {
                            interceptFound = true;
                            float magnitude = sqrt(pow(I.x-p.x,2) +
                                                      pow(I.y-p.y,2));
                            if (magnitude > bestMagnitude) {   // loop through all points to find closest intercept
                                bestMagnitude = magnitude;
                                InterceptPt->x = I.x;
                                InterceptPt->y = I.y;
                                InterceptPt->z = I.z;
                            } // magnitude comparison
                        } // intercept not on image edge
                    } // intercept coords are real numbers
                }   // r > 0, ray points toward contour line facet
            }  // normal not parallel to outer contour plane.
        }

        i++;
    } // i

    return interceptFound;
} // intercept


/**
 @function interceptRadiometerAOI
 @date 29-January-2021
 @note This version of the intercept method uses the point (p) and normal (nv)
       instead of the N class variable for determining the intercept point. This
       version also permits intercepts with the edge of the image, as the radiometer
       AOI can extend to the image edge.
 @param p origination point of vector for intercept test
 @param nv perpendicular to line fitted to leading edge of the fireline
 @param fireline vector of points that describe the active burning region
 @param InterceptPt point to discovered intercept, if found
 @returns bool true if intercept found, false otherwise
 */
bool Normal::interceptRadiometerAOI(point3F p, Vector3F nv, std::vector<cv::Point> fireline,
                       point3F *InterceptPt)
{
    int i=0, j;
    float r, s, t, a, b, D;
    float uu, uv, vv, wu, wv;
    Vector3F u, v, n, w0, w;
    point3F A, B, C, I;
    bool interceptFound = false;
    float bestMagnitude = -1000000000; //absurdly small initial value...

    while (i < fireline.size()) {
        // Assume right angle facet which extends into Z dimension 1 unit from the first point of the contour.
        I.x = 0.0;
        I.y = 0.0;
        I.z = 0.0;

        A.x = fireline[i].x;
        A.y = fireline[i].y;
        A.z = 0.0;

        j = i + 1;
        if (j >= fireline.size()) // if at end of list, connect back to start of list
            j = 0;
        B.x = fireline[j].x;
        B.y = fireline[j].y;
        B.z = 0.0;

        C.x = A.x;
        C.y = A.y;
        C.z = 1.0;

        u.x = B.x - A.x;
        u.y = B.y - A.y;
        u.z = B.z - A.z;

        v.x = C.x - A.x;
        v.y = C.y - A.y;
        v.z = C.z - A.z;

        n.x = u.y*v.z - u.z*v.y;
        n.y = u.z*v.x - u.x*v.z;
        n.z = u.x*v.y - u.y*v.x;

        // test that facet is not degenerate -- if so, don't handle this case
        if (!(n.x==0 && n.y==0 && n.z==0)) {
            w0.x = p.x - A.x;
            w0.y = p.y - A.y;
            w0.z = p.z - A.z;

            a = -(n.x*w0.x + n.y*w0.y + n.z*w0.z); // a = -( n dot w0)
            b = n.x*nv.x + n.y*nv.y +n.z*nv.z;  // b = n dot dir

            if (fabs(b) > SMALL_NUM) {
                // determine intersect point of ray with facet
                r = a/b;
                if (r > 0.0) { // ray points toward contour line facet.
                    I.x = p.x + r*nv.x;
                    I.y = p.y + r*nv.y;
                    I.z = p.z + r*nv.z;

                    uu = u.x*u.x + u.y*u.y + u.z*u.z;  // u dot u
                    uv = u.x*v.x + u.y*v.y + u.z*v.z;  // u dot v
                    vv = v.x*v.x + v.y*v.y + v.z*v.z;  // v dot v

                    w.x = I.x - A.x;
                    w.y = I.y - A.y;
                    w.z = I.z - A.z;

                    wu = w.x*u.x + w.y*u.y + w.z*u.z;  // w dot u
                    wv = w.x*v.x + w.y*v.y + w.z*v.z;  // w dot v

                    D = uv*uv - uu*vv;
                    s = (uv*wv - vv*wu) / D;
                    t = (uv*wu - uu*wv) / D;

                    if (!(s<0.0 || s>1.0) && !(t<0.0 || s+t>1.0) &&
                        !isnan(I.x) && !isnan(I.y) && !isnan(I.z)) {
                        interceptFound = true;
                        float magnitude = sqrt(pow(I.x-p.x,2) +
                                                  pow(I.y-p.y,2));
                        if (magnitude > bestMagnitude) {   // loop through all points to find closest intercept
                            bestMagnitude = magnitude;
                            InterceptPt->x = I.x;
                            InterceptPt->y = I.y;
                            InterceptPt->z = I.z;
                        } // magnitude comparison
                    } // intercept coords are real numbers
                }   // r > 0, ray points toward contour line facet
            }  // normal not parallel to outer contour plane.
        }

        i++;
    } // i

    return interceptFound;
} // interceptRadiometerAOI


 /*
 *   method: intercept
 * modified: 26-August-2020
 *  purpose: Revised normal / fire line contour intercept method based on ray-tracing approach.
 */
bool Normal::intercept(ContourRegion OuterFireContour) {
    int i=0, j;
    float r, s, t, a, b, D;
    float uu, uv, vv, wu, wv;
    Vector3F u, v, n, w0, w;
    point3F Intercept, A, B, C;

    while (i < OuterFireContour.getPointCount()) {
        // Assume right angle facet which extends into Z dimension 1 unit from the first point of the contour.
        Intercept.x = 0.0;
        Intercept.y = 0.0;
        Intercept.z = 0.0;

        A.x = OuterFireContour.getPoint(i).x;
        A.y = OuterFireContour.getPoint(i).y;
        A.z = 0.0;

        j = i + 1;
        if (j >= OuterFireContour.getPointCount()) // if at end of list, connect back to start of list
            j = 0;
        B.x = OuterFireContour.getPoint(j).x;
        B.y = OuterFireContour.getPoint(j).y;
        B.z = 0.0;

        C.x = A.x;
        C.y = A.y;
        C.z = 1.0;

        u.x = B.x - A.x;
        u.y = B.y - A.y;
        u.z = B.z - A.z;

        v.x = C.x - A.x;
        v.y = C.y - A.y;
        v.z = C.z - A.z;

        n.x = u.y*v.z - u.z*v.y;
        n.y = u.z*v.x - u.x*v.z;
        n.z = u.x*v.y - u.y*v.x;

        // test that facet is not degenerate -- if so, don't handle this case
        if (!(n.x ==0 && n.y ==0 && n.z ==0)) {
            w0.x = N.Pt.x - A.x;
            w0.y = N.Pt.y - A.y;
            w0.z = N.Pt.z - A.z;

            a = -(n.x*w0.x + n.y*w0.y + n.z*w0.z); // a = -( n dot w0)
            b = n.x*N.D.x + n.y*N.D.y +n.z*N.D.z;  // b = n dot dir

            if (fabs(b) > SMALL_NUM) {
                // determine intersect point of ray with facet
                r = a/b;
                if (r > 0.0) { // ray points toward contour line facet.
                     Intercept.x = N.Pt.x + r*N.D.x;
                     Intercept.y = N.Pt.y + r*N.D.y;
                     Intercept.z = N.Pt.z + r*N.D.z;

                     uu = u.x*u.x + u.y*u.y + u.z*u.z;  // u dot u
                     uv = u.x*v.x + u.y*v.y + u.z*v.z;  // u dot v
                     vv = v.x*v.x + v.y*v.y + v.z*v.z;  // v dot v

                     w.x = Intercept.x - A.x;
                     w.y = Intercept.y - A.y;
                     w.z = Intercept.z - A.z;

                     wu = w.x*u.x + w.y*u.y + w.z*u.z;  // w dot u
                     wv = w.x*v.x + w.y*v.y + w.z*v.z;  // w dot v

                     D = uv*uv - uu*vv;
                     s = (uv*wv - vv*wu) / D;
                     t = (uv*wu - uu*wv) / D;

                      if (!(s<0.0 || s>1.0) && !(t<0.0 || s+t>1.0) &&
                         !isnan(Intercept.x) && !isnan(Intercept.y) && !isnan(Intercept.z)) {

                         this->normalInterceptPt.x = Intercept.x;
                         this->normalInterceptPt.y = Intercept.y;

                         // update the Normal Vector to reflect the intersection point
                         this->N.D.x = Intercept.x - this->N.Pt.x;
                         this->N.D.y = Intercept.y - this->N.Pt.y;

                         return true; // Intercept point found!
                     }
                }   // r > 0, ray points toward contour line facet
            }  // normal not parallel to outercontour plane.
        }

        i++;
    } // i

    //puts("No Intercept");

    return false; // No Intercept point found!
} // intercept




/*
 *   method: innerContourIntercept
 * modified: 27-August-2020
 *  purpose: Tests to see if the normal intercepts the inner, or originating, fire contour line.
 *           if it does then the inner contour is regarded as being convoluted to the degree
 *           that no forward progression to the outer fire line can be accurately determined at
 *           this point.
 *     Note: This method is identical to the intercept function except that it moves the normal
 *           originating point forward such that it is just off the inner contour.  This is to
 *           false positive intercept tests which will occur if the true originating point is used.
 */
bool Normal::innerContourIntercept(ContourRegion innerContour){
    int i=0, j;
    float r, s, t, a, b, D;
    float uu, uv, vv, wu, wv;
    Vector3F u, v, n, w0, w;
    point3F Intercept, A, B, C;
    Ray3F N_prime;

    N_prime.D.x = N.D.x; N_prime.D.y = N.D.y; N_prime.D.z = N.D.z;
    N_prime.Pt.x = N.Pt.x + 4.0 * N.D.x; // 1.0 == 1 pixel
    N_prime.Pt.y = N.Pt.y + 4.0 * N.D.y;
    N_prime.Pt.z = N.Pt.z;

    int pointCount = innerContour.getPointCount();

    while (i < pointCount) {
        // Assume right angle facet which extends into Z dimension 1 unit from the first point of the contour.
        Intercept.x = 0.0;
        Intercept.y = 0.0;
        Intercept.z = 0.0;

        A.x = innerContour.getPoint(i).x;
        A.y = innerContour.getPoint(i).y;
        A.z = 0.0;

        j = i + 1;
        if (j>=pointCount)
            j=0;
        B.x = innerContour.getPoint(j).x;
        B.y = innerContour.getPoint(j).y;
        B.z = 0.0;

        C.x = innerContour.getPoint(i).x;
        C.y = innerContour.getPoint(i).y;
        C.z = 1.0;

        u.x = B.x - A.x;
        u.y = B.y - A.y;
        u.z = B.z - A.z;

        v.x = C.x - A.x;
        v.y = C.y - A.y;
        v.z = C.z - A.z;

        n.x = u.y*v.z - u.z*v.y;
        n.y = u.z*v.x - u.x*v.z;
        n.z = u.x*v.y - u.y*v.x;

        // see if facet is degenerate -- if so, don't handle this case
        if (!(n.x==0 && n.y==0 && n.z==0)) {
            w0.x = N_prime.Pt.x - A.x;
            w0.y = N_prime.Pt.y - A.y;
            w0.z = N_prime.Pt.z - A.z;

            a = -(n.x*w0.x + n.y*w0.y + n.z*w0.z); // a = -( n dot w0)
            b = n.x*N_prime.D.x + n.y*N_prime.D.y +n.z*N_prime.D.z;  // b = n dot dir

            if (fabs(b) > SMALL_NUM) {
                // determine intersect point of ray with facet
                r = a/b;
                if (r > 0.0) {
                     Intercept.x = N_prime.Pt.x + r*N_prime.D.x;
                     Intercept.y = N_prime.Pt.y + r*N_prime.D.y;
                     Intercept.z = N_prime.Pt.z + r*N_prime.D.z;

                     uu = u.x*u.x + u.y*u.y + u.z*u.z;  // u dot u
                     uv = u.x*v.x + u.y*v.y + u.z*v.z;  // u dot v
                     vv = v.x*v.x + v.y*v.y + v.z*v.z;  // v dot v

                     w.x = Intercept.x - A.x;
                     w.y = Intercept.y - A.y;
                     w.z = Intercept.z - A.z;

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
        i++;
    } // i

    return false; // No Intercept point found!
}


point2F Normal::getInterceptPt() {
    return this->normalInterceptPt;
}


/*
 *   method: CalculateNormal
 * modified: 20-August-2020
 *  purpose: Calculates the normal from the midpoint (N.pt) towards the next FireLine. DirectionHint
 *           is used to determine the normal direction towards the progression of the fireline.
 */
Ray3F Normal::CalculateNormal(point2D directionHint) {
    Vector2F V1, V2; // 2 potential normal vectors

    V1.x = this->FittedLineVector.y;
    V1.y = -1.0 * this->FittedLineVector.x;
    float magnitudeV1 = sqrt(pow(V1.x,2) + pow(V1.y,2));
    // Normalize V1
    V1.x = V1.x / magnitudeV1;
    V1.y = V1.y / magnitudeV1;

    V2.x = -1.0 * this->FittedLineVector.y;
    V2.y = this->FittedLineVector.x;
    float magnitudeV2 = sqrt(pow(V2.x,2) + pow(V2.y,2));
    // Normalize V2
    V2.x = V2.x / magnitudeV2;
    V2.y = V2.y / magnitudeV2;

    // given the normal vectors V1 & V2 pointing in opposite directions, perpendicular to
    // the fittedLine vector.  Determine which one is pointing toward the directionHint
    // point.

    float distanceV1, distanceV2;
    point2F projectedV1, projectedV2;

    projectedV1.x = this->N.Pt.x + V1.x*5;
    projectedV1.y = this->N.Pt.y + V1.y*5;
    distanceV1 = sqrt(pow(((float) directionHint.x) - projectedV1.x, 2) +
                      pow(((float) directionHint.y) - projectedV1.y, 2));

    projectedV2.x = this->N.Pt.x + V2.x*5;
    projectedV2.y = this->N.Pt.y + V2.y*5;
    distanceV2 = sqrt(pow(((float) directionHint.x) - projectedV2.x, 2) +
                      pow(((float) directionHint.y) - projectedV2.y, 2));

    // the point that is closest to DirectionHint is a product of the
    // normal vector that represents the direction of fire progression.

    if (distanceV1 < distanceV2) {
        this->N.D.x = V1.x;
        this->N.D.y = V1.y;        
        //puts("V1 is Normal");
    } else {
        this->N.D.x = V2.x;
        this->N.D.y = V2.y;
        //puts("V2 is Normal");
    }

    this->N.D.z = 0.0;

/*
    printf("CalcNormal: hint = %d, %d\n", directionHint.x, directionHint.y);
    printf("              V1 = %f, %f\n", V1.x, V1.y);
    printf("            V1.p = %f, %f\n", projectedV1.x,  projectedV1.y);
    printf("              V2 = %f, %f\n", V2.x, V2.y);
    printf("            V2.p = %f, %f\n", projectedV2.x,  projectedV2.y);
    this->print();
*/

    return N;
}// CalculateNormals

/**
 * @function CalculateNormals
 * @date 29-December-2020
 * @note
 * @param Ray3F V1
 *        Ray3F V2
 *        V1 and V2 are the normal vectors (2D) to the fitted line AB
 * @return nothing
 */
void Normal::CalculateNormals(Vector3F *V1, Vector3F *V2) {
    V1->x = this->FittedLineVector.y;
    V1->y = -1.0 * this->FittedLineVector.x;
    float magnitudeV1 = sqrt(pow(V1->x,2) + pow(V1->y,2));
    // Normalize V1
    V1->x = V1->x / magnitudeV1;
    V1->y = V1->y / magnitudeV1;
    V1->z = 0.0;

    V2->x = -1.0 * this->FittedLineVector.y;
    V2->y = this->FittedLineVector.x;
    float magnitudeV2 = sqrt(pow(V2->x,2) + pow(V2->y,2));
    // Normalize V2
    V2->x = V2->x / magnitudeV2;
    V2->y = V2->y / magnitudeV2;
    V2->z = 0.0;
}// CalculateNormal

