#ifndef NORMAL_H
#define NORMAL_H

#include "contourregion.h"
#include "opencv2/core.hpp"


#define SMALL_NUM            0.00000001  // anything that avoids division overflow
#define FIT_POINT_MAX        7

typedef struct floatPointStruct {
    float x;
    float y;
} point2F;

typedef struct float3DPointStruct {
    float x;
    float y;
    float z;
} point3F;

typedef struct Float2DVectorStruct {
    float x;
    float y;
} Vector2F;

typedef struct Float3DVectorStruct {
    float x;
    float y;
    float z;
} Vector3F;

typedef struct Ray3FStruct {
    point3F Pt;
    Vector3F D;
} Ray3F;

class Normal
{
    private:
        point2F ptVector[FIT_POINT_MAX]; // points used to determine line on the contour
        point2F A, B; // points on the fitted line
        int pointCount;
        int imageMaxX, imageMaxY; // Image bounds
        Vector2F FittedLineVector;
        Vector2F FireProgressionVector;
        point2F  normalInterceptPt;
        Ray3F N;          // N is the normal vector originating from the midpoint of the line segment
                          // described by ptVector.

    public:
        Normal();
        Normal(point2D fittingPoints[FIT_POINT_MAX], int imageMaxX, int imageMaxY);
        Normal(point2F fittingPoints[FIT_POINT_MAX], int imageMaxX, int imageMaxY);
        Normal(cv::Point A, cv::Point B, int imageMaxX, int imageMaxY);
        void print();
        Ray3F getNormal();
        point2F CalculateMidPoint();
        point2F getInterceptPt();
        void FitBestLine();
        Ray3F CalculateNormal(point2D directionHint);
        void CalculateNormals(Vector3F *V1, Vector3F *V2);
        bool FindNormalRegionIntersect(ContourRegion OuterFireContour);
        bool intercept(ContourRegion OuterFireContour);
        bool intercept(point3F p, Vector3F nv, std::vector<cv::Point> fireline, point3F *InterceptPt);
        bool interceptRadiometerAOI(point3F p, Vector3F nv, std::vector<cv::Point> fireline, point3F *InterceptPt);
        bool innerContourIntercept(ContourRegion OuterFireContour);
};

#endif // NORMAL_H
