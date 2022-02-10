#ifndef FIRELINE_H
#define FIRELINE_H

#define MAXCOORDS 1024   // maximum number of coords that can be contained in a fire contour
#define MINCOORDS  200   // threshold for critical number of points in a fire contour
#define MINAREA   2000.0 // Critical size threshold for area in a contour region segment of a fireline

#include "contourregion.h"
#include "normal.h"

#define MAXNORMALS 8092 // maxnumber of normals generated per fireline contour

class FireLine
{
private:
    int RegionCount; // number of contour regions in the fireline definition
    ContourRegion *region;
    int max_Image_Y, max_Image_X;  // Image height and width, used to avert processing fireline contours on image edge.
    Ray3F normals[MAXNORMALS];
    bool exclude[MAXNORMALS];
    int intersections[MAXNORMALS];
    int nvCount = 0;  // # of normal vectors connected to this fireline from nested fireline

public:
    FireLine(int regionCount, int maxX, int maxY);
    FireLine();
    void setRegionCount(int regionCount);
    int getRegionCount();
    ContourRegion initializeContourRegion(int index, int pointCount);
    ContourRegion getContourRegion(int index);
    void addContourRegions(ContourRegion cr[]);
    point2D getPoint2D(int ContourID, int PointID);
    void addNextPoint(int ContourID, float x, float y);
    void progressionAnalysis(FireLine nextFireLine);
    void filterRegions();
    void DetermineCenterPoints();
    int MinContourPointThreshold();
    Ray3F getNormal(int index);
    int getNormalCount();
    void setImageY(int y);
    void setImageX(int x);
    int getImageY();
    int getImageX();
    void print();
    void FilterNormals();
    bool NormalToNormalIntersection(Ray3F A, Ray3F B);
};

#endif // FIRELINE_H
