#ifndef CONTOURREGION_H
#define CONTOURREGION_H

typedef struct intPointStruct {
    int x;
    int y;
} point2D;


class ContourRegion
{
    private:
        float area;
        int pointCount;
        point2D center;
        point2D *fireline;
        int index;
        bool encapsulated = false;
        bool analysisExcluded = false;

    public:
        ContourRegion();
        ContourRegion(int pointCount);
        void copy(point2D *p, float area);
        void setArea(float area);
        void setEncapsulationFlag(bool flag);
        int getPointCount();
        float getArea();
        void addPoint(int x, int y);
        point2D getPoint(int index);
        point2D *getAllPoints();
        point2D CalculateCenterPoint();
        point2D getCenter();
        void setCenter(point2D cntr);
        point2D getPreviousPoint(int index);
        point2D getNextPoint(int index);
        int getNearestPointIndex(point2D p);
        void PrintCoordinates();
        void PrintSummary(char *message);
        bool EncapsulationTest(ContourRegion cr);
        bool isEncapsulated();
        bool isAnalysisExcluded();
        void setAnalysisExcluded(bool b);
};

#endif // CONTOURREGION_H
