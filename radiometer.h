//
// Created by thomas on 1/25/21.
//

#ifndef MULTICONTOUR_RADIOMETER_H
#define MULTICONTOUR_RADIOMETER_H

class radiometer {
private:
    float diameterMeters;
    int   diameterPixels;
    bool *AOI;
    bool initialized;
    int AerialImgWidth, AerialImgHeight;

public:
    radiometer();
    radiometer(float diameterMeters);
    void setPixelDiameter(float DiameterMeters);
    double getArea();
    void setupAOI(int AerialImgWidth, int AerialImgHeight, float AerialFOV);
    bool isInAOI(int x, int y);
    bool isInitialized();
};

#endif //MULTICONTOUR_RADIOMETER_H
