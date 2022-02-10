//
// Created by thomas on 1/25/21.
//

#include "radiometer.h"
#include <math.h>
#include <stdio_ext.h>

radiometer::radiometer() {
    this->diameterMeters = 0;
    this->initialized = false;
}

radiometer::radiometer(float diameter) {
    this->diameterMeters = diameter;
    this->initialized = false;
}

void radiometer::setupAOI(int ImgWidth, int ImgHeight, float AerialFOV) {
    float Xe, Xe1, Xe2, Ye;
    this->AerialImgWidth = ImgWidth;
    this->AerialImgHeight = ImgHeight;

    this->diameterPixels = this->diameterMeters * (ImgWidth/AerialFOV);

    printf("radiometer::setupAOI  diameter in meters: %f\n", this->diameterMeters);
    printf("radiometer::setupAOI  diameter in pixels: %d\n", this->diameterPixels);

    float a = AerialImgWidth / 2.0;
    float b = AerialImgHeight / 2.0;

    this->AOI = (bool * ) malloc(sizeof(bool) * AerialImgWidth * AerialImgHeight);
    bool *AOIptr;
    for (int y = 0; y < AerialImgHeight; y++) {
        Ye = y - b;
        Xe = sqrt(pow(this->diameterPixels/2.0,2.0)*(1-(pow(Ye,2.0)/pow(this->diameterPixels/2.0,2.0))));
        Xe1 = a - Xe;
        Xe2 = a + Xe;

        AOIptr = this->AOI + y * this->AerialImgWidth;
        for (int x=0; x < AerialImgWidth; x++) {
            if (x >= Xe1 && x <= Xe2)
                *AOIptr = true; // inside AOI
            else
                *AOIptr = false; // outside AOI
            AOIptr++;
        } // x
    } // y

    this->initialized = true;
} // setupAOI;


bool radiometer::isInAOI(int x, int y) {
    bool *AOIPtr;
    AOIPtr = this->AOI + y * this->AerialImgWidth + x;

    if (*AOIPtr)
        return true;

    return false;
}

void radiometer::setPixelDiameter(float diameter) {
    this->diameterMeters = diameter;
}

bool radiometer::isInitialized() {
   return this->initialized;
}

double radiometer::getArea() {
    double area;
    area = pow((diameterMeters / 2.0), 2.0) * 3.141592653;

    return area;
}
