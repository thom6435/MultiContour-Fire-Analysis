//
// Created by thomas on 1/25/21.
//

#include "thermalArea.h"
#include <math.h>
#include <stdio.h>


thermalArea::thermalArea() {
    this->initialized = false;
}

void thermalArea::setupAOI(int ImgWidth, int ImgHeight, float AerialFOV) {
    float Xe, Xe1, Xe2, Ye;

    this->AerialImgWidth = ImgWidth;
    this->AerialImgHeight = ImgHeight;

    float a = this->AerialImgWidth / 2.0;
    float b = this->AerialImgHeight / 2.0;

    this->AOI = (bool * ) malloc(sizeof(bool) * AerialImgWidth * AerialImgHeight);
    bool *AOIptr;

    for (int y = 0; y < AerialImgHeight; y++) {
        Ye = y - b;
        Xe = sqrt(pow(a, 2.0) * (1 - (pow(Ye, 2.0) / pow(b, 2.0))));
        Xe1 = a - Xe;
        Xe2 = a + Xe;

        AOIptr = this->AOI + y * this->AerialImgWidth;
        for (int x = 0; x < AerialImgWidth; x++) {
            if (x >= Xe1 && x <= Xe2)
                *AOIptr = true; // inside AOI
            else
                *AOIptr = false; // outside AOI
            AOIptr++;
        } // x
    } // y

    this->initialized = true;
}

bool thermalArea::isInAOI(int x, int y) {
    bool *AOIPtr;
    AOIPtr = this->AOI + y * this->AerialImgWidth + x;

    if (*AOIPtr)
        return true;

    return false;
    return false;
}

bool thermalArea::isInitialized() {
    return this->initialized;
}
