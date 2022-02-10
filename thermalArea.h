//
// Created by thomas on 1/25/21.
//

#ifndef MULTICONTOUR_THERMALAREA_H
#define MULTICONTOUR_THERMALAREA_H

class thermalArea {
private:
    bool *AOI;
    bool initialized;
    int AerialImgWidth, AerialImgHeight;

public:
    thermalArea();
    void setupAOI(int AerialImgWidth, int AerialImgHeight, float AerialFOV);
    bool isInAOI(int x, int y);
    bool isInitialized();
};

#endif //MULTICONTOUR_THERMALAREA_H
