#ifndef MCA_CONFIG_H
#define MCA_CONFIG_H

// Define the max number of images that will be analyzed in a fireline progression analysis
#define MAX_IMAGE_SERIES  750
#define FILE_PATH_SIZE    150


class MCA_Config
{
    private:
        char  analysisName[FILE_PATH_SIZE];
        char  inputDirectory[FILE_PATH_SIZE];
        char  outputDirectory[FILE_PATH_SIZE];
        int   imageryWidth  = 760;                       // targeted scaled pixel width of thermal and rgb imagery
        int   imageryHeight = 1024;                      // targeted scaled pixel height of thermal and rgb imagery
        int   fileCount;                                 // total number of images in the series
        float XScale;                                    // coverage width of imagery in meters
        float YScale;                                    // coverage height of imagery in meters
        float FOV;                                       // horizontal field of view in meters
        float RadiometerFOV;                             // Radiometer pixel diameter in meters
        float elapsedTime;                               // elapsed time between photos
        int   fireLine_MinContourPoints;                   // Min. contour point count required for a fire line width analysis
        float fireLine_MinContourArea;                   // Min. contour area required for a fire line width analysis
        float setSize;                                   // # of images averaged together to provide single datapoint
        int   normalInterval;                              // # of adjacent points used on a contour level to fit a line and its normal vector
        char rgb_File[MAX_IMAGE_SERIES][FILE_PATH_SIZE];
        char ir_file[MAX_IMAGE_SERIES][FILE_PATH_SIZE];
        char configFileName[FILE_PATH_SIZE];

    public:
        MCA_Config();
        int Read_MCA_Config(char *configFile);
        char *getAnalysisName();
        char *getInputDirectory();
        char *getOutputDirectory();
        int  getFileCount();
        char *getRGB_Input(int i);
        char *getIR_Input(int i);
        char *getConfigFileName();

        float getFireLine_MinContourPoints();
        float getFireLine_MinContourArea();
        float getNormalInterval();
        float getSetSize();
        float getElapsedTime();
        float getFOV();
        float getRadiometerFOV();
        int getImageryWidth();
        int getImageryHeight();
        //float getXScale();
        //float getYScale();
        void printConfig();
};

#endif // MCA_CONFIG_H
