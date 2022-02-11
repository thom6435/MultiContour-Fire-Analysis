#ifndef MCA_CONFIG_H
#define MCA_CONFIG_H

// Define the max number of images that will be analyzed in a fireline progression analysis
#define MAX_IMAGE_SERIES  100
#define FILE_PATH_SIZE    128


class MCA_Config
{
    private:
        char analysisName[FILE_PATH_SIZE];
        char inputDirectory[FILE_PATH_SIZE];
        char outputDirectory[FILE_PATH_SIZE];
        int  fileCount;                                  // total number of images in the series
        float XScale;                                    // coverage width of imagery in meters
        float YScale;                                    // coverage height of imagery in meters
        float elapsedTime;                               // elapsed time between photos 
        float setSize;                                   // # of images averaged together to provide single datapoint
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
        float getSetSize();
        float getElapsedTime();
        float getXScale();
        float getYScale();
        void printConfig();
};

#endif // MCA_CONFIG_H
