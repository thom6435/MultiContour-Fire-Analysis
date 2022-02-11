#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mca_config.h"


MCA_Config::MCA_Config()
{
    this->fileCount = 0;
}


/*
 * function: Read_MCA_Config
 * modified: 3-August-2020
 * returns:  0 == read error/failure
 *           1 == success
 */
int MCA_Config::Read_MCA_Config(char *configFile)
{
    FILE *configStream;
    char inputStr[256], identifier[32];
    bool done;
    unsigned long int i, j, k;

    strcpy(this->getConfigFileName(), configFile);

    if (!(configStream = fopen(configFile,"r"))) {
        printf("Error Unable to open setup file %s\n", configFile);
        return 0;
    }

    while(!feof(configStream)) {
        fgets(inputStr, 256, configStream);
        if (inputStr[0] != '#') {

            // Build Identifier Name

            done = false;
            i = 0;
            while (!done && i < strlen(inputStr)) {
                if (inputStr[i] ==':') {
                    done = true;
                    identifier[i] = 0;
                } else {
                    identifier[i] = inputStr[i];
                }
                i++;
            } // !done

            if (strcmp(identifier,"ANALYSIS_NAME") == 0) {
                k=0;
                for (j=i;j<strlen(inputStr); j++){
                    if (inputStr[j] > 32)
                        this->analysisName[k++]= inputStr[j];
                }
                this->analysisName[k]=0;
            }

            if (strcmp(identifier,"OUTPUT_DIRECTORY") == 0) {
                k=0;
                for (j=i;j<strlen(inputStr); j++){
                    if (inputStr[j] > 32)
                        this->outputDirectory[k++]= inputStr[j];
                }
                this->outputDirectory[k]=0;
            }

            if (strcmp(identifier,"INPUT_DIRECTORY") == 0) {
                k=0;
                for (j=i;j<strlen(inputStr); j++){
                    if (inputStr[j] > 32)
                        this->inputDirectory[k++]= inputStr[j];
                }
                this->inputDirectory[k]=0;
            }

            if (strcmp(identifier,"X_SCALE") == 0) {
                char buffer[16];
                k=0;
                for (j=i;j<strlen(inputStr); j++){
                    if (inputStr[j] > 32)
                       buffer[k++]= inputStr[j];
                }
                buffer[k]=0;

                this->XScale = atof(buffer);
            }

            if (strcmp(identifier,"Y_SCALE") == 0) {
                char buffer[16];
                k=0;
                for (j=i;j<strlen(inputStr); j++){
                    if (inputStr[j] > 32)
                       buffer[k++]= inputStr[j];
                }
                buffer[k]=0;

                this->YScale = atof(buffer);
            }

	    if (strcmp(identifier,"ELAPSED_TIME") == 0) {
	      char buffer[16];
	      k=0;
	      for (j=i;j<strlen(inputStr); j++){
		if (inputStr[j] > 32)
		  buffer[k++]= inputStr[j];
	      }
	      buffer[k]=0;
	      
	      this->elapsedTime = atof(buffer);
	    }

	    if (strcmp(identifier,"SET_SIZE") == 0) {
	      char buffer[16];
	      k=0;
	      for (j=i;j<strlen(inputStr); j++){
		if (inputStr[j] > 32)
		  buffer[k++]= inputStr[j];
	      }
	      buffer[k]=0;
	      
	      this->setSize = atof(buffer);
	    }
	    
            if (strcmp(identifier,"INPUT_FILE_SERIES") == 0) {
                puts("Parsing input file series...");
                k=0;
                while(!feof(configStream)) {
                    fgets(inputStr, 256, configStream);
                    sscanf(inputStr,"%s %s", this->rgb_File[k], this->ir_file[k]);
                    k++;
                }
                this->fileCount = k;
            }
        }
    } // !feof

    return 1;
} //Read_MCA_Config


char *MCA_Config::getAnalysisName() {
    return this->analysisName;
}


char *MCA_Config::getInputDirectory() {
    return this->inputDirectory;
}


char *MCA_Config::getOutputDirectory() {
    return this->outputDirectory;
}


int  MCA_Config::getFileCount() {
    return this->fileCount;
}


char *MCA_Config::getRGB_Input(int i){
    return this->rgb_File[i];
}


char *MCA_Config::getIR_Input(int i) {
    return this->ir_file[i];
}


char *MCA_Config::getConfigFileName() {
    return configFileName;
}

float MCA_Config::getSetSize() {
  return this->setSize;
}


float MCA_Config::getElapsedTime() {
  return this->elapsedTime;
}


float MCA_Config::getXScale() {
  return this->XScale;
}


float MCA_Config::getYScale() {
  return this->YScale;
}


void MCA_Config::printConfig() {
    int i;

    printf(">>--------> Configuration Name: %s\n", this->configFileName);
    printf("                 Analysis Name: %s\n", this->analysisName);
    printf("               Input Directory: %s\n", this->inputDirectory);
    printf("              Output Directory: %s\n", this->outputDirectory);
    printf("             Image Scale Width: %6.1f\n", this->XScale);
    printf("            Image Scale Height: %6.1f\n", this->YScale);
    printf("   Elapsed time between Images: %6.1f\n", this->elapsedTime);
    printf("                Image Set Size: %2.0f\n", this->setSize);
    printf("    Number of Images in Series: %3d\n", this->fileCount);
    printf("Input File Series:\n");
    for (i=0; i<this->fileCount; i++){
        printf("   %3d:  RGB: %s\n", i, this->rgb_File[i]);
        printf("          IR: %s\n", this->ir_file[i]);
    }

} // printConfig
