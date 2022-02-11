MergeImagery: MergeImagery.o contourregion.o fireline.o imageanalysis.o mca_config.o normal.o
	g++ -o MergeImagery MergeImagery.o contourregion.o fireline.o imageanalysis.o mca_config.o normal.o -lm -lopencv_highgui -lopencv_core -lopencv_xphoto -lopencv_imgcodecs -lopencv_imgproc
MergeImagery.o: MergeImagery.cpp MergeImagery.hpp fireline.h contourregion.h normal.h mca_config.h imageanalysis.h
	g++ -Wall -c MergeImagery.cpp 
contourregion.o: contourregion.cpp contourregion.h
	g++ -Wall -c contourregion.cpp
fireline.o: fireline.cpp fireline.h
	g++ -Wall -c fireline.cpp
imageanalysis.o: imageanalysis.cpp imageanalysis.h
	g++ -Wall -c imageanalysis.cpp
mca_config.o: mca_config.cpp mca_config.h
	g++ -Wall -c mca_config.cpp
normal.o: normal.cpp normal.h contourregion.h
	g++ -Wall -c normal.cpp 
clean:
	rm *.o MergeImagery
 
