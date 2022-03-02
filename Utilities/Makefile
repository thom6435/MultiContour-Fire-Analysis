clean:
	rm *.o ClipAlign rgbClipAlign thermalClipAlign
ClipAlign: ClipAlign.o
	g++ -o ClipAlign ClipAlign.o -lopencv_highgui -lopencv_core -lopencv_imgcodecs -lopencv_imgproc
ClipAlign.o: ClipAlign.cpp
	g++ -Wall -c ClipAlign.cpp
rgbClipAlign: rgbClipAlign.o
	g++ -o rgbClipAlign rgbClipAlign.o -lopencv_highgui -lopencv_core -lopencv_imgcodecs -lopencv_imgproc
rgbClipAlign.o: rgbClipAlign.cpp
	g++ -Wall -c rgbClipAlign.cpp
thermalClipAlign: thermalClipAlign.o
	g++ -o thermalClipAlign thermalClipAlign.o -lopencv_highgui -lopencv_core -lopencv_imgcodecs -lopencv_imgproc
thermalClipAlign.o: thermalClipAlign.cpp
	g++ -Wall -c thermalClipAlign.cpp
thermalImageScale: thermalImageScale.o
	g++ -o thermalImageScale thermalImageScale.o -lopencv_highgui -lopencv_core -lopencv_imgcodecs -lopencv_imgproc
thermalImageScale.o: thermalImageScale.cpp
	g++ -Wall -c thermalImageScale.cpp
ThermalImageScaleRotate: ThermalImageScaleRotate.o
	g++ -o ThermalImageScaleRotate ThermalImageScaleRotate.o -lopencv_highgui -lopencv_core -lopencv_imgcodecs -lopencv_imgproc
ThermalImageScaleRotate.o: ThermalImageScaleRotate.cpp
	g++ -Wall -c ThermalImageScaleRotate.cpp
