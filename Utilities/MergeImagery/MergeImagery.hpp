#ifndef MULTICONTOUR_HPP
#define MULTICONTOUR_HPP

// Contour parameters
#define MAXCONTOURS       2
#define CONTOUR_INTERVAL 30
#define MAXREGIONS       32
#define MAX_AVG_SET       2
#define PATH_SIZE       256

// Contour drawing color definitions
// ColorSpace if Blue, Green, Red
#define RED1	       cv::Scalar(0x00, 0x00, 0xcc)
#define RED2	       cv::Scalar(0x11, 0x11, 0xcc)
#define ORANGE1        cv::Scalar(0x00, 0x66, 0xcc)
#define ORANGE2        cv::Scalar(0x33, 0x99, 0xcc)
#define YELLOWORANGE1  cv::Scalar(0x00, 0xcc, 0xff)
#define YELLOWORANGE2  cv::Scalar(0x33, 0xff, 0xff)
#define YELLOW1        cv::Scalar(0x00, 0xcc, 0xcc)
#define YELLOW2        cv::Scalar(0x66, 0xff, 0xcc)
#define YELLOWGREEN1   cv::Scalar(0x66, 0xff, 0x99)
#define YELLOWGREEN2   cv::Scalar(0x66, 0xcc, 0x99)
#define GREEN1         cv::Scalar(0x00, 0xcc, 0x66)
#define GREEN2         cv::Scalar(0x66, 0xcc, 0x66)
#define VIOLET1        cv::Scalar(0x99, 0x00, 0x99)
#define VIOLET2        cv::Scalar(0x99, 0x33, 0x99)
#define VIOLETBLUE1    cv::Scalar(0x99, 0x33, 0x66)
#define VIOLETBLUE2    cv::Scalar(0x99, 0x66, 0x66)
#define VIOLETBLUE3    cv::Scalar(0xcc, 0x66, 0x66)
#define VIOLETBLUE4    cv::Scalar(0xcc, 0x99, 0x66)
#define BLUE1          cv::Scalar(0xff, 0x99, 0x66)
#define BLUE2          cv::Scalar(0xff, 0x33, 0x66)
#define BLUE3          cv::Scalar(0x99, 0x33, 0x33)
#define BLUE4          cv::Scalar(0x66, 0x33, 0x33)
#define BLUEBLACK      cv::Scalar(0x33, 0x00, 0x00)
#define BLACK          cv::Scalar(0x00, 0x00, 0x00)

#endif // MULTICONTOUR_HPP
