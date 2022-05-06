#pragma once
#define EPSILON 1e-5f
#define PI 3.141592653589793238
#define FLOAT_ZERO 1e-18



#define TO_RADIANS(degrees) ((PI / 180) * (degrees))
#define TO_DEGREES(radians) ((180 / PI) * (radians))

#define LINE_SIZE 256
#define PATH_SIZE 256

#define UNUSED_VAR(x) ((void)(x))
#define ARRAY_SIZE(array) (sizeof(array) / sizeof((array)[0]))



#define MAX_DIMENTION  10
//#define MAIN_DATA_PRECISION float
#define Radians(x) ((x)*PI/180.0)
#define MAX(a,b)    (((a)>(b)) ?  (a):(b))
#define MIN(a,b)    (((a)>(b)) ?  (b):(a))

#define DEFAULT_WINDOW_HEIGHT 960
#define DEFAULT_WINDOW_WIDTH 640
