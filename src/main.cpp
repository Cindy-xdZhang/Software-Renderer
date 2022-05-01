#include "matrix.h"
#include "Base.h"
#include "viewer.h"
#include <string>
static const char* const WINDOW_TITLE = "Viewer";
//static const int WINDOW_WIDTH = 800;
//static const int WINDOW_HEIGHT = 600;
//static const int CAMERA_POSITION[] = { 0, 0, 1.5f };
//static const int CAMERA_TARGET[] = { 0, 0, 0 };
//static const float LIGHT_THETA = TO_RADIANS(45);
//static const float LIGHT_PHI = TO_RADIANS(45);
//static const float LIGHT_SPEED = PI;
//static const float CLICK_DELAY = 0.25f;
//




//========================================================================
int main() {
	Viewer Mviewer;
	Mviewer.launch_init(WINDOW_TITLE,320,480, false, false, false);
	Mviewer.launch_rendering();
	Mviewer.launch_shut();
	return EXIT_SUCCESS;
}

