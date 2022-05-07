#include "matrix.h"
#include "viewer.h"
#include <string>
#include "Macros.h"
static const char* const WINDOW_TITLE = "Viewer";


//========================================================================
int main() {
	Viewer Mviewer;
	Mviewer.launch_init(WINDOW_TITLE, DEFAULT_WINDOW_WIDTH, DEFAULT_WINDOW_HEIGHT, false, false);
	Mviewer.launch_rendering<true>();
	Mviewer.launch_shut();
	return EXIT_SUCCESS;
}
