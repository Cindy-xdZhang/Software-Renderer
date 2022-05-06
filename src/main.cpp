#include "matrix.h"
#include "viewer.h"
#include <string>
#include "Macros.h"
static const char* const WINDOW_TITLE = "Viewer";


//========================================================================
int main() {
	Viewer Mviewer;
	Mviewer.launch_init(WINDOW_TITLE,DEFAULT_WINDOW_HEIGHT, DEFAULT_WINDOW_WIDTH, false, false);
	Mviewer.launch_rendering();
	Mviewer.launch_shut();
	return EXIT_SUCCESS;
}

