#include "ofMain.h"
#include "ofApp.h"
#include "Base.h"
#include "ofAppGlutWindow.h"
//#include "matrix.h"
#include "objparser.h"


//========================================================================
int main() {
	ofSetupOpenGL(WINDOW_W, WINDOW_H, OF_WINDOW);		// <-------- setup the GL context Sets up the window aspect and mode.	

	// this kicks off the running of my app
	// can be OF_WINDOW or OF_FULLSCREEN
	// pass in width and height too:
	ofRunApp(new ofApp());//Begins the openGL cycle of the application.
	//It's only called once from main function in main.cpp after setting the window with ofSetupOpenGL.


}
//int main() {
//	OBjReader objFilereader; 
//	string Path = string("C:\\Users\\8\\source\\repos\\assignment3\\assignment3\\KAUST_Beacon_Normal.obj");
//	string outPath = string("C:\\Users\\8\\source\\repos\\assignment3\\assignment3\\KAUST_Beacon_NormalTexture.obj");
//	//objFilereader.read(Path);
//	//objFilereader.Load_CalculateNormal_Store(Path, outPath);
//	objFilereader.Load_CalculateTextureCoordinate_Store(Path, outPath);
//	//RenderObject Sample(objFilereader.TrianglesIdx, objFilereader.Vertexes);
//	//GraphicsPipeline GP;
//	//GP.Render(Sample);
//
//	cout << "";
//}

