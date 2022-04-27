#ifndef OF_MAIN_H
#define OF_MAIN_H

//--------------------------
// utils
//#include "ofConstants.h"
//#include "ofFileUtils.h"
//#include "ofLog.h"
//#include "ofSystemUtils.h"
//
//#include "ofURLFileLoader.h"
//
//#include "ofUtils.h"

//#if !defined(TARGET_EMSCRIPTEN)
//#include "ofThread.h"
//#include "ofThreadChannel.h"
//#endif

//#include "ofFpsCounter.h"
//#include "ofJson.h"
//#include "ofXml.h"

//--------------------------
// types
//#include "ofGraphicsBaseTypes.h"
//#include "ofTypes.h"
//#include "ofColor.h"
//#include "ofPoint.h"
//#include "ofRectangle.h"
//#include "ofParameter.h"
//#include "ofParameterGroup.h"

//--------------------------
// app
#include "ofBaseApp.h"
#include "ofAppRunner.h"
#include "ofAppBaseWindow.h"
#include "ofWindowSettings.h"
#include "ofMainLoop.h"
#if !defined( TARGET_OF_IOS ) & !defined(TARGET_ANDROID) & !defined(TARGET_EMSCRIPTEN) & !defined(TARGET_RASPBERRY_PI_LEGACY)
	#include "ofAppGLFWWindow.h"
	#if !defined( TARGET_LINUX_ARM )
		#include "ofAppGlutWindow.h"
	#endif
#endif


//--------------------------
using namespace std;

#endif
