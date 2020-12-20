#include "ofApp.h"
#include "Base.h"
int gcount = 1;
bool  UsePerspectiveProjection = true;
extern mVec3 lightPos;
//extern  mVec3 eyepos;
extern  int MousePosition_x;
extern int MousePosition_y;
extern Object_List objlist;
extern int MouseSelectObjID ;
extern Camera* MyCamera;
int animation_time=0;
float speed = 0.1;
bool key_ctrl_pressed = false;
bool key_F1_pressed = false;
bool animation = false;
//--------------------------------------------------------------
void ofApp::setup() {
	w = WINDOW_W;
	h = WINDOW_H;
	colorPixels.allocate(w, h, OF_PIXELS_RGB);
	lightPos = mVec3(-5, 0, 0);
	setupSence();
}


//--------------------------------------------------------------
void ofApp::update() {
	if (gcount == 1) {
		lightPos = lightPos + mVec3(-0.05, 0.2, 0);
		if (lightPos.x < -7) gcount = 0;
	}
	else {
		lightPos = lightPos - mVec3(-0.05, 0.2, 0);
		if (lightPos.x > -3) gcount = 1;
	}
	if (animation) {
		float angle = animation_time * 2.0f;	
		animation_time++;
		MyCamera->position = mVec3(0, 5*sin(Radians(angle)), -5*cos( Radians(angle) )    );
		MyCamera->UpdateTarget(mVec3(0, 0, 0));
		auto objp = (Surface*)objlist.GetPointer(0);
		mVec3 translation(sin(Radians(angle*5))*0.15 , 0, 0);
		objp->updateTranslateMatrix(translateMatrix(translation));

		 objp = (Surface*)objlist.GetPointer(5);
		 Matrix4 r = rotateMatrix(2, speed * 50);
		 objp->updateRotateMatrix(r);
	}
}

//--------------------------------------------------------------
void ofApp::draw() {
	FloatRGB* framebuffer = new FloatRGB[w*h];
   multithread_simple_ray_tracing(framebuffer, UsePerspectiveProjection,lightPos);
	//simple_ray_tracing(framebuffer, true, lightPos);
	// color pixels, use x and y to control red and green
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			colorPixels.setColor(x, y, ofColor(framebuffer[y*w + x].r * 255, framebuffer[y*w + x].g * 255, framebuffer[y*w + x].b * 255));
		}
	}

	texColor.allocate(colorPixels);
	ofSetHexColor(0xffffff);
	texColor.draw(0, 0, w, h);
	delete framebuffer;
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {
	if (key == 't') {
		UsePerspectiveProjection = !UsePerspectiveProjection;
	}
	if (key == 'u') {
		if (MouseSelectObjID != -1) {
			auto objp=(Surface*)objlist.GetPointer(MouseSelectObjID);
			mVec3 translation(-1,0,0);
			translation = translation * speed;
			objp->updateTranslateMatrix(translateMatrix(translation));
		}
	}	if (key == 'j') {
		if (MouseSelectObjID != -1) {
			auto objp = (Surface*)objlist.GetPointer(MouseSelectObjID);
			mVec3 translation(1, 0, 0);
			translation = translation * speed;
			objp->updateTranslateMatrix(translateMatrix(translation));
		}
	}
	if (key == 'w') {
		if (key_F1_pressed) {
			MyCamera->UpdatePitchAngle(speed * -20);
		}
		else if (MouseSelectObjID != -1) {
				auto objp = (Surface*)objlist.GetPointer(MouseSelectObjID);
				mVec3 translation(0, 0, 1);
				translation = translation * speed;
				objp->updateTranslateMatrix(translateMatrix(translation));
		}
	}	if (key == 's') {
		if (key_F1_pressed) {
			MyCamera->UpdatePitchAngle(speed * 20);
		}
		else if (MouseSelectObjID != -1) {
				auto objp = (Surface*)objlist.GetPointer(MouseSelectObjID);
				mVec3 translation(0, 0, -1);
				translation = translation * speed;
				objp->updateTranslateMatrix(translateMatrix(translation));
			
		}
	}
	if (key == 'a') {
		if (key_F1_pressed) {
			MyCamera->UpdateYawAngle(speed * 20);
		}
		else if (MouseSelectObjID != -1) {
			auto objp = (Surface*)objlist.GetPointer(MouseSelectObjID);
			mVec3 translation(0, -1, 0);
			translation = translation * speed;
			objp->updateTranslateMatrix(translateMatrix(translation));
		}
	}	if (key == 'd') {
		if (key_F1_pressed) {
			MyCamera->UpdateYawAngle(speed * -20);
		}
		else if (MouseSelectObjID != -1) {
			auto objp = (Surface*)objlist.GetPointer(MouseSelectObjID);
			mVec3 translation(0, 1, 0);
			translation = translation * speed;
			objp->updateTranslateMatrix(translateMatrix(translation));
		}
	}
	if (key == 'x') {
		if (MouseSelectObjID != -1) {
			if (!key_ctrl_pressed) {
				auto objp = (Surface*)objlist.GetPointer(MouseSelectObjID);
				Matrix4 r = rotateMatrix(0, speed * 50);
				objp->updateRotateMatrix(r);
			}
			else {
				auto objp = (Surface*)objlist.GetPointer(MouseSelectObjID);
				Matrix4 s = scaleMatrix(1.1f, 1.0f, 1.0f);
				objp->updateScaleMatrix(s);
			}
		
		}
	}
	if (key == 'y') {
		if (!key_ctrl_pressed) {
			auto objp = (Surface*)objlist.GetPointer(MouseSelectObjID);
			Matrix4 r = rotateMatrix(1, speed * 50);
			objp->updateRotateMatrix(r);
		}
		else {
			auto objp = (Surface*)objlist.GetPointer(MouseSelectObjID);
			Matrix4 s = scaleMatrix(1.0f, 1.1f, 1.0f);
			objp->updateScaleMatrix(s);
		}
	}
	if (key == 'z') {
		if (!key_ctrl_pressed) {
			auto objp = (Surface*)objlist.GetPointer(MouseSelectObjID);
			Matrix4 r = rotateMatrix(2, speed * 50);
			objp->updateRotateMatrix(r);
		}
		else {
			auto objp = (Surface*)objlist.GetPointer(MouseSelectObjID);
			Matrix4 s = scaleMatrix(1.0f, 1.0f, 1.1f);
			objp->updateScaleMatrix(s);
		}
	}
	if (key == 'r') {
		if (MouseSelectObjID != -1) {
			auto objp = (Surface*)objlist.GetPointer(MouseSelectObjID);
			objp->ResetModelMatrix();
		}
	}
	if (key == OF_KEY_CONTROL) {
		key_ctrl_pressed = !key_ctrl_pressed;
		if (key_ctrl_pressed)cout << " scale!" << endl;
		else cout << " rotate!" << endl;
	}
	if (key == OF_KEY_LEFT) {
		key_ctrl_pressed = !key_ctrl_pressed;
		MyCamera->UpdatePos(mVec3(0,-1,0)*speed);
			
	}
	if (key == OF_KEY_RIGHT) {
		key_ctrl_pressed = !key_ctrl_pressed;
		MyCamera->UpdatePos(mVec3(0, 1, 0)*speed);
	
	}
	if (key == OF_KEY_UP) {
		key_ctrl_pressed = !key_ctrl_pressed;
		MyCamera->UpdatePos(mVec3(-1, 0, 0)*speed);
	
	}
	if (key == OF_KEY_DOWN) {
		key_ctrl_pressed = !key_ctrl_pressed;
		MyCamera->UpdatePos(mVec3(1, 0, 0)*speed);
	}
	if (key == OF_KEY_PAGE_UP) {
		key_ctrl_pressed = !key_ctrl_pressed;
		MyCamera->position += (mVec3(0, 0, 1)*speed);

	}
	if (key == OF_KEY_PAGE_DOWN) {
		key_ctrl_pressed = !key_ctrl_pressed;
		MyCamera->position += (mVec3(0, 0, -1)*speed);

	}
	if (key == OF_KEY_F1) {
		key_F1_pressed = !key_F1_pressed;
		if (key_F1_pressed)cout << "CAMERA ROTATE!" << endl;
		else cout << "  OBJ MOTION!" << endl;

	}
	if (key == 'v') {
		animation = true;
	}
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button) {

}
void ofApp::mouseScrolled(int x, int y, float scrollX, float scrollY) {
	cout << "MouseScrolled: x=" << y << " y=" << x << " scrollX=" << scrollX << " scrollY=" << scrollY << endl;//鼠标的x,y坐标
	MyCamera->UpdateFov(scrollY);

}
//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button) {
	cout <<"x="<< y<< " y=" << x << " button=" << button<<endl;//鼠标的x,y坐标
	MousePosition_x = y;
	MousePosition_y = x;
	if (button == 2) {
		animation = false;
		MyCamera->Reset(); 
	
	}
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h) {

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg) {

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo) {

}
