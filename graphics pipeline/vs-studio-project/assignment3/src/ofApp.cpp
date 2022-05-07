#include "obj.hpp"
#include "objparser.h"
#include "ofApp.h"
#include "Base.h"

#include <malloc.h>


GraphicsPipeline GP;
RenderObject Sample;


//--------------------------------------------------------------
void ofApp::setup() {
	w = WINDOW_W;
	h = WINDOW_H;
	colorPixels.allocate(w, h, OF_PIXELS_RGB);
	OBjReader objFilereader;
	string Path = string("C:\\Users\\8\\source\\repos\\assignment3\\assignment3\\KAUST_Beacon_NormalTexture.obj");
	objFilereader.read(Path);
	
	Sample=RenderObject(objFilereader.TrianglesIdx, objFilereader.Vertexes, objFilereader.VertexesNormal, objFilereader.VertexesTexture);
	GP = GraphicsPipeline(h, w);
}

//--------------------------------------------------------------
void ofApp::update() {
	static int a=0;
	float angle = Radians(a) ;
	GP.InitLightPos = { 125*cos(angle) ,55 ,125 * sin(angle) };
	a += 5;
	a = a % 360;
}

//--------------------------------------------------------------
void ofApp::draw() {
	mVec3* framebuffer = (mVec3*)malloc(sizeof(mVec3) * w*h);
	GP.Render(Sample, framebuffer);
	// color pixels, use x and y to control red and green
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			colorPixels.setColor(x, y, ofColor(framebuffer[y*w + x].x * 255, framebuffer[y*w + x].y * 255, framebuffer[y*w + x].z * 255));
		}
	}

	texColor.allocate(colorPixels);
	ofSetHexColor(0xffffff);
	texColor.draw(0, 0, w, h);
	free(framebuffer);
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {

	if (key == 'd') {
		m_g_pip->_Camera.UpdatePos(mVec3(5,0,0));
	}	
	if (key == 'a') {
		m_g_pip->_Camera.UpdatePos(mVec3(-5, 0, 0));
	}
	if (key == 's') {
		m_g_pip->_Camera.UpdatePos(mVec3(0, -5, 0));
	}
	if (key == 'w') {
		m_g_pip->_Camera.UpdatePos(mVec3(0, 5, 0));
	}
	if (key == 'q') {
		GP.InitLightPos.z -=5;
	}
	if (key == 'e') {//go outside 
		m_g_pip->_Camera.UpdatePos(mVec3(0, 0, -5));
	}
	if (key == 'r') {//go into 
		m_g_pip->_Camera.UpdatePos(mVec3(0, 0, 5));
	}
	if (key == 'b') {//go into 
		GP.UsePhongShading = !GP.UsePhongShading;
	}
	if (key == 'n') {//go into 
		GP.Texture = (GP.Texture+1)%3;
	}
	/*if (key == 'x') {
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
	}*/
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y) {

}

//--------------------------------------------------------------
static bool dragging = false;
void ofApp::mouseReleased(int x, int y, int button) {
	dragging = false;
}
void ofApp::mouseDragged(int x, int y, int button) {
	static Matrix4 LastRot = eye(4);
	static Matrix4 ThisRot = eye(4);
	static int lastX = 0, lastY = 0;
	static mVec3 Op1, Op2;
	static ArcBallControler ArcControl = ArcBallControler(w, h);
	if (button == 2) {
		LastRot = eye(4);
		ThisRot = eye(4);
		Sample.ModleMatrix = ThisRot * translateMatrix(mVec3(-125, -125, -125));
	}
	else if (button == 0) {
		if (dragging) {
			cout << "mouseDragged: y=" << y - lastY << " x=" << x - lastX << " button=" << button << endl;//鼠标的x,y坐标
			Op2 = ArcControl.GetArcBallPositionVector(x, y);
			ThisRot = ArcControl.GetArcBallrotateMatrix(Op1, Op2);
			ThisRot = ThisRot * LastRot;
			Sample.ModleMatrix = ThisRot * translateMatrix(mVec3(-125, -125, -125));


		}
		else {
			lastX = x;
			lastY = y;
			Op1 = ArcControl.GetArcBallPositionVector(lastX, lastY);
			dragging = true;
			LastRot = ThisRot;
		}
	}
}
void ofApp::mouseScrolled(int x, int y, float scrollX, float scrollY) {
	cout << "MouseScrolled: x=" << y << " y=" << x << " scrollX=" << scrollX << " scrollY=" << scrollY << endl;//鼠标的x,y坐标
	m_g_pip->_Camera.UpdateFov(-scrollY);

}
//--------------------------------------------------------------

void ofApp::mousePressed(int x, int y, int button) {

}

//--------------------------------------------------------------

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
