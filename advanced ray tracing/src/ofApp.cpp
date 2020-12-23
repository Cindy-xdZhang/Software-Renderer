#include "ofApp.h"
#include "Base.h"


//extern mVec3 lightPos;
extern  int MousePosition_x;
extern int MousePosition_y;
extern Object_List objlist;
int MouseSelectObjID ;
extern Camera* MyCamera;
extern BvhTopLevelStructure* my;
extern bool TureisUsePredefinedNormalFalseIsUseAverageNormal;
extern bool UseDistributedRaySSP ;
extern bool UseDistributedRaySoftShadow;
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
	setupSence();
	MouseSelectObjID = -1;

}


//--------------------------------------------------------------
void ofApp::update() {
	std::stringstream strm;
	strm << "fps: " << ofGetFrameRate();
	ofSetWindowTitle(strm.str());

	
		//lightPos = lightPos + mVec3(-0.05, 0.0, 0);


	
}

//--------------------------------------------------------------
void ofApp::draw() {
	mVec3* framebuffer = new mVec3[w*h];
   multithread_simple_ray_tracing(framebuffer);
  // simple_ray_tracing(framebuffer, true, lightPos);
	// color pixels, use x and y to control red and green
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			colorPixels.setColor(x, y, ofColor(framebuffer[y*w + x].x * 255, framebuffer[y*w + x].y * 255, framebuffer[y*w + x].z* 255));
		}
	}

	texColor.allocate(colorPixels);
	ofSetHexColor(0xffffff);
	texColor.draw(0, 0, w, h);
	delete[] framebuffer;
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {
	if (key == 't') {
		TureisUsePredefinedNormalFalseIsUseAverageNormal = !TureisUsePredefinedNormalFalseIsUseAverageNormal;
		if(TureisUsePredefinedNormalFalseIsUseAverageNormal)
		cout << "Use Predefined Normal!" << endl;
		else		cout << "Use Perface Normal!" << endl;
	}
	if (key == 'm') {
		UseDistributedRaySSP = !UseDistributedRaySSP;
	}
	if (key == 'b') {
		UseDistributedRaySoftShadow = !UseDistributedRaySoftShadow;
	}
	if (key == 'u') {
		if (MouseSelectObjID != -1) {
			auto objp = (Surface*)&(my->Meshes[MouseSelectObjID]);
			mVec3 translation(-1,0,0);
			translation = translation * speed;
			objp->updateModelMatrix(translateMatrix(translation));
			my->UpdateBoundingBOx();
		}
	}	
	if (key == 'j') {
		if (MouseSelectObjID != -1) {
			auto objp = (Surface*)&(my->Meshes[MouseSelectObjID]);
			mVec3 translation(1, 0, 0);
			translation = translation * speed;
			objp->updateModelMatrix(translateMatrix(translation));
			my->UpdateBoundingBOx();
		}
	}
	if (key == 'w') {
		if (key_F1_pressed) {
			MyCamera->UpdatePitchAngle(speed * -5);
		}
		else if (MouseSelectObjID != -1) {
			auto objp = (Surface*)&(my->Meshes[MouseSelectObjID]);
				mVec3 translation(0, 0, 1);
				translation = translation * speed;
				objp->updateModelMatrix(translateMatrix(translation));
				my->UpdateBoundingBOx();
		}
	}	if (key == 's') {
		if (key_F1_pressed) {
			MyCamera->UpdatePitchAngle(speed * 5);
		}
		else if (MouseSelectObjID != -1) {
			auto objp = (Surface*)&(my->Meshes[MouseSelectObjID]);
				mVec3 translation(0, 0, -1);
				translation = translation * speed;
				objp->updateModelMatrix(translateMatrix(translation));
				my->UpdateBoundingBOx();
		}
	}
	if (key == 'a') {
		if (key_F1_pressed) {
			MyCamera->UpdateYawAngle(speed * 5);
		}
		else if (MouseSelectObjID != -1) {
			auto objp = (Surface*)&(my->Meshes[MouseSelectObjID] );
			mVec3 translation(0, -1, 0);
			translation = translation * speed;
			objp->updateModelMatrix(translateMatrix(translation));
			my->UpdateBoundingBOx();
		}
	}	if (key == 'd') {
		if (key_F1_pressed) {
			MyCamera->UpdateYawAngle(speed * -5);
		}
		else if (MouseSelectObjID != -1) {
			auto objp = (Surface*)&(my->Meshes[MouseSelectObjID]);
			mVec3 translation(0, 1, 0);
			translation = translation * speed;
			objp->updateModelMatrix(translateMatrix(translation));
			my->UpdateBoundingBOx();
		}
	}
	if (key == 'x') {
		if (MouseSelectObjID != -1) {
			if (!key_ctrl_pressed) {
				auto objp = (Surface*)&(my->Meshes[MouseSelectObjID]);
				Matrix4 r = rotateMatrix(0, speed * 50);
				Matrix4 realM = objp->pModleMatrix  * r * objp->pIModleMatrix;
				objp->updateModelMatrix(realM);
				my->UpdateBoundingBOx();
			}
			else {
				auto objp = (Surface*)&(my->Meshes[MouseSelectObjID]);
				Matrix4 s = scaleMatrix(1.1f, 1.0f, 1.0f);
				Matrix4 realM = objp->pModleMatrix  * s * objp->pIModleMatrix;
				objp->updateModelMatrix(realM);
				my->UpdateBoundingBOx();
			}
		
		}
	}
	if (key == 'y') {
		if (!key_ctrl_pressed) {
			auto objp = (Surface*)&(my->Meshes[MouseSelectObjID]);
			Matrix4 r = rotateMatrix(1, speed * 50);
			Matrix4 realM = objp->pModleMatrix  * r * objp->pIModleMatrix;
			objp->updateModelMatrix(realM);
			my->UpdateBoundingBOx();
		}
		else {
			auto objp = (Surface*)&(my->Meshes[MouseSelectObjID]);
			Matrix4 s = scaleMatrix(1.0f, 1.1f, 1.0f);
			Matrix4 realM = objp->pModleMatrix  * s * objp->pIModleMatrix;
			objp->updateModelMatrix(realM);
			my->UpdateBoundingBOx();
		}
	}
	if (key == 'z') {
		if (!key_ctrl_pressed) {
			auto objp = (Surface*)&(my->Meshes[MouseSelectObjID]);
			Matrix4 r = rotateMatrix(2, speed * 50);
			Matrix4 realM = objp->pModleMatrix  * r * objp->pIModleMatrix;
			objp->updateModelMatrix(realM);
			my->UpdateBoundingBOx();
		}
		else {
			auto objp = (Surface*)&(my->Meshes[MouseSelectObjID]);
			Matrix4 s = scaleMatrix(1.0f, 1.0f, 1.1f);
			Matrix4 realM = objp->pModleMatrix  * s * objp->pIModleMatrix;
			objp->updateModelMatrix(realM);
			my->UpdateBoundingBOx();
		}
	}
	if (key == 'r') {
		if (MouseSelectObjID != -1) {
			auto objp = (Surface*)&(my->Meshes[MouseSelectObjID]);
			objp->initModelMatrix();
			my->UpdateBoundingBOx();
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
	if (key == OF_KEY_F2) {
	   cout << "Save file!" << endl;
	   my->storeModelMatrix();
	}
	if (key == OF_KEY_F3) {
		cout << "Load file!" << endl;
		my->LoadModelMatrix();
	}
	ofstream outFile;

	
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
	if (button == 0) {
		if (MouseSelectObjID != -1) {
			auto Cobjp = &my->Meshes[MouseSelectObjID];
			Surface* Cobj = (Surface*)Cobjp;
			Cobj->ResetMaterial();
		}
		MouseSelectObjID = (MouseSelectObjID+1)% my->Meshes.size();
		auto Cobjp = &my->Meshes[MouseSelectObjID];
		Surface* Cobj = (Surface*)Cobjp;
		Cobj->updateMaterial(mVec3(1, 0.02, 1) * 0.01, mVec3(1, 0.02, 1)*1.2, mVec3(1, 1, 1) * 1.0);

	}
	if (button == 1) {
		

		auto Cobjp = &my->Meshes[MouseSelectObjID];
		Surface* Cobj = (Surface*)Cobjp;
		Cobj->ResetMaterial();
		MouseSelectObjID = -1;
	   

	}
	if (button == 2) {

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
