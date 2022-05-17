//#include <chrono>
#include "viewer.h"
#include <filesystem>
#include <assert.h>
#include <string>

using namespace std;
namespace fs = std::filesystem;

static const std::string AssetsNames[] = {
    "KAUST_Beacon_NormalTexture.obj",
    "crab.obj",
	"crab.bmp",
     "Dream.bmp",
    "brick1.bmp",
	 "sky.bmp"
};
static MyTimer timer;

static bool mouse_dragging = false;
static int control_ob_id= 0;
extern mVec3f InitLightPos;

static Cube lightbox;


Viewer::Viewer() {
    m_framebuffer = std::make_unique<framebuffer_t>();
    m_objreader = std::make_unique<OBjReader>();
    m_g_pip = nullptr;
	m_mutable_objects = make_unique<std::vector<RenderableObject>>();
	m_FixObjects = make_unique<std::vector<RenderableObject>>();
}

Viewer::~Viewer() {

}

int  Viewer::initAssetsPath () {

    fs::path Curdir = fs::absolute(fs::current_path());
    fs::path Curdir_p = fs::absolute(Curdir.parent_path());
    fs::path Curdir_pp = fs::absolute(Curdir_p.parent_path());
    fs::path Curdir_ppp = fs::absolute(Curdir_pp.parent_path());

    fs::path assetdir0=Curdir / "assets";
    fs::path assetdir1 = Curdir_p/ "assets";
    fs::path assetdir2 = Curdir_pp/ "assets";
    fs::path assetdir3 = Curdir_ppp / "assets";
    if (fs::exists(assetdir0)) {
        m_projectRootDir = Curdir;
        return EXIT_SUCCESS;
    }
	else if (fs::exists(assetdir1)) {
        m_projectRootDir = Curdir_p;
		return EXIT_SUCCESS;
    }
	else if (fs::exists(assetdir2)) {
        m_projectRootDir = Curdir_pp;
        return EXIT_SUCCESS;
	}
	else if (fs::exists(assetdir3)) {
        m_projectRootDir = Curdir_ppp;
		return EXIT_SUCCESS;
	}
    else {
        return EXIT_FAILURE;
    }

}



void Viewer::init_callbacks( ) {
	static bool shift_is_pressed = false;
    if (windowHandle!=nullptr)
    {

		auto keypressed_lambda = [this](window_t* window, unsigned char key, int pressed)->void {
			if (pressed == 1)//for key pressed
			{
				if (key == unsigned char(16)) {
					shift_is_pressed = true;
				}
				if (key == unsigned char(27)) {
					//press esc to close
					windowHandle->should_close=true;
				}
				if (key == 'D') {
					if (shift_is_pressed) m_mutable_objects->at(control_ob_id).movePos(mVec3f(1, 0, 0));
					else m_g_pip->_Camera.UpdatePos(mVec3f(1, 0, 0));
				}
				if (key == 'A') {
					if (shift_is_pressed) m_mutable_objects->at(control_ob_id).movePos(mVec3f(-1, 0, 0));
					else m_g_pip->_Camera.UpdatePos(mVec3f(-1, 0, 0));
				}
				if ( key == 'S') {
					if (shift_is_pressed) m_mutable_objects->at(control_ob_id).movePos(mVec3f(0, -1, 0));
					else m_g_pip->_Camera.UpdatePos(mVec3f(0, -1, 0));
				}
				if ( key == 'W') {
					if (shift_is_pressed) m_mutable_objects->at(control_ob_id).movePos(mVec3f(0, 1, 0));
					else m_g_pip->_Camera.UpdatePos(mVec3f(0, 1, 0));
				}
				if (key == 'Q') {
					if (shift_is_pressed) m_mutable_objects->at(control_ob_id).movePos(mVec3f(0, 0, -1));
					else m_g_pip->_Camera.UpdatePos(mVec3f(0, 0, -1));
				}
				if ( key == 'E') {//go outside 
					if (shift_is_pressed) m_mutable_objects->at(control_ob_id).movePos(mVec3f(0, 0, 1));
					else m_g_pip->_Camera.UpdatePos(mVec3f(0, 0, 1));
				}
				if ( key == 'R') {//go into 
					if (shift_is_pressed) m_mutable_objects->at(control_ob_id).setModelMatrix(eye(4));
					else 	control_ob_id = (control_ob_id + 1) % (m_mutable_objects->size());

				}	
				if (key == 'C'&& shift_is_pressed) {//go into 
					m_g_pip->Clip = !m_g_pip->Clip;
					if (m_g_pip->Clip) std::cout << "clipping is on .";
					else std::cout << "clipping is off .";
				}
				

				if (key == 'y' || key == 'Y') {
					m_g_pip->_Camera.UpdateYawAngle(1.0);
				}
				if (key == 'u' || key == 'U') {
					m_g_pip->_Camera.UpdateYawAngle(-1.0);
				}
				if (key == 'p' || key == 'P') {//go outside 
					m_g_pip->_Camera.UpdatePitchAngle(1.0);
				}
				if (key == 'O' || key == 'o') {//go outside 
					m_g_pip->_Camera.UpdatePitchAngle(-1.0);
				}

				if (key == 'K' || key == 'k') {//go into 
					if (pressed == 1)//only for key pressed
						m_mutable_objects->at(control_ob_id).texture_scaling *= 0.5f;
				}
				if (key == 'b' || key == 'B') {//go into 
					if (pressed == 1)//only for key pressed
						m_g_pip->UsePhongShading = !m_g_pip->UsePhongShading;
				}
				if (key == 'l' || key == 'L') {//go into 
					if (pressed == 1)//only for key pressed
						m_mutable_objects->at(control_ob_id).texture_scaling *= 2;
				}
				//next texture mode
				if (key == 'n' || key == 'N') {//go into 
					int activeTexCount = m_g_pip->getActiveTextureNumber();
					m_mutable_objects->at(control_ob_id).texturechannelID = (m_mutable_objects->at(control_ob_id).texturechannelID + 1) % (activeTexCount + 2);
					std::cout << "Texture channel=" << m_mutable_objects->at(control_ob_id).texturechannelID << ". Active texture channel= " << activeTexCount << "\n";
				}


			}


			if (pressed == 0) {
				if (key == unsigned char(16)) {
					shift_is_pressed = false;
				}
			}

		};


		auto scroll_lambda = [this](window_t* window, float scrollY)->void {
            m_g_pip->_Camera.UpdateFov(-scrollY);
		};

		auto mouseDraggedLambda= [this](window_t* window, button_t button, int pressed)->void {
			static Matrix4 LastRot = eye(4);
			static Matrix4 ThisRot = eye(4);
			static int lastX = 0, lastY = 0;
			static mVec3 Op1, Op2;
	
			if (pressed==1)
			{
				if (button == BUTTON_R) {
					/*ThisRot = eye(4);
					m_objects->at(control_ob_id).setModelMatrix(ThisRot);*/
				}
				else if (button == BUTTON_L)
				{
					float x, y;
					input_query_cursor(windowHandle, &x, &y);
					lastX = x;
					lastY = y;
					Op1 = mArcControl.GetArcBallPositionVector(lastX, lastY);
					mouse_dragging = true;
					LastRot = ThisRot;
				}
			}
			else if (pressed==3&& button == BUTTON_L) {//drag
					if (mouse_dragging) {
						float x, y;
						input_query_cursor(windowHandle, &x, &y);
						//cout << "mouseDragged: y=" << y - lastY << " x=" << x - lastX << " button=" << button << endl;
						Op2 = mArcControl.GetArcBallPositionVector(x, y);
						ThisRot = mArcControl.GetArcBallrotateMatrix(Op1, Op2);
						Op1 = Op2;
						ThisRot = ThisRot * LastRot;
						auto preModelMat = m_mutable_objects->at(control_ob_id).getModelMatrix();
						m_mutable_objects->at(control_ob_id).setModelMatrix(preModelMat* ThisRot);
						
					}
			}
			else
			{
				LastRot = eye(4);
				ThisRot = eye(4);
				mouse_dragging = false;
			}
			
			
		};

		auto WindowResizeLambda = [this](window_t* window, const RECT& NEW_REC	) {

			m_g_pip->window_height = NEW_REC.bottom - NEW_REC.top + 1;
			m_g_pip->window_width= NEW_REC.right- NEW_REC.left+ 1;
			m_framebuffer->resize_buffer(m_g_pip->window_width, m_g_pip->window_height);
			mArcControl.updateWH(m_g_pip->window_width, m_g_pip->window_height);

		};


		//init callbacks
        windowHandle->callbacks.key_callback =lambda2FuncPtr<0>(keypressed_lambda);
        windowHandle->callbacks.scroll_callback = lambda2FuncPtr<1>(scroll_lambda);
        windowHandle->callbacks.button_callback = lambda2FuncPtr<2>(mouseDraggedLambda);
		windowHandle->callbacks.window_resize_callback= lambda2FuncPtr<3>(WindowResizeLambda);

    }

}


void Viewer::init_demo_scene() {
	//init_land by marching cube
	MarchingCubesDrawer tmp;
	tmp.randomGenrateLandmass({ 256,6,256 }, 0.25);

	RenderableObject land;
	land.BuildLikGLBegin(*tmp.getVertices(), *tmp.getNormals(), MeshBuildConvention::LIKE_GL_TRIANGLE);
	land.setModelMatrix(translateMatrix({240,0,0}) * scaleMatrix(480.0f,120.0f,240.0f));
	land.updateMaterial({ mVec3f(1.39,0.69,0.19) * 0.001, mVec3f(139,69,19)*0.002, mVec3f(0.0f) ,8 });
	land.fixit();
	Cube origin;
	land.setModelMatrix(scaleMatrix(0.25f));
	origin.fixit();
	Cube skybox;
	skybox.setModelMatrix(scaleMatrix(620.0f));
	skybox.texturechannelID = 3;

	skybox.fixit();
	m_FixObjects->push_back(std::move(skybox));
	m_FixObjects->push_back(std::move(origin));
	m_FixObjects->push_back(std::move(land));

	lightbox.updateMaterial({ mVec3f(2,1,1) * 0.005, mVec3f(2,1,1) * 0.0001, mVec3f(0) ,8 });

	fs::path  Path_obj_beacon = m_projectRootDir / "assets" / AssetsNames[0];
	fs::path  Path_obj = m_projectRootDir / "assets" / AssetsNames[1];
	fs::path  Path_texture = m_projectRootDir / "assets" / AssetsNames[2];
	fs::path  Path_texture_sky = m_projectRootDir / "assets" / AssetsNames[5];
	m_objreader->readObjFile(Path_obj.string());
	auto obj_crab = RenderableObject(m_objreader->TrianglesIdx, m_objreader->Vertexes, m_objreader->VertexesNormal, m_objreader->VertexesTexture);
	m_objreader->clear();
	m_objreader->readObjFile(Path_obj_beacon.string());
	auto obj_beacon = RenderableObject(m_objreader->TrianglesIdx, m_objreader->Vertexes, m_objreader->VertexesNormal, m_objreader->VertexesTexture);
	obj_beacon.updateMaterial({ mVec3f(0,0.69,0.99) * 0.005, mVec3f(0,69,99) * 2, mVec3f(0,69,99) ,8 });
	obj_crab.updateMaterial({ mVec3f(1,0.69,0.99) * 0.005, mVec3f(1,0.5,0.99) * 2, mVec3f(0.1,0.2,0.2) ,8 });
	obj_beacon.setModelMatrix(translateMatrix({ -55,-45,-395 }));
	/*m_objects->push_back(obj);
	obj.ModleMatrix *= translateMatrix(mVec3f(50, -0, -0));*/

	m_mutable_objects->push_back(std::move(obj_crab));
	m_mutable_objects->push_back(std::move(obj_beacon));

	m_g_pip->LoadTexture(Path_texture);
	m_g_pip->LoadTexture(Path_texture_sky);

}







void Viewer::launch_init(const char* title, int width /*= 0*/, int height /*= 0*/, bool resizable /*= true*/ , bool fullscreen /*= false*/ )
{

	
    
    try
    {
		m_framebuffer->resize_buffer(width, height);
		platform_initialize();
		windowHandle = window_create(title, width, height);
        init_callbacks();
		if (initAssetsPath() != EXIT_SUCCESS) {
			throw("Fail: can't find assets folder!");
		}

		m_g_pip = std::make_unique<GraphicsPipeline>(height, width);
		mArcControl = { DEFAULT_WINDOW_WIDTH,DEFAULT_WINDOW_HEIGHT };
		init_demo_scene();
    }
    catch (std::string errorInfo)
    {

        std::cout << "Exception Info: " << errorInfo;
		std::cout << "root Dir:" << m_projectRootDir.string() << "\n";
        exit(1);
    }
    catch (...) {
        exit(1);
    }
    
}


// string to wstring
static auto string2wchars (std::string str, std::wstring& szDst) {
	std::string temp = str;
	int len = MultiByteToWideChar(CP_ACP, 0, (LPCSTR)temp.c_str(), -1, NULL, 0);
	wchar_t* wszUtf8 = new wchar_t[len + 1];
	memset(wszUtf8, 0, len * 2 + 2);
	MultiByteToWideChar(CP_ACP, 0, (LPCSTR)temp.c_str(), -1, (LPWSTR)wszUtf8, len);
	szDst = wszUtf8;
	std::wstring r = wszUtf8;
	delete[] wszUtf8;
	return;
};


void Viewer::pre_draw() {
	static int a = 0;
	float angle = Radians(a);
	InitLightPos = {  cos(angle) ,InitLightPos.y ,sin(angle)+5 };
	lightbox.setModelMatrix(translateMatrix(InitLightPos) * scaleMatrix(0.2f) * rotateMatrix({ 0.45,0.45,0 }, 45));
	a += 1;
	a = a % 360;
    m_framebuffer->framebuffer_fast_clear();
	m_FixObjects->push_back(lightbox);
}

void Viewer::post_draw() {
	m_FixObjects->pop_back();
}

void Viewer::draw(bool first) {
    pre_draw();
	
    m_g_pip->Render(*m_FixObjects,*m_mutable_objects, m_framebuffer.get());
    
    window_draw_buffer(windowHandle, m_framebuffer.get());
    post_draw();
}


template< bool VERBOSE>
void Viewer::launch_rendering()
{

	
	// Rendering loop
	bool first_frame = true;
	int frame_counter = 0;
	while (!windowHandle->should_close)
	{
		timer.begin();
		draw(first_frame);
		first_frame = false;
		input_poll_events();
		timer.end();

		double fps = 1 / timer.interval;
		std::string output;
		if (m_g_pip->UsePhongShading)
			output = "PhongShading FPS= " + to_string(fps) + " .";
		else
			output = "GauroudShading FPS= " + to_string(fps) + " .";
		std::wstring wszDest;
		string2wchars(output, wszDest);
		SetWindowTextW(windowHandle->handle, wszDest.c_str());
		/*if constexpr (VERBOSE) {
		}*/
	}
}

template void Viewer::launch_rendering<true>();
template void Viewer::launch_rendering<false>();


void Viewer::launch_shut()
{
    platform_terminate();
	window_destroy(windowHandle);
    // core().shut(); // Doesn't do anything
    //shutdown_plugins();
    //glfwDestroyWindow(window);
    //glfwTerminate();
}
