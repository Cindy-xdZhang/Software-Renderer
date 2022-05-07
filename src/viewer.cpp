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
    "brick1.bmp"
};
static MyTimer timer;

static bool mouse_dragging = false;
static int control_ob_id= 0;
extern mVec3f InitLightPos;
extern float texture_scaling;

Viewer::Viewer() {
    m_framebuffer = std::make_unique<framebuffer_t>();
    m_objreader = std::make_unique<OBjReader>();
    m_g_pip = nullptr;
    m_objects = make_unique<std::vector<RenderableObject>>();
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
    if (windowHandle!=nullptr)
    {

		mArcControl = {DEFAULT_WINDOW_WIDTH,DEFAULT_WINDOW_HEIGHT};

		auto keypressed_lambda = [this](window_t* window, unsigned char key, int pressed)->void {
        
			if (key == 'd' || key == 'D') {
				m_g_pip->_Camera.UpdatePos(mVec3f(5, 0, 0));
			}
			if (key == 'a'|| key == 'A') {
				m_g_pip->_Camera.UpdatePos(mVec3f(-5, 0, 0));
			}
			if (key == 's' || key == 'S') {
				m_g_pip->_Camera.UpdatePos(mVec3f(0, -5, 0));
			}
			if (key == 'w' || key == 'W') {
				m_g_pip->_Camera.UpdatePos(mVec3f(0, 5, 0));
			}
			if (key == 'q' || key == 'Q') {
				InitLightPos.z -= 5;
			}
			if (key == 'e' || key == 'E') {//go outside 
				m_g_pip->_Camera.UpdatePos(mVec3f(0, 0, -5));
			}
			if (key == 'r' || key == 'R') {//go into 
				m_g_pip->_Camera.UpdatePos(mVec3f(0, 0, 5));
			}
			if (key == 'b' || key == 'B') {//go into 
				if (pressed == 1)//only for key pressed
				    m_g_pip->UsePhongShading = !m_g_pip->UsePhongShading;
			}
			if (key == 'l' || key == 'L') {//go into 
				if (pressed == 1)//only for key pressed
					texture_scaling *= 2;
			}
			if (key == 'K' || key == 'k') {//go into 
				if (pressed == 1)//only for key pressed
					texture_scaling *= 0.5f;
			}
			//next texture mode
			if (key == 'n' || key == 'N') {//go into 
                if (pressed == 1)//only for key pressed
                {
                    m_g_pip->TextureMode = (m_g_pip->TextureMode + 1) % 3;
                    std::cout << "Texture Mode=" << m_g_pip->TextureMode << " texture channel=default.\n";
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
					ThisRot = eye(4);
					m_objects->at(control_ob_id).ModleMatrix = ThisRot;
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
						cout << "mouseDragged: y=" << y - lastY << " x=" << x - lastX << " button=" << button << endl;
						Op2 = mArcControl.GetArcBallPositionVector(x, y);
						ThisRot = mArcControl.GetArcBallrotateMatrix(Op1, Op2);
						ThisRot = ThisRot * LastRot;
						m_objects->at(control_ob_id).ModleMatrix = ThisRot ;
					}
			}
			else
			{
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
	tmp.randomGenrateLandmass({ 32,32,5}, 0.1);

	RenderableObject land;
	land.BuildLikGLBegin(*tmp.getVertices(), *tmp.getNormals(), MeshBuildConvention::LIKE_GL_TRIANGLE);
	land.ModleMatrix = scaleMatrix(30.0f) * rotateMatrix({1,0,0}, 90);
	land.updateMaterial({ mVec3f(1.39,1.17,0) * 0.02, mVec3f(1.39,1.17,0) * 60, mVec3f(1.39,1.17,0) * 2,8 });



	fs::path  Path_obj = m_projectRootDir / "assets" / AssetsNames[1];
	fs::path  Path_texture = m_projectRootDir / "assets" / AssetsNames[2];

	m_objreader->readObjFile(Path_obj.string());
	auto obj = RenderableObject(m_objreader->TrianglesIdx, m_objreader->Vertexes, m_objreader->VertexesNormal, m_objreader->VertexesTexture);
	/*m_objects->push_back(obj);
	obj.ModleMatrix *= translateMatrix(mVec3f(50, -0, -0));*/

	m_objects->push_back(std::move(obj));

	m_g_pip->LoadTexture(Path_texture);

	m_objects->push_back(std::move(land));
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
	InitLightPos = { 125 * cos(angle) ,55 ,125 * sin(angle) };
	a += 5;
	a = a % 360;
    m_framebuffer->framebuffer_fast_clear();
}

void Viewer::post_draw() {

}

void Viewer::draw(bool first) {
    pre_draw();
    for (const auto& m_object: (* m_objects))
    {
        m_g_pip->Render(m_object, m_framebuffer.get());
    }

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
