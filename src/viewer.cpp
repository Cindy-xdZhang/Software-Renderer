//#include <chrono>
#include "viewer.h"
#include <filesystem>
#include <assert.h>

using namespace std;
namespace fs = std::filesystem;

static const std::string AssetsNames[] = {
    "KAUST_Beacon.obj",
    "KAUST_Beacon_NormalTexture.obj",
    "brick1.bmp"
};


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

void Viewer::launch_init(const char* title, int width /*= 0*/, int height /*= 0*/, bool resizable /*= true*/ , bool fullscreen /*= false*/ )
{
    try
    {
		m_framebuffer->resize_buffer(width, height);
		platform_initialize();
		if (initAssetsPath() != EXIT_SUCCESS) {
			throw("Fail: can't find assets folder!");
		}
		m_g_pip = std::make_unique<GraphicsPipeline>(height, width);
		windowHandle = window_create(title, width, height);
		fs::path  Path_obj = m_projectRootDir / "assets" / AssetsNames[1];
		fs::path  Path_texture = m_projectRootDir / "assets" / AssetsNames[2];

		m_objreader->read(Path_obj);
		auto obj = RenderableObject(m_objreader->TrianglesIdx, m_objreader->Vertexes, m_objreader->VertexesNormal, m_objreader->VertexesTexture);
		m_objects->push_back(obj);
		obj.ModleMatrix *= translateMatrix(mVec3f(50, -0, -0));

		m_objects->push_back(std::move(obj));
		
		m_g_pip->LoadTexture(Path_texture);

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


void Viewer::pre_draw() {
    static int a = 0;
    float angle = Radians(a);
    m_g_pip->InitLightPos = { 125 * cos(angle) ,55 ,125 * sin(angle) };
    a += 5;
    a = a % 360;
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


void Viewer::launch_rendering(bool loop)
{
    // Rendering loop
    bool first_frame = true;
    int frame_counter = 0;
    while (!shouldcloseWindow)
    {
        double tic = std::chrono::duration<double>(
            std::chrono::system_clock::now().time_since_epoch()).count();

        draw(first_frame);

        first_frame = false;
        //glfwSwapBuffers(window);
        
		// In microseconds
	/*	double toc = std::chrono::duration<double>(
			std::chrono::system_clock::now().time_since_epoch()).count();
		double duration = 1000000. * (toc - tic);
		double fps = 1 / duration;*/
        if (!loop)
        {
            return;
        }

    }
}



void Viewer::launch_shut()
{
    platform_terminate();
    // core().shut(); // Doesn't do anything
    //shutdown_plugins();
    //glfwDestroyWindow(window);
    //glfwTerminate();
}
