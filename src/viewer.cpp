#include <chrono>
#include "viewer.h"
#include "objparser.h"
using namespace std;
Viewer::Viewer() {
    m_framebuffer = std::make_unique<framebuffer_t>();
    m_g_pip = nullptr;
}
Viewer::~Viewer() {

}
int Viewer::launch_init(const char* title, int windowWidth, int windowHeight, bool resizable , bool fullscreen , bool maximize )
{
    m_framebuffer->resize_buffer(windowWidth,windowHeight);
    platform_initialize();

    OBjReader objFilereader;
    string Path = string("C:\\Users\\zhanx0o\\Documents\\sources\\Software-Renderer\\assets\\KAUST_Beacon_NormalTexture.obj");
    objFilereader.read(Path);
    m_object = std::make_unique < RenderObject>(objFilereader.TrianglesIdx, objFilereader.Vertexes, objFilereader.VertexesNormal, objFilereader.VertexesTexture);
    m_g_pip = std::make_unique<GraphicsPipeline>(320, 480);

    windowHandle =window_create(title, 320, 480);

    return 0;
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

    m_g_pip->Render(*m_object, m_framebuffer.get());
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
		double toc = std::chrono::duration<double>(
			std::chrono::system_clock::now().time_since_epoch()).count();
		double duration = 1000000. * (toc - tic);
        double fps = 1 / duration;
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
