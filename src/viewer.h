#pragma once
#include <string>
#include "platforms/framebuffer.h"
#include "platforms/platform.h"
#include <functional>
#include <memory>
#include "GraphicsPipeline.h"
#include "objparser.h"
class Viewer
{

public:

public:
    // UI Enumerations
    enum class MouseButton
    {
        Left, Middle, Right
    };
    enum class MouseMode
    {
        None, Rotation, Zoom, Pan, Translation
    };
	//    Eigen::Vector3f down_translation;

    void launch_init(
               const char* title, int width = 0, int height = 0, bool resizable = true, bool fullscreen = false );

    template<bool SET_MAX_FPS = true>
    void launch_rendering();
    void launch_shut();


    [[maybe_unused]] void init_plugins();
    [[maybe_unused]] void shutdown_plugins();

    Viewer();
    ~Viewer();

    // Draw everything
    void pre_draw();
    void draw(bool first);
    void post_draw();

    // resize
    void resize(int w, int h); // explicitly set window size
    void post_resize(int w, int h); // external resize due to user interaction

    ////////////////////////////////////////////////////////////////////////////
    /// Member variables
    ////////////////////////////////////////////////////////////////////////////
    // List of registered plugins
    //std::vector<ViewerPlugin*> plugins;

    // Temporary data stored when the mouse button is pressed
    MouseMode mouse_mode = Viewer::MouseMode::None;
    //Eigen::Quaternionf down_rotation;

   
private:

    ArcBallControler mArcControl;

    void init_callbacks();

    int  initAssetsPath();
    std::filesystem::path  m_projectRootDir;

    std::unique_ptr< framebuffer_t> m_framebuffer;
    std::unique_ptr< OBjReader> m_objreader;

    std::unique_ptr< GraphicsPipeline> m_g_pip;
    std::unique_ptr< std::vector<RenderableObject>> m_objects;


    window_t* windowHandle=nullptr;
    //std::unordered_map<int, keycallbackFunctionType> mBindKeycallbacks;



};

using keycallbackFunctionType = void (*)(window_t* window, unsigned char key, int pressed);
using scrollcallbackFunctionType = void (*)(window_t* window, unsigned char key, int pressed);


template<int CALLBACK_TYPE,typename F>
auto lambda2FuncPtr(F lambda) {
    static auto lambdaBank = lambda;
    if constexpr (CALLBACK_TYPE == 0) {
    return [](window_t* window, unsigned char key, int pressed) ->void{ lambdaBank(window,key,pressed); };

    }
    else if constexpr (CALLBACK_TYPE == 1) {
    return [](window_t* window, float scroll) ->void { lambdaBank(window, scroll); };

	}
	else if constexpr (CALLBACK_TYPE == 2) {

		return [](window_t* window, button_t button, int pressed) ->void { lambdaBank(window, button, pressed); };

	}
	else if constexpr (CALLBACK_TYPE == 3) {

		return [](window_t* window, const RECT& rec) ->void { lambdaBank(window, rec); };

	}
}

