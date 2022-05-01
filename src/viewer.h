#pragma once
#include <string>
#include "platforms/framebuffer.h"
extern "C" {
#include "platforms/platform.h"
}
#include <functional>
#include <memory>
#include "obj.hpp"

class Viewer
{

public:
    // C++-style functions
//
// Returns **true** if action should be cancelled.
//std::function<bool(Viewer& viewer)> callback_init;
//std::function<bool(Viewer& viewer)> callback_pre_draw;
//std::function<bool(Viewer& viewer)> callback_post_draw;
//std::function<bool(Viewer& viewer, int w, int h)> callback_post_resize;
//std::function<bool(Viewer& viewer, int button, int modifier)> callback_mouse_down;
//std::function<bool(Viewer& viewer, int button, int modifier)> callback_mouse_up;
//std::function<bool(Viewer& viewer, int mouse_x, int mouse_y)> callback_mouse_move;
//std::function<bool(Viewer& viewer, float delta_y)> callback_mouse_scroll;
//std::function<bool(Viewer& viewer, unsigned int key, int modifiers)> callback_key_pressed;
//std::function<bool(Viewer& viewer, unsigned int key, int modifiers)> callback_key_down;
//std::function<bool(Viewer& viewer, unsigned int key, int modifiers)> callback_key_up;
//std::function<bool(Viewer& viewer, unsigned int key, int modifiers)> callback_key_repeat;
//std::function<bool(Viewer& viewer)> InitAction;
//std::function<void(void)>DrawAction;   //RenderAframe

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

    void launch_init(
               const char* title, int width = 0, int height = 0, bool resizable = true, bool fullscreen = false, bool maximize = false );

    void launch_rendering(bool loop = true);
    void launch_shut();


    [[maybe_unused]] void init_plugins();
    [[maybe_unused]] void shutdown_plugins();

    Viewer();
    ~Viewer();

    // Callbacks
    bool key_pressed(unsigned int key, int modifier);
    bool key_down(int key, int modifier);
    bool key_up(int key, int modifier);
    bool key_repeat(int key, int modifier);
    bool mouse_down(MouseButton button, int modifier);
    bool mouse_up(MouseButton button, int modifier);
    bool mouse_move(int mouse_x, int mouse_y);
    bool mouse_scroll(float delta_y);


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
    int current_mouse_x = 0;
    int current_mouse_y = 0;
    int down_mouse_x = 0;
    int down_mouse_y = 0;
    float down_mouse_z = 0;
    //    Eigen::Vector3f down_translation;
    bool down;
    bool hack_never_moved;

    // Keep track of the global position of the scroll wheel
    float scroll_position;

    // Pointers to per-callback data
    /*void* callback_init_data;
    void* callback_pre_draw_data;
    void* callback_post_draw_data;
    void* callback_mouse_down_data;
    void* callback_mouse_up_data;
    void* callback_mouse_move_data;
    void* callback_mouse_scroll_data;
    void* callback_key_pressed_data;
    void* callback_key_down_data;
    void* callback_key_up_data;
    void* callback_key_repeat_data;*/
private:
    bool shouldcloseWindow = false;


    std::unique_ptr< framebuffer_t> m_framebuffer;
    std::unique_ptr< GraphicsPipeline> m_g_pip;
    std::unique_ptr< RenderObject> m_object;
    window* windowHandle;




};
