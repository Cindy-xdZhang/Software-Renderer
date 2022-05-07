#pragma once
#include"framebuffer.h"
#include "../Macros.h"
extern "C" {
#include "image.h"
#include <windows.h>
}

typedef enum { KEY_A, KEY_D, KEY_S, KEY_W, KEY_Q,KEY_B, KEY_T, KEY_N, KEY_E, KEY_R, KEY_SPACE, KEY_NUM } keycode_t;
typedef enum { BUTTON_L, BUTTON_R, BUTTON_NUM } button_t;

typedef struct window window_t;

typedef struct {
	void (*key_callback)(window_t* window, unsigned char key, int pressed);
	void (*button_callback)(window_t* window, button_t button, int pressed);
	void (*scroll_callback)(window_t* window, float offset);
	void (*window_resize_callback)(window_t* window, const RECT& newRect);
} callbacks_t;

struct window {
	HWND handle;
	HDC memory_dc;
	image_t* surface;
	/* common data */
	int should_close;
	char keys[KEY_NUM];
	char buttons[BUTTON_NUM];
	callbacks_t callbacks;
	void* userdata;
};





extern "C" {


	/* platform initialization */
	void platform_initialize(void);
	void platform_terminate(void);

	/* window related functions */
	window_t* window_create(const char* title, int width, int height);
	void window_destroy(window_t* window);
	int window_should_close(window_t* window);
	void window_set_userdata(window_t* window, void* userdata);
	void* window_get_userdata(window_t* window);
	void window_draw_buffer(window_t* window, framebuffer_t* buffer);

	/* input related functions */
	void input_poll_events(void);
	int input_key_pressed(window_t* window, keycode_t key);
	int input_button_pressed(window_t* window, button_t button);
	void input_query_cursor(window_t* window, float* xpos, float* ypos);
	void input_set_callbacks(window_t* window, callbacks_t callbacks);

	/* misc platform functions */
	float platform_get_time(void);



}

