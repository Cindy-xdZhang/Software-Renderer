#include "framebuffer.h"
#include<memory>
#include<assert.h>

#define float_to_uchar(x)  (unsigned char)(x*255)

/* framebuffer management */
framebuffer::framebuffer(int width, int height) :width(width), height(height)
{
	int color_buffer_size = width * height * 4;
	int depth_buffer_size = sizeof(float) * width * height;
	float default_color[4] = { 0, 0, 0, 1 };
	float default_depth = 1;
	assert(width > 0 && height > 0);

	framebuffer_t* framebuffer = this;
	framebuffer->width = width;
	framebuffer->height = height;
	framebuffer->color_buffer = (unsigned char*)malloc(color_buffer_size);
	framebuffer->depth_buffer = (float*)malloc(depth_buffer_size);

	framebuffer_clear_color(default_color);
	framebuffer_clear_depth(default_depth);

}

framebuffer::~framebuffer() {
	if (depth_buffer)
	{
		free(depth_buffer);
		depth_buffer = nullptr;
	}
	if (color_buffer)
	{
		free(color_buffer);
		color_buffer = nullptr;
	}
}



void release_buffer(framebuffer_t* fb) {
		free(fb->depth_buffer);
		fb->depth_buffer = nullptr;

		free(fb->color_buffer);
		fb->color_buffer = nullptr;
	
}

void framebuffer::resize_buffer(int width, int height) {
	int color_buffer_size = width * height * 4;
	int depth_buffer_size = sizeof(float) * width * height;
	float default_color[4] = { 0, 0, 0, 1 };
	float default_depth = 1;
	assert(width > 0 && height > 0);
	if (depth_buffer)
		free(depth_buffer);
	if (color_buffer)
		free(color_buffer);

	framebuffer_t* framebuffer = this;
	framebuffer->width = width;
	framebuffer->height = height;

	framebuffer->color_buffer = (unsigned char*)malloc(color_buffer_size);
	framebuffer->depth_buffer = (float*)malloc(depth_buffer_size);

	framebuffer_clear_color(default_color);
	framebuffer_clear_depth(default_depth);

}




void framebuffer::framebuffer_clear_color(const float color[4]) {
	framebuffer_t* framebuffer = this;
    int num_pixels = framebuffer->width * framebuffer->height;
    int i;
    for (i = 0; i < num_pixels; i++) {
        framebuffer->color_buffer[i * 4 + 0] = float_to_uchar(color[0]);
        framebuffer->color_buffer[i * 4 + 1] = float_to_uchar(color[1]);
        framebuffer->color_buffer[i * 4 + 2] = float_to_uchar(color[2]);
        framebuffer->color_buffer[i * 4 + 3] = float_to_uchar(color[3]);
    }
}

void framebuffer::framebuffer_clear_depth( float depth) {
	framebuffer_t* framebuffer = this;
    int num_pixels = framebuffer->width * framebuffer->height;
    int i;
    for (i = 0; i < num_pixels; i++) {
        framebuffer->depth_buffer[i] = depth;
    }
}

void framebuffer::framebuffer_fast_clear()
{
	size_t SizeBuff= static_cast<size_t>(this->width) * static_cast<size_t>(this->height);
	SizeBuff *= sizeof(float);
	memset(this->color_buffer, 0, SizeBuff);
	
}

