#pragma once

//framebuffers
 struct framebuffer {
	int width, height;
	unsigned char* color_buffer;
	float* depth_buffer;

	framebuffer():width(-1),height(-1),color_buffer(nullptr), depth_buffer(nullptr)
	{
		
	}
	framebuffer(int width, int height);

	~framebuffer();

	void resize_buffer(int width, int height);
	void framebuffer_clear_color( const float color[4]);
	void framebuffer_clear_depth(float depth);


	inline void setvalue(int r, int c, unsigned char color_rgba[4]) {

		int  dst_index = (r * width + c) * 4;
		unsigned char* dst_pixel = &color_buffer[dst_index];
		dst_pixel[0] = color_rgba[0];  /* red */
		dst_pixel[1] = color_rgba[1];  /* green */
		dst_pixel[2] = color_rgba[2];  /* blue */
		dst_pixel[3] = color_rgba[3];  /* blue */
	}

} ;


 void release_buffer(framebuffer* fb);

 typedef framebuffer framebuffer_t;
/* framebuffer management */

