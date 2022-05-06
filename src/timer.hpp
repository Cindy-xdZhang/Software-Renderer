#pragma   once  
#ifndef TIMER
#define TIMER
#include <iostream>
#include <windows.h>   
class MyTimer
{
public:
	inline MyTimer(void);
	inline ~MyTimer(void);

private:
	LARGE_INTEGER start_time{};

	LARGE_INTEGER end_time{};

	LARGE_INTEGER CPUfrenquency = { };

public:
	double interval = 0;

public:
	inline void begin();
	inline void end();
};

inline MyTimer::MyTimer(void)
{
	QueryPerformanceFrequency(&CPUfrenquency);
}
inline MyTimer::~MyTimer(void)
{
}
inline void MyTimer::begin()
{
	interval = 0;
	QueryPerformanceCounter(&start_time);
}

inline void MyTimer::end()
{
	QueryPerformanceCounter(&end_time);

	interval = ((double)end_time.QuadPart - (double)start_time.QuadPart) / (double)CPUfrenquency.QuadPart;
}
#endif