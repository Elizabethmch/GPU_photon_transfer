#include "MyCudaToolkit.h"
#include <cuda_runtime.h>
#include <iostream>

#ifdef _WIN32
int gettimeofday(struct timeval* tp, void* tzp)
{
    time_t clock;
    struct tm tm;
    SYSTEMTIME wtm;
    GetLocalTime(&wtm);
    tm.tm_year = wtm.wYear - 1900;
    tm.tm_mon = wtm.wMonth - 1;
    tm.tm_mday = wtm.wDay;
    tm.tm_hour = wtm.wHour;
    tm.tm_min = wtm.wMinute;
    tm.tm_sec = wtm.wSecond;
    tm.tm_isdst = -1;
    clock = mktime(&tm);
    tp->tv_sec = clock;
    tp->tv_usec = wtm.wMilliseconds * 1000;
    return (0);
}
#endif

double cpuSecond() {
    LARGE_INTEGER t, tc;
    QueryPerformanceFrequency(&tc);
    QueryPerformanceCounter(&t);
    return ((double)t.QuadPart / (double)tc.QuadPart);

}


//void initialData_int(int* ip, int size) {
//    for (int i = 0; i < size; i++) {
//        ip[i] = 1;
//    }
//}


void initDevice(int devNum)
{
    int dev = devNum;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp, dev));
    printf("Using device %d: %s\n", dev, deviceProp.name);
    CHECK(cudaSetDevice(dev));

}
