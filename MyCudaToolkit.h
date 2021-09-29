#ifndef MYCUDATOOLKIT_H
#define MYCUDATOOLKIT_H

#define CHECK(call)\
{\
  const cudaError_t error=call;\
  if(error!=cudaSuccess)\
  {\
      printf("ERROR: %s:%d,",__FILE__,__LINE__);\
      printf("code:%d,reason:%s\n",error,cudaGetErrorString(error));\
      exit(1);\
  }\
}


#include <time.h>
#ifdef _WIN32
#	include <windows.h>
#else
#	include <sys/time.h>
#endif
#ifdef _WIN32
int gettimeofday(struct timeval* tp, void* tzp);
//{
//    time_t clock;
//    struct tm tm;
//    SYSTEMTIME wtm;
//    GetLocalTime(&wtm);
//    tm.tm_year = wtm.wYear - 1900;
//    tm.tm_mon = wtm.wMonth - 1;
//    tm.tm_mday = wtm.wDay;
//    tm.tm_hour = wtm.wHour;
//    tm.tm_min = wtm.wMinute;
//    tm.tm_sec = wtm.wSecond;
//    tm.tm_isdst = -1;
//    clock = mktime(&tm);
//    tp->tv_sec = clock;
//    tp->tv_usec = wtm.wMilliseconds * 1000;
//    return (0);
//}
#endif
//double cpuSecond()
//{
//  struct timeval tp;
//  gettimeofday(&tp,NULL);
//  return((double)tp.tv_sec+(double)tp.tv_usec*1e-6);
//
//}

double cpuSecond(); /*{
    LARGE_INTEGER t, tc;
    QueryPerformanceFrequency(&tc);
    QueryPerformanceCounter(&t);
    return ((double)t.QuadPart / (double)tc.QuadPart);

}*/


//void initialData_int(int* ip, int size) {
//    for (int i = 0; i < size; i++) {
//        ip[i] = 1;
//    }
//}


void initDevice(int devNum);
//{
//    int dev = devNum;
//    cudaDeviceProp deviceProp;
//    CHECK(cudaGetDeviceProperties(&deviceProp, dev));
//    printf("Using device %d: %s\n", dev, deviceProp.name);
//    CHECK(cudaSetDevice(dev));
//
//}


#endif  //MYCUDATOOLKIT_H
