#ifndef MTRLTINOPENCL_H
#define MTRLTINOPENCL_H

#ifndef CPU_MINING_ONLY
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif
#endif

#include "pthread.h"


#ifndef CPU_MINING_ONLY
struct cl_mem_settings
{
	size_t size;
	cl_mem_flags flags;
	cl_mem_settings():size(0),flags(0){}
	cl_mem_settings(size_t size_, cl_mem_flags flags_): size(size_),flags(flags_){}
	~cl_mem_settings(){}
};

struct _clState
{
	cl_context context;
	map<string,cl_kernel> kernels;
	cl_command_queue commandQueue;
	cl_program program;
	map<string,cl_mem> buffers;
	map<string,string> defines;
	vector<string> kernelnames;
	cl_uint device_id;
	cl_device_id device;
	pthread_mutex_t device_mutex;
	
	uint thread_id;
	uint offset;
	
	bool in_use;

	pthread_t thread;

	bool shares_available;
	deque<Share> shares;
	pthread_mutex_t share_mutex;

	ullint hashes;
};

extern vector<_clState> GPUstates;
#endif

class OpenCL
{
private:
#ifndef CPU_MINING_ONLY
	cl_int status;
	vector<cl_device_id> devices;
#endif
	
public:
	OpenCL():inited(false){}
	~OpenCL(){}

	void Init(vector<string> kernelnames, map<string,cl_mem_settings> buffersettings, map<string,string> defines);
	void Quit();
	
	
#ifndef CPU_MINING_ONLY
	bool inited;

	cl_platform_id GetPlatform(uint id);
	vector<cl_device_id> GetDeviceIDs(cl_platform_id platform, cl_uint devicetype);
	string GetSourceFilename(uint device_id);
	string GetSource(uint device_id);
	string GetBinaryFilename(cl_device_id device, uint device_id, map<string,string> defines);
	vector<uchar> GetBinary(string filename);
	void Compile(_clState& GPUstate);
	void SaveBinary(_clState& GPUstate, string binaryfilename);
	void LoadBinary(_clState& GPUstate, vector<uchar>& binaryfile);
	void RebuildKernels(_clState& GPUstate, map<string,string> defines);

	void RunKernel(uint device_num, string name, size_t globalsize, vector<std::pair<string,string> > arguments);
	void WriteBuffer(uint device_num, string buffername, void* data, size_t data_length);
	void WriteBufferPattern(uint device_num, string buffername, size_t data_length, void* pattern, size_t pattern_length);
	void ReadBuffer(uint device_num, string buffername, void* data, size_t data_length);

	void WriteBufferGlobal(string buffername, void* data, size_t data_length);
	void ReadBufferGlobal(string buffername, void* data, size_t data_length);
#endif
};


#include "Util.h"
#include "App.h"


#endif
