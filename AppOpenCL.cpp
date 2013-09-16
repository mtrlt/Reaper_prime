#include "Global.h"
#include "AppOpenCL.h"
#include "Util.h"

#include <algorithm>

#ifndef CPU_MINING_ONLY
vector<_clState> GPUstates;
#endif

extern pthread_mutex_t current_work_mutex;
extern Work current_work;

#include <ctime>

extern ullint shares_hwinvalid;
extern ullint shares_hwvalid;

#include <deque>
using std::deque;

string VectorToHexString(vector<uchar> vec);

#ifndef CPU_MINING_ONLY

_clState clState;


cl_platform_id OpenCL::GetPlatform(uint id)
{
	cl_uint numPlatforms;
	cl_platform_id platform = NULL;

	status = clGetPlatformIDs(0, NULL, &numPlatforms);
	if(status != CL_SUCCESS)
		throw string("Error getting OpenCL platforms");

	if(numPlatforms > 0)
	{   
		cl_platform_id* platforms = new cl_platform_id[numPlatforms];

		status = clGetPlatformIDs(numPlatforms, platforms, NULL);
		if(status != CL_SUCCESS)
			throw string("Error getting OpenCL platform IDs");

		unsigned int i;
		cout << "List of platforms:" << endl;
		for(i=0; i < numPlatforms; ++i)
		{   
			char pbuff[100];

			status = clGetPlatformInfo( platforms[i], CL_PLATFORM_NAME, sizeof(pbuff), pbuff, NULL);
			if(status != CL_SUCCESS)
			{   
				delete [] platforms;
				throw string("Error getting OpenCL platform info");
			}

			cout << "\t" << i << "\t" << pbuff << endl;
			if (globalconfs.platform == id)
			{
				platform = platforms[id];
			}
		}   
		delete [] platforms;
	}
	else
	{
		throw string("No OpenCL platforms found");
	}
	if (platform == NULL)
		throw string("Chosen platform number does not exist");
	cout << "Using platform number " << globalconfs.platform << endl;

	return platform;
}

vector<cl_device_id> OpenCL::GetDeviceIDs(cl_platform_id platform, cl_uint devicetype)
{
	cl_uint numDevices;
	status = clGetDeviceIDs(platform, devicetype, 0, NULL, &numDevices);
	if(status != CL_SUCCESS)
	{
		throw string("Error getting OpenCL device IDs");
	}

	if (numDevices == 0)
		throw string("No OpenCL devices found");
	
	vector<cl_device_id> devices;
	cl_device_id* devicearray = new cl_device_id[numDevices];

	status = clGetDeviceIDs(platform, devicetype, numDevices, devicearray, NULL);
	if(status != CL_SUCCESS)
		throw string("Error getting OpenCL device ID list");

	for(uint i=0; i<numDevices; ++i)
		devices.push_back(devicearray[i]);
		
	delete [] devicearray;

	return devices;
}

string OpenCL::GetSourceFilename(uint device_id)
{
	string sourcefilename = config.GetCombiValue<string>("device", device_id, "kernel");
	if (sourcefilename == "")
		sourcefilename = config.GetValue<string>("kernel");
	sourcefilename = globalconfs.coin.protocol + "-" + sourcefilename;
	return sourcefilename;
}

string OpenCL::GetSource(uint device_id)
{
	string source;
	string sourcefilename;
	{
		sourcefilename = GetSourceFilename(device_id);
		FILE* filu = fopen(sourcefilename.c_str(), "rb");
		if (filu == NULL)
		{
			throw string("Couldn't find kernel file ") + sourcefilename;
		}
		fseek(filu, 0, SEEK_END);
		uint size = ftell(filu);
		fseek(filu, 0, SEEK_SET);
		size_t readsize = 0;
		for(uint i=0; i<size; ++i)
		{
			char c;
			readsize += fread(&c, 1, 1, filu);
			source.push_back(c);
		}
		if (readsize != size)
		{
			cout << "Read error while reading kernel source " << sourcefilename << endl;
		}
	}
	return source;
}

string OpenCL::GetBinaryFilename(cl_device_id device, uint device_id, map<string,string> defines)
{
	char pbuff[512];
	status = clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(pbuff), pbuff, NULL);
	if(status != CL_SUCCESS)
		throw string("Error getting OpenCL device info");
	
	string filebinaryname;
	for(char*p = &pbuff[0]; *p != 0; ++p)
	{
		//get rid of unwanted characters in filenames
		if (*p >= 33 && *p < 127 && *p != '\\' && *p != ':' && *p != '/' && *p != '*' && *p != '<' && *p != '>' && *p != '"' && *p != '?' && *p != '|')
			filebinaryname += *p;
	}
	for(map<string,string>::iterator it=defines.begin(); it!=defines.end(); ++it)
		filebinaryname += string("-") + it->first.substr(0,2) + it->second;
	string sourcefilename = GetSourceFilename(device_id);
	filebinaryname = sourcefilename.substr(0,sourcefilename.size()-3) + REAPER_VERSION + "." + filebinaryname + ".bin";
	return filebinaryname;
}

vector<uchar> OpenCL::GetBinary(string binaryfilename)
{
	vector<uchar> binaryfile;
	FILE* filu = fopen(binaryfilename.c_str(), "rb");
	if (filu != NULL)
	{
		fseek(filu, 0, SEEK_END);
		uint size = ftell(filu);
		fseek(filu, 0, SEEK_SET);
		if (size > 0)
		{
			binaryfile.assign(size,0);
			size_t readsize = fread(&binaryfile[0], size, 1, filu);
			if (readsize != 1)
				cout << "Read error while reading binary" << endl;
		}
		fclose(filu);
	}
	return binaryfile;
}

void OpenCL::Compile(_clState& GPUstate)
{
	string source = GetSource(GPUstate.device_id);
	vector<size_t> sourcesizes;
	sourcesizes.push_back(source.length());

	const char* see = source.c_str();
	cout << "Compiling kernel... this could take up to 2 minutes." << endl;
	char pbuff[512];
	status = clGetDeviceInfo(devices[GPUstate.device_id], CL_DEVICE_EXTENSIONS, sizeof(pbuff), pbuff, NULL);
	vector<string> extensions = Explode(string(pbuff),' ');

	string compile_options;
	for(map<string,string>::iterator it = GPUstate.defines.begin(); it!=GPUstate.defines.end(); ++it)
	{
		compile_options += " -D " + it->first + (!it->second.empty()?"="+it->second:"");
	}
	GPUstate.program = clCreateProgramWithSource(clState.context, 1, (const char **)&see, &sourcesizes[0], &status);
	if(status != CL_SUCCESS) 
		throw string("Error creating OpenCL program from source");

	status = clBuildProgram(GPUstate.program, 1, &devices[GPUstate.device_id], compile_options.c_str(), NULL, NULL);
	if(status != CL_SUCCESS) 
	{   
		size_t logSize;
		status = clGetProgramBuildInfo(GPUstate.program, devices[GPUstate.device_id], CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);

		char* log = new char[logSize];
		status = clGetProgramBuildInfo(GPUstate.program, devices[GPUstate.device_id], CL_PROGRAM_BUILD_LOG, logSize, log, NULL);
		cout << log << endl;
		delete [] log;
		throw string("Error building OpenCL program");
	}
}

void OpenCL::SaveBinary(_clState& GPUstate, string binaryfilename)
{
	uint device_amount;
	clGetProgramInfo(GPUstate.program, CL_PROGRAM_NUM_DEVICES, sizeof(uint), &device_amount, NULL);

	size_t* binarysizes = new size_t[device_amount];
	uchar** binaries = new uchar*[device_amount];
	for(uint curr_binary = 0; curr_binary<device_amount; ++curr_binary)
	{
		clGetProgramInfo(GPUstate.program, CL_PROGRAM_BINARY_SIZES, device_amount*sizeof(size_t), binarysizes, NULL);
		binaries[curr_binary] = new uchar[binarysizes[curr_binary]];
	}
	clGetProgramInfo(GPUstate.program, CL_PROGRAM_BINARIES, sizeof(uchar*)*device_amount, binaries, NULL);
	for(uint binary_id = 0; binary_id < device_amount; ++binary_id)
	{
		if (binarysizes[binary_id] == 0)
			continue;

		cout << "Saving binary, size: " << binarysizes[binary_id] << " bytes" << endl;
		FILE* filu = fopen(binaryfilename.c_str(), "wb");
		fwrite(binaries[binary_id], binarysizes[binary_id], 1, filu);
		fclose(filu);
	}

	cout << "Program built from source." << endl;
	delete [] binarysizes;
	for(uint binary_id=0; binary_id < device_amount; ++binary_id)
		delete [] binaries[binary_id];
	delete [] binaries;
}

void OpenCL::LoadBinary(_clState& GPUstate, vector<uchar>& binaryfile)
{
	cl_int binary_status, errorcode_ret;
	size_t filebinarysize = binaryfile.size();
	{
		const uchar* pointer1 = &binaryfile[0];
		const uchar** pointer2 = &pointer1;
		GPUstate.program = clCreateProgramWithBinary(clState.context, 1, &GPUstate.device, &filebinarysize, pointer2, &binary_status, &errorcode_ret);
	}
	if (binary_status != CL_SUCCESS)
		cout << "Binary status error code: " << binary_status << endl;
	if (errorcode_ret != CL_SUCCESS)
		cout << "Binary loading error code: " << errorcode_ret << endl;
	status = clBuildProgram(GPUstate.program, 1, &GPUstate.device, NULL, NULL, NULL);
	if (status != CL_SUCCESS)
		cout << "Error while building from binary: " << status << endl;

	cout << "Program built from saved binary." << endl;	
}

void OpenCL::WriteBuffer(uint device_num, string buffername, void* data, size_t data_length)
{
	_clState& GPUstate = GPUstates[device_num];
	if (GPUstate.buffers[buffername] == NULL)
		cout << "Buffer " << buffername << " not found on GPU #" << device_num << endl;
	cl_int status = clEnqueueWriteBuffer(GPUstate.commandQueue, GPUstate.buffers[buffername], true, 0, data_length, data, 0, NULL, NULL);
	if (globalconfs.coin.config.GetValue<bool>("opencldebug"))
		cout << "Write buffer " << buffername << ", " << data_length << " bytes. Status: " << status << endl;
}
void OpenCL::WriteBufferPattern(uint device_num, string buffername, size_t data_length, void* pattern, size_t pattern_length)
{
	_clState& GPUstate = GPUstates[device_num];
	if (GPUstate.buffers[buffername] == NULL)
		cout << "Buffer " << buffername << " not found on GPU #" << device_num << endl;
	cl_int status = clEnqueueFillBuffer(GPUstate.commandQueue, GPUstate.buffers[buffername], pattern, pattern_length, 0, data_length, 0, NULL, NULL);
	if (globalconfs.coin.config.GetValue<bool>("opencldebug"))
		cout << "Write buffer pattern " << buffername << ", " << pattern_length << " bytes. Status: " << status << endl;
}
void OpenCL::ReadBuffer(uint device_num, string buffername, void* data, size_t data_length)
{
	_clState& GPUstate = GPUstates[device_num];
	if (GPUstate.buffers[buffername] == NULL)
		cout << "Buffer " << buffername << " not found on GPU #" << device_num << endl;
	cl_int status = clEnqueueReadBuffer(GPUstate.commandQueue, GPUstate.buffers[buffername], true, 0, data_length, data, 0, NULL, NULL);
	if (globalconfs.coin.config.GetValue<bool>("opencldebug"))
		cout << "Read buffer " << buffername << ", " << data_length << " bytes. Status: " << status << endl;
}

void OpenCL::WriteBufferGlobal(string buffername, void* data, size_t data_length)
{
	for(uint i=0; i<GPUstates.size(); ++i)
		WriteBuffer(i,buffername,data,data_length);
}
void OpenCL::ReadBufferGlobal(string buffername, void* data, size_t data_length)
{
	for(uint i=0; i<GPUstates.size(); ++i)
		ReadBuffer(i,buffername,data,data_length);
}

void OpenCL::RunKernel(uint device_num, string name, size_t globalsize, vector<std::pair<string,string> > arguments)
{
	while(!inited);
	
	_clState* state = &GPUstates[device_num];
	
	cl_kernel kernel = state->kernels[name];
	if (kernel == NULL)
		cout << "Kernel " << name << " not found." << endl;
	
	for(uint i=0; i<arguments.size(); ++i)
	{
		if (arguments[i].first == "buffer")
		{
			if (state->buffers.find(arguments[i].second) == state->buffers.end())
				cout << "Buffer " << arguments[i].second << " not found." << endl;
			clSetKernelArg(kernel, i, sizeof(cl_mem), &state->buffers[arguments[i].second]);
		}
		else
		{
			clSetKernelArg(kernel, i, arguments[i].second.size(), &arguments[i].second[0]);
		}
	}
	
	size_t localsize = globalconfs.coin.local_worksize;
	size_t base=0;
	
	//TODO: figure out what this implies
	globalsize -= globalsize%localsize;
	
	cl_int status = clEnqueueNDRangeKernel(state->commandQueue, kernel, 1, &base, &globalsize, &localsize, 0, NULL, NULL);
	if (globalconfs.coin.config.GetValue<bool>("opencldebug"))
		cout << "Kernel " << name << " run, globalsize " << globalsize << ". Status: " << status << endl;
}

void OpenCL::RebuildKernels(_clState& GPUstate, map<string,string> defines)
{
	bool changed=false;
	for(map<string,string>::iterator it = defines.begin(); it != defines.end(); ++it)
	{
		if (GPUstate.defines[it->first] != it->second)
		{
			GPUstate.defines[it->first] = it->second;
			changed=true;
		}
	}
	if (!changed)
		return;
		
	string binaryfilename = GetBinaryFilename(GPUstate.device, GPUstate.device_id, GPUstate.defines);
	vector<uchar> binaryfile;
	if (globalconfs.save_binaries)
		binaryfile = GetBinary(binaryfilename);
	
	if (binaryfile.empty())
	{
		Compile(GPUstate);
		
		if (globalconfs.save_binaries)
			SaveBinary(GPUstate, binaryfilename);
	}
	else
		LoadBinary(GPUstate, binaryfile);

	vector<string>& kernelnames = GPUstate.kernelnames;
		
	for(uint i=0; i<kernelnames.size(); ++i)
	{
		if (GPUstate.kernels[kernelnames[i]] != NULL)
			clReleaseKernel(GPUstate.kernels[kernelnames[i]]);
		cout << "Building kernel " << kernelnames[i].c_str() << endl;
		cl_kernel kernel = clCreateKernel(GPUstate.program, kernelnames[i].c_str(), &status);
		if(status != CL_SUCCESS)
		{
			cout << "Kernel build not successful: " << status << endl;
			throw string("Error creating OpenCL kernel");
		}
		GPUstate.kernels[kernelnames[i]] = kernel;
	}
}

#endif

#include "Config.h"
extern Config config;
void OpenCL::Init(vector<string> kernelnames, map<string,cl_mem_settings> buffersettings, map<string,string> defines)
{
#ifdef CPU_MINING_ONLY
	if (globalconfs.coin.config.GetValue<bool>("use_gpu"))
	{
		cout << "This binary was built with CPU mining support only." << endl;
	}
#else
	if (globalconfs.coin.config.GetValue<bool>("use_gpu") == false)
	{
		cout << "GPU not used." << endl;
		return;
	}

	cl_uint devicetype = CL_DEVICE_TYPE_GPU;
	
	cl_platform_id platform = GetPlatform(globalconfs.platform);
	devices = GetDeviceIDs(platform,devicetype);
	uint numDevices = devices.size();

	cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0 };
	cl_int status = 0;
	clState.context = clCreateContextFromType(cps, devicetype, NULL, NULL, &status);
	if(status != CL_SUCCESS) 
		throw string("Error creating OpenCL context");

	cout << endl;
	if (globalconfs.devices.empty())
	{
		for(uint i=0; i<numDevices; ++i)
			globalconfs.devices.push_back(i);
	}
	{
		cout << "Using device" << (globalconfs.devices.size()==1?"":"s") << " ";
		for(uint i=0; i<globalconfs.devices.size(); ++i)
		{
			cout << globalconfs.devices[i];
			if (i+1 < globalconfs.devices.size())
			{
				cout << ", ";
			}
		}
		cout << endl;
	}
	for(uint h=0; h<globalconfs.devices.size(); ++h)
	{
		_clState GPUstate;
		uint device_id = globalconfs.devices[h];
		GPUstate.device_id = device_id;
		GPUstate.device = devices[device_id];
		GPUstate.kernelnames = kernelnames;
		cout << "OpenCL device " << device_id << "........" << endl;

		//default values for defines
		defines["WORKSIZE"] = ToString(globalconfs.coin.local_worksize);
		defines["NCHAINLENGTH"] = ToString(5);
		defines["LOCAL_MEM_USED"] = ToString(4096);
		
		RebuildKernels(GPUstate,defines);
		

		{
			for(map<string,cl_mem_settings>::iterator it = buffersettings.begin(); it!=buffersettings.end(); ++it)
			{
				cl_mem tmp = clCreateBuffer(clState.context, it->second.flags, it->second.size, NULL, &status);
				if (globalconfs.coin.config.GetValue<bool>("opencldebug"))
					cout << "Created buffer " << it->first << ". size " << it->second.size << ", flags " << it->second.flags << ". Status: " << status << endl;					
				GPUstate.buffers[it->first] = tmp;
			}
			
			GPUstate.commandQueue = clCreateCommandQueue(clState.context, devices[device_id], CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &status);
			if (status != CL_SUCCESS)
				GPUstate.commandQueue = clCreateCommandQueue(clState.context, devices[device_id], 0, &status);
			if(status != CL_SUCCESS)
				throw string("Error creating OpenCL command queue");

			if(status != CL_SUCCESS)
			{
				cout << status << endl;
				throw string("Error creating OpenCL buffer");
			}

			GPUstate.offset = 0x100000000ULL/numDevices*(device_id);

			pthread_mutex_init(&GPUstate.share_mutex,NULL);
			GPUstate.shares_available = false;
			pthread_mutex_init(&GPUstate.device_mutex,NULL);

			GPUstate.thread_id = device_id;
			GPUstate.in_use = false;
			GPUstates.push_back(GPUstate);
		}
	}

	if (GPUstates.empty())
	{
		cout << "No GPUs selected." << endl;
		return;
	}
	cout << "done" << endl;
#endif
}

void OpenCL::Quit()
{
	
}
