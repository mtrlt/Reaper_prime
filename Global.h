#ifndef GLOBAL_H
#define GLOBAL_H

#include "CMakeConf.h"

#define _CRT_SECURE_NO_WARNINGS

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

#include <vector>
#include <iostream>
#include <string>
#include <deque>
#include <utility>
using namespace std;

typedef unsigned long long int ullint;
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char uchar;

#include "Config.h"

struct Coin
{
	string name;
	string protocol;
	Config config;
	uint global_worksize;
	uint local_worksize;
	uint threads_per_gpu;
	bool max_aggression;
	double sharekhs;

	uint cputhreads;
	string cpu_algorithm;
	string host,port,user,pass,proxy;

	bool check_shares;
	uint nfactor;
	bool nfactor_overridden;
};

struct GlobalConfs
{
	bool save_binaries;
	vector<uint> devices;
	uint platform;
	Coin coin;
	bool restart;
};

enum SHARETEST_VALUE
{
	ST_HNOTZERO,
	ST_MORETHANTARGET,
	ST_GOOD,
};

extern GlobalConfs globalconfs;
extern bool shutdown_now;

const uint KERNEL_INPUT_SIZE = 128;
const uint COPPERLARK_INPUT_SIZE = 31*8;
const uint KERNEL_OUTPUT_SIZE = 256;

const uint YACOIN_PRECALC_SIZE = sizeof(ullint)*30;

const uint WORK_EXPIRE_TIME_SEC = 120;
const uint SHARE_THREAD_RESTART_THRESHOLD_SEC = 20;

const uint TARGET_RUNTIME_MS = 320;
const uint TARGET_RUNTIME_ALLOWANCE_MS = 25;
const uint RUNTIMES_SIZE = 16;

const uint TEMPLATE_ARRAY_SIZE = 16;

const uint CPU_BATCH_SIZE = 4096;

const uint MAX_SIEVE_AMOUNT = 500;

const uint MAX_FERMAT_NUMBERS = 256*1024;

#define foreachgpu() for(vector<_clState>::iterator it = GPUstates.begin(); it != GPUstates.end(); ++it)
#define foreachcpu() for(vector<Reap_CPU_param>::iterator it = CPUstates.begin(); it != CPUstates.end(); ++it)

#if defined(_M_X64) || defined(__x86_64__)
#define REAPER_PLATFORM "64-bit"
#else
#define REAPER_PLATFORM "32-bit"
#endif

#include "App.h"

extern App app;

#endif
