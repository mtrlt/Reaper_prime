#ifndef CPUMINER_H
#define CPUMINER_H

#include "pthread.h"

struct Reap_CPU_param
{
	uint thread_id;

	pthread_t thread;

	bool shares_available;
	deque<Share> shares;
	pthread_mutex_t share_mutex;

	ullint hashes;
};

class CPUMiner
{
private:
public:
	void Init();
	void Quit();
};

#endif
