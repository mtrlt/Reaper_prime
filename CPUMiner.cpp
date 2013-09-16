#include "Global.h"
#include "CPUMiner.h"
#include "Util.h"
#include "pthread.h"
#include "Sieve.h"

void CPU_Got_share(Reap_CPU_param* state, Work& tempwork,vector<uchar>& auxdata)
//uchar* tempdata, vector<uchar>& target, uint serverid, uint dataid, vector<uchar>& auxdata)
{
	Share s(vector<uchar>(),vector<uchar>(),tempwork.server_id);
	
	s.data = tempwork.data;
	s.auxdata = auxdata;
	s.dataid = tempwork.dataid;
	s.templateid = tempwork.templateid;
	pthread_mutex_lock(&state->share_mutex);
	state->shares_available = true;
	state->shares.push_back(s);
	pthread_mutex_unlock(&state->share_mutex);
}

bool CPU_Hash_Below_Target(uchar* hash, uchar* target)
{
	cout << "D'ONT CALL THIS FUNCTION YOU BASTARD" << endl;
	for(int i=31; i>=0; --i)
	{
		if (hash[i] > target[i])
			return false;
		else if (hash[i] < target[i])
			return true;
	}
	return true;
}

#include "CPUAlgos.h"
#include "Primes.h"

vector<Reap_CPU_param> CPUstates;

#include <map>
#include "Primes.h"

void GeneratePrimeTable();

void CPUMiner::Init()
{
	if (globalconfs.coin.cputhreads == 0)
	{
#ifdef CPU_MINING_ONLY
		cout << "Config warning: cpu_mining_threads 0" << endl;
#endif
		return;
	}
	const string DEFAULT_ALGORITHM = "hp7";
	
	//hp7 init functions
	GeneratePrimeTable();
	
	map<string,void*(*)(void*)> funcs;
	if(globalconfs.coin.protocol == "primecoin")
	{
		funcs["hp7"] = Reap_CPU_XPM_hp7;
	}
	cout << "Available CPU mining algorithms: ";
	for(map<string,void*(*)(void*)>::iterator it=funcs.begin(); it!=funcs.end(); ++it)
	{
		if (it!=funcs.begin())
			cout << ", ";
		cout << it->first;
	}
	cout << endl;

	void*(*cpualgo)(void*) = funcs[globalconfs.coin.cpu_algorithm];

	if (cpualgo == NULL)
	{
		if (globalconfs.coin.cpu_algorithm != "")
			cout << "CPU algorithm " << globalconfs.coin.cpu_algorithm << " not found." << endl;
		cout << "Using default: " << DEFAULT_ALGORITHM << endl;
		cpualgo = funcs[DEFAULT_ALGORITHM];
	}
	else
		cout << "Using CPU algorithm: " << globalconfs.coin.cpu_algorithm << endl;
	
	for(uint i=0; i<globalconfs.coin.cputhreads; ++i)
	{
		Reap_CPU_param state;
		pthread_mutex_t initializer = PTHREAD_MUTEX_INITIALIZER;

		state.share_mutex = initializer;
		state.shares_available = false;

		state.hashes = 0;

		state.thread_id = i;

		CPUstates.push_back(state);
	}
	
	cout << "Creating " << CPUstates.size() << " CPU thread" << (CPUstates.size()==1?"":"s") << "." << endl;
	for(uint i=0; i<CPUstates.size(); ++i)
	{
		cout << i+1 << "...";
		pthread_attr_t attr;
	    pthread_attr_init(&attr);
		int schedpolicy;
		pthread_attr_getschedpolicy(&attr, &schedpolicy);
		int schedmin = sched_get_priority_min(schedpolicy);
		int schedmax = sched_get_priority_max(schedpolicy);
		if (i==0 && schedmin == schedmax)
		{
			cout << "Warning: can't set thread priority" << endl;
		}
		sched_param schedp;
		schedp.sched_priority = schedmin;
		pthread_attr_setschedparam(&attr, &schedp);

		pthread_create(&CPUstates[i].thread, &attr, cpualgo, (void*)&CPUstates[i]);
		pthread_attr_destroy(&attr);
	}
	cout << "done" << endl;
}

void CPUMiner::Quit()
{
}
