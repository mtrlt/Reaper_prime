#include "Global.h"
#include "App.h"
#include "Util.h"

#include "SHA256.h"

uchar HexToChar(char data)
{
	if (data <= '9')
		return data-'0';
	else if (data <= 'Z')
		return data-'7';
	else
		return data-'W';
}

#include "AppOpenCL.h"

#ifndef CPU_MINING_ONLY
extern vector<_clState> GPUstates;
#endif
extern vector<Reap_CPU_param> CPUstates;

uchar HexToChar(char h, char l)
{
	return HexToChar(h)*16+HexToChar(l);
}

const char* hextable[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f"};
string CharToHex(uchar c)
{
	return string(hextable[c/16]) + string(hextable[c%16]);
}

vector<uchar> HexStringToVector(string str)
{
	vector<uchar> ret;
	ret.assign(str.length()/2, 0);
	for(uint i=0; i<str.length(); i+=2)
	{
		ret[i/2] = HexToChar(str[i+0], str[i+1]);
	}
	return ret;
}
string VectorToHexString(vector<uchar> vec)
{
	string ret;
	for(uint i=0; i<vec.size(); i++)
	{
		ret += CharToHex(vec[i]);
	}
	return ret;
}

void DoubleSHA256(uint* output2, uint* workdata, uint* midstate);

void LineClear()
{
	cout << "\r                                                                      \r";
}

#ifdef WIN32
void Wait_ms(uint n);
#undef SetPort
#else
void Wait_ms(uint n);
#endif

#include "Config.h"

Config config;
GlobalConfs globalconfs;

ullint shares_valid = 0;
ullint shares_invalid = 0;
ullint shares_hwvalid = 0;
ullint shares_hwinvalid = 0;
ullint cpu_shares_hwvalid = 0;
ullint cpu_shares_hwinvalid = 0;
clock_t current_work_time = 0;

extern Work current_work;

SHARETEST_VALUE ShareTest_BTC(uint* workdata, uint* target);

bool getwork_now = false;

unsigned char GetNfactor(long long int nTimestamp);

#define dump(x) cout << #x << ": " << x << endl;

void SubmitShare(Curl& curl, Share& w)
{
	try
	{
		string str;
		{
			blktemplate_t* tmpl = app.templates[w.templateid];
			uint NONCE = EndianSwap(*(uint*)&w.data[76]);
			
			json_t* readyblock = blkmk_submit_jansson(tmpl, &w.data[0], w.dataid, NONCE, &w.auxdata[0], w.auxdata.size());
			char *s = json_dumps(readyblock, JSON_INDENT(2));
			str = s;
			free(s);
		}
		string ret = curl.TestWork(servers[w.server_id],str);
	}
	catch(std::exception s)
	{
		cout << "(3) Error: " << s.what() << endl;
	}
}

bool sharethread_active;
void* ShareThread(void* param)
{
	cout << "Share thread started" << endl;
	Curl* pcurl;

	while(!shutdown_now)
	{
		sharethread_active = true;
		Wait_ms(50);
#ifndef CPU_MINING_ONLY
		foreachgpu()
		{
			sharethread_active = true;
			if (!it->shares_available)
				continue;
			Share s;
			pthread_mutex_lock(&it->share_mutex);
			if (it->shares.empty())
			{
				it->shares_available = false;
				pthread_mutex_unlock(&it->share_mutex);
				continue;
			}
			s = it->shares.front();
			it->shares.pop_front();
			if (it->shares.empty())
				it->shares_available = false;
			pthread_mutex_unlock(&it->share_mutex);
			uint currentserverid = s.server_id;
			{
				uint tiimm = ticker();
				SubmitShare(*pcurl, s);
				tiimm = ticker()-tiimm;
				if (tiimm > 5000)
					cout << "Share submit took " << tiimm/1000.0 << " s!  " << endl;
			}
		}
#endif
		foreachcpu()
		{
			sharethread_active = true;
			if (!it->shares_available)
				continue;
			Share s;
			pthread_mutex_lock(&it->share_mutex);
			if (it->shares.empty())
			{
				pthread_mutex_unlock(&it->share_mutex);
				continue;
			}
			s = it->shares.front();
			it->shares.pop_front();
			if (it->shares.empty())
				it->shares_available = false;
			pthread_mutex_unlock(&it->share_mutex);
			uint currentserverid = s.server_id;
			{
				uint tiimm = ticker();
				SubmitShare(*pcurl, s);
				tiimm = ticker()-tiimm;
				if (tiimm > 5000)
					cout << "Share submit took " << tiimm/1000.0 << " s!  " << endl;
			}
		}
	}
	pthread_exit(NULL);
	return NULL;
}

extern string longpoll_url;
extern bool longpoll_active;

struct LongPollThreadParams
{
	Curl* curl;
	App* app;
};

void* LongPollThread(void* param)
{
	if(servers.size() != 1)
	{
		cout << "Long polling disabled when using many servers" << endl;
		pthread_exit(NULL);
		return NULL;
	}
	cout << "Long polling disabled." << endl;
	pthread_exit(NULL);
	return NULL;
	LongPollThreadParams* p = (LongPollThreadParams*)param; 

	Curl* curl = p->curl;
	
	string LP_url = longpoll_url;
	string LP_path;

	cout << "Long polling URL: [" << LP_url << "]. trying to parse." << endl;
	clock_t lastcall = 0;

	{//parsing LP address
		vector<string> exploded = Explode(LP_url, '/');
		if (exploded.size() >= 2 && exploded[0] == "http:")
		{
			vector<string> exploded2 = Explode(exploded[1], ':');
			if (exploded2.size() != 2)
				goto couldnt_parse;
			cout << "LP Host: " << exploded2[0] << endl;
			//curl.SetHost(exploded2[0]);
			cout << "LP Port: " << exploded2[1] << endl;
			//curl.SetPort(exploded2[1]);
			if (exploded.size() <= 2)
				LP_path = '/';
			else
				LP_path = "/" + exploded[2];
			cout << "LP Path: " << LP_path << endl;
		}
		else if (LP_url.length() > 0 && LP_url[0] == '/')
		{
			LP_path = LP_url;
			cout << "LP Path: " << LP_path << endl;
		}
		else
		{
			goto couldnt_parse;
		}
	}

	while(!shutdown_now)
	{
		clock_t ticks = ticker();
		if (ticks-lastcall < 5000)
		{
			Wait_ms(ticks-lastcall);
		}
		lastcall = ticks;
		string r = curl->GetWork(servers[0], LP_path, 60);
		cout << "Got LP" << endl;
		p->app->Parse(r);
	}
	pthread_exit(NULL);
	return NULL;

couldnt_parse:
	cout << "Couldn't parse long polling URL [" << LP_url << "]. turning LP off." << endl;
	pthread_exit(NULL);
	return NULL;
}

bool shutdown_now=false;
void* ShutdownThread(void* param)
{
	cout << "Press [Q] and [Enter] to quit" << endl;
	while(shutdown_now == false)
	{
		string s;
		std::cin >> s;
		if (s == "q" || s == "Q")
			shutdown_now = true;
	}
	cout << "Quitting." << endl;
	pthread_exit(NULL);
	return NULL;
}

#include <sstream>
using std::stringstream;

void App::SetupCurrency()
{
	map<string,Coin> coins;
	{
		Coin c;
		c.name = config.GetValue<string>("mine");
		c.config.Load(c.name + ".conf");
		c.protocol = c.config.GetValue<string>("protocol");
		c.local_worksize = c.config.GetValue<uint>("worksize");
		{
			if (c.config.GetValue<string>("aggression") == "max")
			{
				c.global_worksize = 1<<11;
				c.max_aggression = true;
			}
			else
			{
				c.global_worksize = (1<<c.config.GetValue<uint>("aggression"));
				c.max_aggression = false;
			}
			c.global_worksize /= c.local_worksize;
			c.global_worksize *= c.local_worksize;
		}
		c.threads_per_gpu = c.config.GetValue<uint>("threads_per_gpu");
		c.cputhreads = c.config.GetValue<uint>("cpu_mining_threads");
		c.cpu_algorithm = c.config.GetValue<string>("cpu_algorithm");
		c.check_shares = c.config.GetValue<bool>("check_shares");

		c.host = c.config.GetValue<string>("host");
		c.port = c.config.GetValue<string>("port");
		c.user = c.config.GetValue<string>("user");
		c.pass = c.config.GetValue<string>("pass");
		c.proxy = c.config.GetValue<string>("proxy");
		if (c.local_worksize > c.global_worksize)
			c.global_worksize = c.local_worksize;

		coins[c.name] = c;
	}
	string minedcoin = config.GetValue<string>("mine");
	if (coins.find(minedcoin) == coins.end())
	{
		cout << "Coin chosen for mining \"" << minedcoin << "\" not found." << endl;
		throw string("");
	}
	else
		cout << "I'm now mining " << minedcoin << "!" << endl;
	globalconfs.coin = coins[minedcoin];
	if (globalconfs.coin.host == "" ||
		globalconfs.coin.port == "" ||
		globalconfs.coin.user == "" ||
		globalconfs.coin.pass == "")
		throw string("Config ") + globalconfs.coin.name + ".conf is missing one of host/port/user/pass.";

	globalconfs.coin.sharekhs = 0;

	if (globalconfs.coin.cputhreads == 0 && globalconfs.coin.threads_per_gpu == 0)
	{
		throw string("No CPU or GPU mining threads.. please set either cpu_mining_threads or threads_per_gpu to something other than 0.");
	}
}

vector<ServerSettings> servers;
void App::LoadServers()
{
	uint servercount = globalconfs.coin.config.GetValueCount("host");
	for(uint i=0; i<servercount; ++i)
	{
		ServerSettings s;
		s.host = globalconfs.coin.config.GetValue<string>("host",i);
		s.port = globalconfs.coin.config.GetValue<unsigned short>("port",i);
		s.user = globalconfs.coin.config.GetValue<string>("user",i);
		s.pass = globalconfs.coin.config.GetValue<string>("pass",i);
		s.proxy = globalconfs.coin.config.GetValue<string>("proxy",i);
		servers.push_back(s);
	}
	current_server_id = 0;
}

#if (BLKMAKER_VERSION < 2)
#error libblkmaker must be at least version 2!
#endif

string SIize(double num)
{
	if (num<0) //NOT IMPLEMENTED YET
		return ToString(num);
	
	if (num<1e3)
		return ToString(num);
	if (num<1e6)
		return ToString(num/1e3)+"k";
	if (num<1e9)
		return ToString(num/1e6)+"M";
	if (num<1e12)
		return ToString(num/1e9)+"G";
	return ToString(num/1e12)+"T";
}


extern ullint chainspersec[20];
extern ullint totalpersec;
extern ullint fermats,gandalfs;

extern std::vector<unsigned int> vPrimes;
extern std::vector<unsigned int> vTwoInverses;

extern unsigned int nSieveSize;

void App::Main(vector<string> args)
{
	cout << "\\|||||||||||||||||||||/" << endl;
	cout << "-  Reaper " << REAPER_VERSION << " " << REAPER_PLATFORM << "  -" << endl;
	cout << "-    PRIME BETA 2     -" << endl;
	cout << "-   coded by mtrlt    -" << endl;
	cout << "/|||||||||||||||||||||\\" << endl;
	cout << endl;
	cout << endl;
	
	string config_name = "reaper.conf";
	if (args.size() > 2)
	{
		cout << "Please use config files to set host/port/user/pass" << endl;
		return;
	}
	if (args.size() == 2)
		config_name = args[1];
	getworks = 0;
	config.Load(config_name);
	blkmk_sha256_impl = mysha256;

	for(uint i=0; i<TEMPLATE_ARRAY_SIZE; ++i)
		templates[i] = NULL;
	SetupCurrency();
	LoadServers();
	Wait_ms(100);
	nickbase = globalconfs.coin.user;
	globalconfs.restart = false;
	globalconfs.save_binaries = config.GetValue<bool>("save_binaries");
	uint numdevices = config.GetValueCount("device");
	for(uint i=0; i<numdevices; ++i)
		globalconfs.devices.push_back(config.GetValue<uint>("device", i));

#ifdef CPU_MINING_ONLY
	if (globalconfs.coin.cputhreads == 0)
	{
		throw string("cpu_mining_threads is zero. Nothing to do, quitting.");
	}
#endif
	globalconfs.platform = config.GetValue<uint>("platform");

	current_work.old = true;
	current_work.time = 0;
	
	template_nonce = 0;

	Curl::GlobalInit();

	pthread_t sharethread;
	pthread_create(&sharethread, NULL, ShareThread, &curl);
	if (globalconfs.coin.config.GetValue<uint>("sharethreads"))
	{
		for(uint i=1; i<globalconfs.coin.config.GetValue<uint>("sharethreads"); ++i)
			pthread_create(&sharethread, NULL, ShareThread, &curl);
	}

	cpuminer.Init();
	{
		map<string,string> defines;
		defines["NSIEVESIZE"] = ToString(nSieveSize);
	
		vector<string> kernelnames;
		kernelnames.push_back("CalculateMultipliers");
		kernelnames.push_back("Sieve");
		kernelnames.push_back("Combine");
		kernelnames.push_back("Fermat");
		
		const unsigned int zSIZE = 4;
		const unsigned int ySIZE = vPrimes.size();
		const unsigned int xSIZE = 10; //max half chain length
		
		unsigned int SieveSizeBytes = nSieveSize/8;
		
		map<string,cl_mem_settings> buffersettings;
		buffersettings["vPrimes"] = cl_mem_settings(vPrimes.size()*sizeof(unsigned int),CL_MEM_READ_ONLY);
		buffersettings["vfCandidates"] = cl_mem_settings(SieveSizeBytes,CL_MEM_READ_WRITE);
		buffersettings["vfCompositeCunningham"] = cl_mem_settings(4*SieveSizeBytes,CL_MEM_READ_WRITE);
		buffersettings["vfCunninghamMultipliers"] = cl_mem_settings(sizeof(uint)*zSIZE*ySIZE*xSIZE,CL_MEM_READ_WRITE);
		buffersettings["FermatNumbers"] = cl_mem_settings(sizeof(uint)*MAX_FERMAT_NUMBERS,CL_MEM_READ_WRITE);
		buffersettings["FermatOutput"] = cl_mem_settings(sizeof(uint)*MAX_FERMAT_NUMBERS,CL_MEM_READ_WRITE);
		buffersettings["Counter"] = cl_mem_settings(sizeof(uint),CL_MEM_READ_WRITE);
		buffersettings["CandidateList"] = cl_mem_settings(sizeof(uint)*(256*1024),CL_MEM_READ_WRITE);
		
		opencl.Init(kernelnames, buffersettings, defines);
		
		opencl.WriteBufferGlobal("vPrimes",&vPrimes[0],vPrimes.size()*sizeof(unsigned int));
		
		opencl.inited = true;
	}
	current_server_id = (current_server_id+1)%servers.size();
	current_templateid = 0;
	Parse(curl.GetWork(servers[current_server_id]));

	int work_update_period_ms = globalconfs.coin.config.GetValue<uint>("getwork_rate");

	if (work_update_period_ms == 0)
		work_update_period_ms = 1500;

	pthread_t longpollthread;
	LongPollThreadParams lp_params;
	if (longpoll_active)
	{
		cout << "Activating long polling." << endl;
		lp_params.app = this;
		lp_params.curl = &curl;
		pthread_create(&longpollthread, NULL, LongPollThread, &lp_params);
	}

	if (config.GetValue<bool>("enable_graceful_shutdown"))
	{
		pthread_t shutdownthread;
		pthread_create(&shutdownthread, NULL, ShutdownThread, NULL);
	}

	clock_t ticks = ticker();
	clock_t starttime = ticker();
	workupdate = ticker();
	clock_t lastprimeprinttime = ticker();

	clock_t sharethread_update_time = ticker();

	while(!shutdown_now)
	{
		Wait_ms(100);
		clock_t timeclock = ticker();
		if (timeclock - current_work_time >= WORK_EXPIRE_TIME_SEC*1000)
		{
			if (!current_work.old)
			{
				cout << humantime() << "Work too old... waiting for getwork.    " << endl;
			}
			current_work.old = true;
		}

		if (sharethread_active)
		{
			sharethread_active = false;
			sharethread_update_time = timeclock;
		}
		if (timeclock-sharethread_update_time >= SHARE_THREAD_RESTART_THRESHOLD_SEC*1000)
		{
			cout << humantime() << "Share thread messed up. Starting another one.   " << endl;
			pthread_create(&sharethread, NULL, ShareThread, &curl);
		}
		if (getwork_now || timeclock - workupdate >= work_update_period_ms)
		{
			uint timmii = ticker();
			current_server_id = (current_server_id+1)%servers.size();
			Parse(curl.GetWork(servers[current_server_id]));
			timmii = ticker()-timmii;
			if (timmii > 5000)
				cout << "Getwork took " << timmii/1000.0 << " s!  " << endl;
			getwork_now = false;
		}
		if (timeclock - ticks >= 1000)
		{
			ullint totalhashesGPU=0;
#ifndef CPU_MINING_ONLY
			foreachgpu()
			{
				totalhashesGPU += it->hashes;
			}
#endif
			ullint totalhashesCPU=0;
			foreachcpu()
			{
				totalhashesCPU += it->hashes;
			}

			ticks += (timeclock-ticks)/1000*1000;
			float stalepercent = 0.0f;
			if (shares_valid+shares_invalid != 0)
				stalepercent = 100.0f*float(shares_invalid)/float(shares_invalid+shares_valid);
				
			if (ticks-lastprimeprinttime >= 4000)
			{
				lastprimeprinttime += 4000;
				
				cout << fermats/double((ticker()-starttime)*0.001) << " fermats/s, " << gandalfs/double((ticker()-starttime)*0.001) << " gandalfs/s.";
				
				cout << " | Per hour:  ";
				for(uint i=2; i<15; ++i)
				{
					if (chainspersec[i] > 0)
					{
						double num = 3600.0*chainspersec[i]/((ticker()-starttime)*0.001);
						string str = SIize(num);
						cout << str << " " << i << "-ch  ";
					}
				}
				cout << endl;
				
				
				//cout << foundprimes/(double(ticker()-starttime)/1000.0) << " primes per sec" << endl;
			}
			/*if (ticks-starttime == 0)
				cout << dec << "   ??? kH/s, shares: " << shares_valid << "|" << shares_invalid << ", stale " << stalepercent << "%, " << (ticks-starttime)/1000 << "s    \r";
			else
			{
				stringstream stream;
				stream.precision(4);
				if (totalhashesGPU != 0)
				{
					stream << "GPU ";
					double spd = double(totalhashesGPU)/(ticks-starttime);
					if (spd >= 1000000.0)
						stream << spd/1000000.0 << "GH/s, ";
					else if (spd >= 1000.0)
						stream << spd/1000.0 << "MH/s, ";
					else
						stream << spd << "kH/s, ";
				}
				if (totalhashesCPU != 0)
				{
					stream << "CPU " << double(totalhashesCPU)/(ticks-starttime) << "kH/s, ";
				}
				stream << "shares: " << shares_valid << "|" << shares_invalid << ", stale " << stalepercent << "%, ";
				float hwpercent = 100.0f*float(shares_hwinvalid)/float(shares_hwinvalid+shares_hwvalid);
				if (shares_hwinvalid+shares_hwvalid > 0)
					stream << "GPU errors: " << hwpercent << "%, ";
				float cpuhwpercent = 100.0f*float(cpu_shares_hwinvalid)/float(cpu_shares_hwinvalid+cpu_shares_hwvalid);
				if (cpu_shares_hwinvalid+cpu_shares_hwvalid > 0)
					stream << "CPU errors: " << cpuhwpercent << "%, ";
				cout << dec << stream.str() << "~" << 1000.0*(shares_hwvalid+cpu_shares_hwvalid)/double(ticks-starttime)*globalconfs.coin.sharekhs << " kH/s, " << (ticks-starttime)/1000 << "s  \r";
				cout.flush();
			}*/
		}
	}
	cpuminer.Quit();
	opencl.Quit();
	Curl::GlobalQuit();
}

bool targetprinted=false;

pthread_mutex_t current_work_mutex = PTHREAD_MUTEX_INITIALIZER;
Work current_work;

void App::Parse(string data)
{
	workupdate = ticker();
	if (data == "")
	{
		cout << humantime() << "Couldn't connect to server. ";
		if (servers.size() > 1)
			cout << "Trying next server in a few seconds... " << endl;
		else
			cout << "Trying again in a few seconds... " << endl;

		return;
	}

	Parse_XPM(data);
}

struct _cbscript_t {
	char *data;
	size_t sz;
};

const char *set_b58addr(const char *arg, struct _cbscript_t *p)
{
	size_t scriptsz = blkmk_address_to_script(NULL, 0, arg);
	if (!scriptsz)
	{
		return "Invalid address";
	}
	char *script = (char*)malloc(scriptsz);
	if (blkmk_address_to_script(script, scriptsz, arg) != scriptsz) 
	{
		free(script);
		return "Failed to convert address to script";
	}
	p->data = script;
	p->sz = scriptsz;
	return NULL;
}

#include "gmp.h"


void Precalc_BTC(Work& work, uint vectors);
vector<uchar> CalculateMidstate(vector<uchar> in);
void App::Parse_XPM(string data)
{
	json_error_t jsone;
	
	json_t* req = json_loads(&data[0], 0, &jsone);
	
	/*
	{
		char *s = json_dumps(req, JSON_INDENT(2));
		cout << "DATA: " << s << endl;
		free(s);
	}
	*/
	
	
	blktemplate_t* tmpl;
	tmpl = blktmpl_create();
	
	const char* err = blktmpl_add_jansson(tmpl, req, time(NULL));
	json_decref(req);
	if (err)
	{
		cout << "Error adding block template: " << err << endl;
	}
	unsigned char header[80] = {}, hash[32] = {};
	unsigned int dataid=0;
	
	//TODO: CBTXN SHIT
	{
		_cbscript_t opt_coinbase_script;
		string address = globalconfs.coin.config.GetValue<string>("primecoin_address");
		if (address == "")
			throw string("Please set primecoin_address in primecoin.conf to the address you want to mine to.");

		const char* ret = set_b58addr(address.c_str(),&opt_coinbase_script);
		if (ret != NULL)
		{
			cout << "RET! " << ret << endl;
		}
		if (opt_coinbase_script.sz)
		{
			bool newcb;
			ullint num = blkmk_init_generation2(tmpl, opt_coinbase_script.data, opt_coinbase_script.sz, &newcb);
			/*if (newcb)
			{
				ssize_t ae = blkmk_append_coinbase_safe(tmpl, &template_nonce, sizeof(template_nonce));
				if (ae < (ssize_t)sizeof(template_nonce))
				{
					cout << "Cannot append template-nonce to coinbase.." << ae << "..you might be wasting hashing! Or priming! Whatever!" << endl;
				}
				++template_nonce;
			}*/
		}
		else
			cout << "There was no coinbase! Or something." << endl;
	}
	
	unsigned int datasz = blkmk_get_data(tmpl, header, 80, time(NULL), NULL, &dataid);
	
	if (templates[current_templateid] != NULL)
	{
		blktmpl_free(templates[current_templateid]);
	}
	templates[current_templateid] = tmpl;

	Work newwork;
	newwork.data = vector<uchar>(header,header+80);
	//cout << VectorToHexString(newwork.data) << endl;
	newwork.old = false;
	newwork.time = ticker();
	newwork.dataid = dataid;
	newwork.templateid = current_templateid;
	newwork.server_id = current_server_id;
	current_templateid = (current_templateid+1)%TEMPLATE_ARRAY_SIZE;
	current_work_time = ticker();
	
	current_work.time = ticker();
	pthread_mutex_lock(&current_work_mutex);
	current_work = newwork;
	pthread_mutex_unlock(&current_work_mutex);
	current_server_id = (current_server_id+1)%servers.size();
	
	
	//json_t* readyblock = blkmk_submit_jansson(work->tmpl, data, work->dataid, le32toh(*((uint32_t*)&work->data[76])), work->sig, work->sigsz);
	//tmpl is the template
	//data is the block header, 80 bytes.
	//dataid is something.
	//letoh is the nonce
	//work->sig is the extra data used by primecoin
	//work->sigsz is the sig size
	/*{
		mpz_t tempnum;
		mpz_init_set_str(tempnum, "100", 10);
		vector<uchar> auxdata = XPM_create_auxdata(&tempnum);
		
		uint NONCE = EndianSwap(*(uint*)&data[76]);
		
		json_t* readyblock = blkmk_submit_jansson(tmpl, &newwork.data[0], dataid, NONCE, &auxdata[0], auxdata.size());

		char *s = json_dumps(readyblock, JSON_INDENT(2));
		cout << "BLOCK: " << s << endl;
		free(s);
	}*/
/*
	Json::Value root, result, error;
	Json::Reader reader;
	bool parsing_successful = reader.parse( data, root );
	if (!parsing_successful)
	{
		goto got_error;
	}

	result = root.get("result", "null");
	error = root.get("error", "null");
	
	if (result.isObject())
	{
		Json::Value::Members members = result.getMemberNames();
		uint neededmembers=0;
		bool has_midstate = false;
		for(Json::Value::Members::iterator it = members.begin(); it != members.end(); ++it)
		{
			if (*it == "data")
				++neededmembers;
			if (*it == "midstate")
				has_midstate = true;
		}
		if (neededmembers != 1 || !result["data"].isString())
		{
			goto got_error;
		}

		++getworks;
		Work newwork;
		if (has_midstate)
			newwork.midstate = HexStringToVector(result["midstate"].asString());
		else
			newwork.midstate = CalculateMidstate(HexStringToVector(result["data"].asString().substr(0,128)));
		newwork.data = HexStringToVector(result["data"].asString().substr(0, 160));
		newwork.old = false;
		newwork.time = ticker();
		current_work_time = ticker();

		if (!targetprinted)
		{
			targetprinted = true;
			cout << result["target"].asString() << endl;
		}

		{
			vector<uchar> targetuchar = HexStringToVector(result["target"].asString());
			newwork.target_share.assign(32,0);
			for(uint i=0; i<8; ++i)
			{
				uint number = (((uint*)(&targetuchar[0]))[i]);
				newwork.target_share[4*i] = number;
				newwork.target_share[4*i+1] = number>>8;
				newwork.target_share[4*i+2] = number>>16;
				newwork.target_share[4*i+3] = number>>24;
			}
			newwork.server_id = current_server_id;
		}

		Precalc_BTC(newwork,opencl.GetVectorSize());
		current_work.time = ticker();
		pthread_mutex_lock(&current_work_mutex);
		current_work = newwork;
		pthread_mutex_unlock(&current_work_mutex);
		current_server_id = (current_server_id+1)%servers.size();
		return;
	}
	else if (!error.isNull())
	{
		cout << error.asString() << endl;
		cout << "Code " << error["code"].asInt() << ", \"" << error["message"].asString() << "\"" << endl;
	}
got_error:
	cout << "Error with server: " << data << endl;
	return;
*/
}
