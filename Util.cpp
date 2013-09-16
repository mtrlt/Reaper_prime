#include "Global.h"
#include "Util.h"

#include <sstream>
using std::stringstream;

template<> 
bool FromString<bool>(string key)
{
	if (key == "true" || key == "on" || key == "yes")
		return true;
	if (key == "false" || key == "off" || key == "no")
		return false;
	if (key == "")
		return false;
	return false;
}

string ToString(bool key, string truestring, string falsestring)
{
	if (key)
		return truestring;
	else
		return falsestring;
}

template<> int FromString<int>(string key)
{
	if (key == "")
		return 0;
	stringstream sstr(key);
	int ret=0;
	sstr >> ret;
	return ret;
}

uint EndianSwap(uint n)
{
	return ((n&0xFF)<<24) | ((n&0xFF00)<<8) | ((n&0xFF0000)>>8) | ((n&0xFF000000)>>24);
}

#include <ctime>

string humantime()
{
	time_t rawtime;
	char formattedtime[100];
	tm* timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(formattedtime, 100, "%Y-%m-%d %H:%M:%S ", timeinfo);
	string ret = formattedtime;
	return ret;
}

#ifdef WIN32
#include "windows.h"
clock_t ticker()
{
	//todo: convert to QueryPerformanceCounter
	return clock()/(CLOCKS_PER_SEC/1000);
}

void Wait_ms(uint n)
{
	Sleep(n);
}
#else
#include <unistd.h>
#include <sys/time.h>
clock_t ticker()
{
	timeval t;
	gettimeofday(&t, NULL);
	unsigned long long l = ((unsigned long long)(t.tv_sec))*1000 + t.tv_usec/1000;
	return l;
}

void Wait_ms(uint n)
{
	timespec ts;
	ts.tv_sec = n/1000;
	ts.tv_nsec = n*1000000;
	nanosleep(&ts, NULL);
}
#endif

vector<string> Explode(string s, char delim)
{
	vector<string> returner;
	returner.clear();
	if (s.length() == 0)
		return returner;
	string temp;
	for(uint i=0; i<s.length(); ++i)
	{
		if (s[i] == delim)
		{
			if (temp.size() > 0)
				returner.push_back(temp);
			temp.clear();
		}
		else
		{
			temp.push_back(s[i]);
		}
	}
	returner.push_back(temp);
	return returner;
}
