#ifndef CURLMUNACPP
#define CURLMUNACPP

#include "ServerSettings.h"

class Curl
{
private:
	enum EXEC_TYPE
	{
		GETWORK,
		GETWORK_LP,
		TESTWORK,
	};

	string Execute(ServerSettings& s, Curl::EXEC_TYPE type, string work, string path, uint timeout);
	void Execute_XPM(void* curl, Curl::EXEC_TYPE type, string work, string path, uint timeout);

public:

	Curl() {}
	~Curl() {}

	static void GlobalInit();
	static void GlobalQuit();

	void* Init();
	void Quit(void* curl);

	string GetWork_LP(ServerSettings& s, string path="", uint timeout = 60);
	string GetWork(ServerSettings& s, string path="", uint timeout = 5);
	string TestWork(ServerSettings& s, string work);
};

#endif
