#ifndef APP_H
#define APP_H

#include "Curl.h"
#include "AppOpenCL.h"
#include "CPUMiner.h"


extern "C" {

#include <blkmaker.h>
#include <blkmaker_jansson.h>
#include <blktemplate.h>

}

class App
{
public:
	struct ServerData
	{
		uint current_id;
		uint last_tried;
	};

	uint current_server_id;
	uint current_templateid;
	string nickbase;

	Curl curl;
	OpenCL opencl;
	CPUMiner cpuminer;

	clock_t workupdate;

	uint getworks;
	
	uint template_nonce;

	void SetupCurrency();

	void Parse_XPM(string data);

	void Restart();
	
	blktemplate_t* templates[TEMPLATE_ARRAY_SIZE];

	void Main(vector<string> args);
	void Parse(string data);
	void LoadServers();
};

#endif
