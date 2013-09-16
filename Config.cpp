#include "Global.h"
#include "Config.h"

#include <fstream>

#include <map>
using std::map;
using std::pair;
using std::ifstream;
using std::ofstream;

#include <cstdio>
#include <algorithm>

void Config::Clear()
{
	config.clear();
}

void Config::Load(string filename, vector<string> included_already)
{
	included_already.push_back(filename);
	{
		FILE* filu = fopen(filename.c_str(), "r");
		if (filu == NULL)
			throw string("Config file " + filename + " not found.");
		fclose(filu);
	}

	map<string, string> config_values;
	config_values["aggression"] = "uint";
	config_values["worksize"] = "uint";
	config_values["threads_per_gpu"] = "uint";
	config_values["device"] = "uint array";
	config_values["kernel"] = "string";
	config_values["save_binaries"] = "bool";
	config_values["cpu_mining_threads"] = "uint";
	config_values["platform"] = "uint";
	config_values["enable_graceful_shutdown"] = "bool";
	config_values["host"] = "string";
	config_values["port"] = "string";
	config_values["user"] = "string";
	config_values["pass"] = "string";
	config_values["proxy"] = "string";
	config_values["long_polling"] = "bool";
	config_values["include"] = "string";
	config_values["cpu_algorithm"] = "string";
	config_values["protocol"] = "string";
	config_values["mine"] = "string";
	config_values["gpu_sharemask"] = "string";
	config_values["vectors"] = "uint";
	config_values["sharethreads"] = "uint";
	config_values["dont_check_shares"] = "bool";
	config_values["lookup_gap"] = "uint";
	config_values["gpu_thread_concurrency"] = "uint";
	config_values["check_shares"] = "bool";
	config_values["nfactor_override"] = "uint";
	config_values["getwork_rate"] = "uint";
	config_values["sieve_size"] = "uint";
	config_values["min_chain_length"] = "uint";
	config_values["sievesize"] = "uint";
	config_values["sievepercentage"] = "uint";
	config_values["primecoin_address"] = "string";
	config_values["debug"] = "bool";
	
	ifstream filu(filename.c_str());
	while(!filu.eof())
	{
		string prop;
		filu >> prop;

		string value;
		filu >> value;

		if (prop == "include")
		{
			if (std::find(included_already.begin(), included_already.end(), value) != included_already.end())
			{
				cout << "Circular include: ";
				for(uint i=0; i<included_already.size(); ++i)
					cout << included_already[i] << " -> ";
				cout << value << endl;
			}
			else
			{
				Config includedconfig;
				includedconfig.Load(value, included_already);
				for(map<string,vector<string> >::iterator it = includedconfig.config.begin(); it != includedconfig.config.end(); ++it)
					for(uint i=0; i<it->second.size(); ++i)
						config[it->first] = it->second;
			}
		}

		if (config_values.find(prop) == config_values.end())
		{
			bool fail = true;
			CombiKey c = GetCombiKey(prop);
			if (c.base != "" && c.id != -1 && c.prop != "")
			{
				fail = false;
			}
			if (fail)
			{
				if (prop != "")
				{
					cout << "Warning: unknown property \"" << prop << "\" in configuration file." << endl;
				}
			}
		}
		config[prop].push_back(value);
	}
	included_already.pop_back();
}
