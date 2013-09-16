#include "Global.h"
#include "App.h"

#include <stdexcept>

#include "AppOpenCL.h"
#include "Util.h"

App app;

int main(int argc, char* argv[])
{
	try
	{
		vector<string> args;
		for(int i=0; i<argc; ++i)
		{
			args.push_back(argv[i]);
		}
		app.Main(args);
	}
	catch(string s)
	{
		cout << humantime() << "Error: " << s << endl;
	}
	catch(std::runtime_error s)
	{
		cout << humantime() << "Runtime error: " << s.what() << endl;
	}
	catch(std::exception s)
	{
		cout << humantime() << "Exception: " << s.what() << endl;
	}
	catch(...)
	{
		cout << humantime() << "Unknown error." << endl;
	}
}
