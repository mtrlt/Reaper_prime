#include "Global.h"
#include "CPUAlgos.h"
#include "CPUAlgos_global.h"

vector<uchar> XPM_create_auxdata(mpz_t* bnChainOrigin)
{
	size_t exportsz;
	uchar* ex_port = (uchar*)mpz_export(NULL, &exportsz, -1, 1, -1, 0, *bnChainOrigin);
	//assert(exportsz < 250);  // FIXME: bitcoin varint
	uchar* result = new uchar[exportsz+1];
	result[0] = exportsz;
	result[1] = '\0';
	memcpy(&result[1], ex_port, exportsz);
	free(ex_port);
	vector<uchar> vec(result,result+exportsz+1);
	delete[]result;
	if (vec.size() > 1 && vec.back() >= 0x80)
	{
		++vec[0];
		vec.push_back(0);
	}
	return vec;
}
