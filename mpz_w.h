#ifndef MPZ_WRAPPER
#define MPZ_WRAPPER

#include <gmp.h>

class Mpz_w
{
public:
	mpz_t n;
	
	Mpz_w(){ mpz_init(n); }
	~Mpz_w(){ mpz_clear(n); }
};

#endif