#ifndef CPUALGOS_GLOBAL_H
#define CPUALGOS_GLOBAL_H

#include "gmp.h"

#include "mpz_w.h"

vector<uchar> XPM_create_auxdata(mpz_t* bnChainOrigin);

void CPU_Got_share(Reap_CPU_param* state, Work& tempwork,vector<uchar>& auxdata);

static inline
void set_mpz_to_hash(mpz_t *hash, const uint8_t *hashb)
{
	mpz_import(*hash, 8, -1, 4, -1, 0, hashb);
}

extern ullint cpu_shares_hwinvalid;
extern ullint cpu_shares_hwvalid;

extern Work current_work;
extern pthread_mutex_t current_work_mutex;

extern std::vector<unsigned int> vPrimes;
extern std::vector<unsigned int> vTwoInverses;
extern unsigned int nSieveSize;
extern unsigned int nSievePercentage;

unsigned int TargetGetLength(unsigned int nBits);

typedef unsigned long long int uint64;
typedef long long int int64;

extern ullint chainspersec[20];
extern ullint totalpersec;

extern uint found0,foundtotal;
extern ullint fermats,gandalfs;

unsigned int int_invert(unsigned int a, unsigned int nPrime);

#endif