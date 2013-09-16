#include "Global.h"
#include "CPUMiner.h"
#include "CPUAlgos.h"
#include "Util.h"
#include "SHA256.h"
#include "Sieve.h"
#include "Primes.h"

#include "CPUAlgos_global.h"

//void CPU_Got_share(Reap_CPU_param* state, uchar* tempdata, vector<uchar>& target, uint serverid);
bool CPU_Hash_Below_Target(uchar* hash, uchar* target);



#define nFractionalBits 24
#define TARGET_FRACTIONAL_MASK ((1u << nFractionalBits) - 1)
#define TARGET_LENGTH_MASK (~TARGET_FRACTIONAL_MASK)

//A bunch of helper functions
static
void TargetIncrementLength(unsigned int *pnBits)
{
    *pnBits += (1 << nFractionalBits);
}

static
unsigned int TargetFromInt(unsigned int nLength)
{
    return (nLength << nFractionalBits);
}

// Check Fermat probable primality test (2-PRP): 2 ** (n-1) = 1 (mod n)
// true: n is probable prime
// false: n is composite; set fractional length in the nLength output
static
bool PrimeTest(mpz_t *n, unsigned int *pnLength)
{
	//fast tests
	/*{
		for(uint i=MAX_SIEVE_AMOUNT; i<MAX_SIEVE_AMOUNT+100; ++i)
		{
			if (mpz_divisible_ui_p(*n,Primes::v[i]))
			{
				return false;
			}
		}
	}*/


	mpz_t a, e, r;
	mpz_init_set_ui(a, 2); // base; Fermat witness
	mpz_init(e);
	mpz_sub_ui(e, *n, 1);
	mpz_init(r);
	
	mpz_powm(r, a, e, *n); // r = (2**(n-1))%n
	mpz_clear(a);
	mpz_clear(e);
	if (!mpz_cmp_ui(r, 1)) //if r is 1, this is a witness!
	{
		mpz_clear(r);
		return true;
	}
	
	// Failed Fermat test, calculate fractional length
	// nFractionalLength = ( (n-r) << nFractionalBits ) / n
	mpz_sub(r, *n, r);
	mpz_mul_2exp(r, r, nFractionalBits);
	mpz_fdiv_q(r, r, *n);
	unsigned int nFractionalLength = mpz_get_ui(r);
	mpz_clear(r);
	
	if (nFractionalLength >= (1 << nFractionalBits))
	{
		cout << "PrimeTest() : fractional assert" << endl;
		return false;
	}
	*pnLength = (*pnLength & TARGET_LENGTH_MASK) | nFractionalLength;
	return false;
}

//This function determines the primality of a number of the form 2p+-1
//Sophie Germain is 2p+1, the other kind is 2p-1.
//Euler-Gandalf-Dipshit test
/*
	Form 2p+1 (Sophie Germain):
		n%8 == 7
		OR
		n%8 == 3
	Form 2p-1
		n%8 == 1
		OR
		n%8 == 5
*/

// disabling this, it's not any faster or more accurate than the Fermat test. It is indeed equivalent to the Fermat test.
/*
static
bool SpecialPrimeTest(mpz_t *n, bool fSophieGermain, unsigned int *pnLength)
{
	mpz_t e, r;
	mpz_init(e);
	mpz_sub_ui(e, *n, 1);
	mpz_fdiv_q_2exp(e, e, 1); //e = (n-1)/2;
	
	mpz_t a; mpz_init_set_ui(a, 2); //a = 2;

	mpz_init(r);
	mpz_powm(r, a, e, *n); //r = (2**((n-1)/2))%n;
	mpz_clear(a);
	mpz_clear(e);

	
	unsigned nMod8 = mpz_fdiv_ui(*n, 8);//nMod8 = n%8;
	bool fPassedTest = false;
	if (fSophieGermain && (nMod8 == 7)) // Euler & Lagrange
		fPassedTest = !mpz_cmp_ui(r, 1);
	else if (nMod8 == (fSophieGermain ? 3 : 5)) // Lifchitz
	{
		mpz_t mp;
		mpz_init_set_ui(mp, 1); //mp=1
		mpz_add(mp, r, mp); //mp = r+1;
		fPassedTest = !mpz_cmp(mp, *n); //r+1==n i.e. r%n == -1
		mpz_clear(mp);
	}
	else if ((!fSophieGermain) && (nMod8 == 1)) // LifChitz
		fPassedTest = !mpz_cmp_ui(r, 1);
	else
	{
		mpz_clear(r);
		cout << "ELLP Test: invalid n%%8 = " << nMod8 << ", " << (fSophieGermain?"First kind":"Second kind") << endl;
		return false;
		//return error("EulerLagrangeLifchitzPrimalityTest() : invalid n %% 8 = %d, %s", nMod8, (fSophieGermain? "first kind" : "second kind"));
	}

	if (fPassedTest)
	{
		mpz_clear(r);
		return true;
	}
	// Failed test, calculate fractional length
	
	// derive Fermat test remainder
	mpz_mul(r, r, r);
	mpz_fdiv_r(r, r, *n);
	
	// nFractionalLength = ( (n-r) << nFractionalBits ) / n
	mpz_sub(r, *n, r);
	mpz_mul_2exp(r, r, nFractionalBits);
	mpz_fdiv_q(r, r, *n);
	unsigned int nFractionalLength = mpz_get_ui(r);
	mpz_clear(r);
	
	if (nFractionalLength >= (1 << nFractionalBits))
	{
		cout << "EulerLagrangeLifchitzPrimalityTest() : fractional assert" << endl;
		return false;
		//return error("EulerLagrangeLifchitzPrimalityTest() : fractional assert");
	}
	*pnLength = (*pnLength & TARGET_LENGTH_MASK) | nFractionalLength;
	return false;
}
*/

// Test Probable Cunningham Chain for: n
// fSophieGermain:
//   true - Test for Cunningham Chain of first kind (n, 2n+1, 4n+3, ...)
//   false - Test for Cunningham Chain of second kind (n, 2n-1, 4n-3, ...)
// Return value:
//   true - Probable Cunningham Chain found (length at least 2)
//   false - Not Cunningham Chain
static
bool CunnChainTest(mpz_t *n, bool fSophieGermain, bool fFermatTest, unsigned int *pnProbableChainLength)
{
	*pnProbableChainLength = 0;
	mpz_t N;
	mpz_init_set(N, *n);
	
	// Fermat test for n first
	if (!PrimeTest(&N, pnProbableChainLength))
	{
		mpz_clear(N);
		return false;
	}

	// Euler-Lagrange-Lifchitz test for the following numbers in chain
	while (true)
	{
		TargetIncrementLength(pnProbableChainLength);
		mpz_add(N, N, N);
		if (fSophieGermain)
			mpz_add_ui(N, N, 1);
		else
			mpz_sub_ui(N, N, 1);
		//disabled, see SpecialPrimeTest() comments for further information
		/*if (fFermatTest)
		{
			if (!PrimeTest(&N, pnProbableChainLength))
				break;
		}
		else
		{
			if (!SpecialPrimeTest(&N, fSophieGermain, pnProbableChainLength))
				break;
		}*/
		if (!PrimeTest(&N, pnProbableChainLength))
			break;
	}
	mpz_clear(N);

#ifdef SUPERDEBUG
	printf("PCCT => %u (%u)\n", TargetGetLength(*pnProbableChainLength), *pnProbableChainLength);
#endif
	return (TargetGetLength(*pnProbableChainLength) >= 2);
}

// Test probable prime chain for: nOrigin
// Return value:
//   true - Probable prime chain found (one of nChainLength meeting target)
//   false - prime chain too short (none of nChainLength meeting target)
static
bool AllChainTest(mpz_t *bnPrimeChainOrigin, unsigned int nBits, bool fFermatTest, unsigned int *pnChainLengthCunningham1, unsigned int *pnChainLengthCunningham2, unsigned int *pnChainLengthBiTwin, unsigned int sievenumber)
{
	*pnChainLengthCunningham1 = 0;
	*pnChainLengthCunningham2 = 0;
	*pnChainLengthBiTwin = 0;
	
	mpz_t mp;
	mpz_init(mp);
	
	// Test for Cunningham Chain of first kind
	if (sievenumber&1)
	{
		mpz_sub_ui(mp, *bnPrimeChainOrigin, 1);
		CunnChainTest(&mp, true, fFermatTest, pnChainLengthCunningham1);
	}
	// Test for Cunningham Chain of second kind
	if (sievenumber&2)
	{
		mpz_add_ui(mp, *bnPrimeChainOrigin, 1);
		CunnChainTest(&mp, false, fFermatTest, pnChainLengthCunningham2);
	}
	
	mpz_clear(mp);
	// Figure out BiTwin Chain length
	// BiTwin Chain allows a single prime at the end for odd length chain
	*pnChainLengthBiTwin = (TargetGetLength(*pnChainLengthCunningham1) > TargetGetLength(*pnChainLengthCunningham2)) ? (*pnChainLengthCunningham2 + TargetFromInt(TargetGetLength(*pnChainLengthCunningham2)+1)) : (*pnChainLengthCunningham1 + TargetFromInt(TargetGetLength(*pnChainLengthCunningham1)));
	
	return (*pnChainLengthCunningham1 >= nBits || *pnChainLengthCunningham2 >= nBits || *pnChainLengthBiTwin >= nBits);
}

ullint chainspersec[20] = {};
ullint totalpersec = 0;

uint found0=0,foundtotal=0;

bool MinePrime(Reap_CPU_param* state, Work& tempwork)
{
	uchar* tempdata = &tempwork.data[0];
	uchar hash[32];
	mysha256(hash,tempdata,80);
	mysha256(hash,hash,32);
	
	//does this need byte flipping?
	uint bits = *(uint*)&tempdata[72];
	
	if (!(hash[31] & 0x80))
		return false; //hash is too small, abort

	Mpz_w hashnum;
	
	set_mpz_to_hash(&hashnum.n, hash);
	if (mpz_fdiv_ui(hashnum.n, 2*3) != 0)
		return false;
	
	bool found=false;

	
	//5431526412865007455
	mpz_t factor; mpz_init_set_str(factor, "5431526412865007455", 10);
	mpz_mul(hashnum.n,hashnum.n,factor);
	
	uint remainders[MAX_SIEVE_AMOUNT] = {};
	for(int i=0; i<MAX_SIEVE_AMOUNT; ++i)
	{
		remainders[i] = mpz_fdiv_ui(hashnum.n,Primes::v[i]);
	}
	
	mpz_t newhashnum; mpz_init(newhashnum);
	for(uint h=1; h<500; ++h)
	{
		mpz_add(newhashnum,newhashnum,hashnum.n);
		uint c1=0,c2=0,tw=0;

		uint sievenumber = 0;
		for(uint i=16; i<MAX_SIEVE_AMOUNT; ++i)
		{
			sievenumber |= Sieve::Get(i,remainders[i]*h%Primes::v[i])^3;
		}
		sievenumber ^= 3;
		//cout << sievenumber << endl;;
		if (sievenumber == 0)
			continue;
		
		//TODO: fix second parameter, it should be bits!
		AllChainTest(&newhashnum, 0, true, &c1, &c2, &tw, sievenumber);

		uint c1_i = TargetGetLength(c1);
		uint c2_i = TargetGetLength(c2);
		uint tw_i = TargetGetLength(tw);

		const int minlength=5;

		/*{
			if (c1_i >= minlength)
				cout << "First kind: " << c1_i << endl;
			if (c2_i >= minlength)
				cout << "Second kind: " << c2_i << endl;
			if (tw_i >= minlength)
				cout << "Twin kind: " << tw_i << endl;
		}*/
			
		
		++chainspersec[c1_i];
		++chainspersec[c2_i];
		++chainspersec[tw_i];
		++totalpersec;

		found = (c1_i >= minlength || c2_i >= minlength || tw_i >= minlength);
		if (found)
		{
			mpz_mul_ui(factor,factor,h);
			vector<uchar> auxdata = XPM_create_auxdata(&factor);
			Share share;
			
			CPU_Got_share(state,tempwork,auxdata);//tempdata,tempwork.target_share,current_server_id,tempwork.dataid,auxdata);
		}
	}

	return found;
}

void* Reap_CPU_XPM_mtrlt(void* param)
{

	Reap_CPU_param* state = (Reap_CPU_param*)param;

	Work tempwork;
	tempwork.time = 13371337;

	//uchar tempdata[80];
	//memset(tempdata, 0, 80);

	uchar finalhash[32];
	uchar temphash[32];
	uchar hash_results[1] = {};

	uint current_server_id;
	
	uint starttime = ticker();
	uint currenttime = starttime;
	
	uint foundprimes=0;

	while(!shutdown_now)
	{
		if (current_work.old)
		{
			Wait_ms(20);
			continue;
		}
		if (tempwork.time != current_work.time)
		{
			pthread_mutex_lock(&current_work_mutex);
			tempwork = current_work;
			pthread_mutex_unlock(&current_work_mutex);
			//memcpy(tempdata, &tempwork.data[0], 80);

			*(uint*)&tempwork.data[76] = state->thread_id<<28;
			current_server_id = tempwork.server_id;

		}


		uint trues=0;
		
		for(uint h=0; h<CPU_BATCH_SIZE; ++h)
		{
			bool result = MinePrime(state,tempwork);
			if (result) 
			{
				foundprimes++;
			}
			
			++*(uint*)&tempwork.data[76];
		}
		//cout << "Every " << double(CPU_BATCH_SIZE)/double(trues) << "th num is true" << endl;

		state->hashes += CPU_BATCH_SIZE;
	}
	pthread_exit(NULL);
	return NULL;
}


/*
TEMPLATE OF A CPU MINER
void* Reap_CPU_V1(void* param)
{
	Reap_CPU_param* state = (Reap_CPU_param*)param;

	Work tempwork;
	tempwork.time = 13371337;

	uchar tempdata[512];
	memset(tempdata, 0, 512);

	uchar finalhash[32];
	uchar hash_results[1] = {};

	uint current_server_id;

	while(!shutdown_now)
	{
		if (current_work.old)
		{
			Wait_ms(20);
			continue;
		}
		if (tempwork.time != current_work.time)
		{
			pthread_mutex_lock(&current_work_mutex);
			tempwork = current_work;
			pthread_mutex_unlock(&current_work_mutex);
			memcpy(tempdata, &tempwork.data[0], 128);
			*(uint*)&tempdata[100] = state->thread_id;
			current_server_id = tempwork.server_id;
		}

		*(ullint*)&tempdata[76] = tempwork.ntime_at_getwork + (ticker()-tempwork.time)/1000;

		for(uint h=0; h<CPU_BATCH_SIZE; ++h)
		{
			BlockHash_1_mine_V1(tempdata, finalhash, hash_results);
			if (hash_results[0])
			{
				BlockHash_1(tempdata, finalhash);
				if (finalhash[30] != 0 || finalhash[31] != 0)
					cpu_shares_hwinvalid++;
				else
					cpu_shares_hwvalid++;
				if (CPU_Hash_Below_Target(finalhash, &tempwork.target_share[0]))
					CPU_Got_share(state,tempdata,tempwork.target_share,current_server_id);
			}
			++*(uint*)&tempdata[108];
		}
		state->hashes += CPU_BATCH_SIZE;
	}
	pthread_exit(NULL);
	return NULL;
}
*/