#include "Global.h"
#include "CPUMiner.h"
#include "CPUAlgos.h"
#include "Util.h"
#include "SHA256.h"
#include "CPUAlgos_global.h"

#include <iomanip>

ullint fermats=0,gandalfs=0;


// Copyright (c) 2013 Primecoin developers
// Distributed under conditional MIT/X11 software license,
// see the accompanying file COPYING

// Copyright (c) 2013 Primecoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

//START HEADER

#include <gmp.h>
#include <gmpxx.h>
#include <bitset>

static const unsigned int nMaxSievePercentage = 100;
static const unsigned int nDefaultSievePercentage = 10;
static const unsigned int nMinSievePercentage = 1;
extern unsigned int nSievePercentage;
static const unsigned int nMaxSieveSize = 1000000000u;
static const unsigned int nDefaultSieveSize = 1000000u;
static const unsigned int nMinSieveSize = 100000u;
extern unsigned int nSieveSize;
//static const uint256 hashBlockHeaderLimit = (uint256(1) << 255);
static const mpz_class mpzOne = 1;
static const mpz_class mpzTwo = 2;
static const mpz_class mpzPrimeMax = (mpzOne << 2000) - 1;
static const mpz_class mpzPrimeMin = (mpzOne << 255);

// Estimate how many 5-chains are found per hour
static const unsigned int nStatsChainLength = 5;

extern unsigned int nTargetInitialLength;
extern unsigned int nTargetMinLength;

// Generate small prime table
void GeneratePrimeTable();
// Get next prime number of p
bool PrimeTableGetNextPrime(unsigned int& p);
// Get previous prime number of p
bool PrimeTableGetPreviousPrime(unsigned int& p);

// Compute primorial number p#
void Primorial(unsigned int p, mpz_class& mpzPrimorial);
// Compute Primorial number p#
// Fast 32-bit version assuming that p <= 23
unsigned int PrimorialFast(unsigned int p);
// Compute the first primorial number greater than or equal to bn
void PrimorialAt(mpz_class& bn, mpz_class& mpzPrimorial);

// Test probable prime chain for: bnPrimeChainOrigin
// fFermatTest
//   true - Use only Fermat tests
//   false - Use Fermat-Euler-Lagrange-Lifchitz tests
// Return value:
//   true - Probable prime chain found (one of nChainLength meeting target)
//   false - prime chain too short (none of nChainLength meeting target)
//bool ProbablePrimeChainTest(const CBigNum& bnPrimeChainOrigin, unsigned int nBits, bool fFermatTest, unsigned int& nChainLengthCunningham1, unsigned int& nChainLengthCunningham2, unsigned int& nChainLengthBiTwin);

static const unsigned int nFractionalBits = 24;
static const unsigned int TARGET_FRACTIONAL_MASK = (1u<<nFractionalBits) - 1;
static const unsigned int TARGET_LENGTH_MASK = ~TARGET_FRACTIONAL_MASK;
static const uint64 nFractionalDifficultyMax = (1llu << (nFractionalBits + 32));
static const uint64 nFractionalDifficultyMin = (1llu << 32);
static const uint64 nFractionalDifficultyThreshold = (1llu << (8 + 32));
static const unsigned int nWorkTransitionRatio = 32;
unsigned int TargetGetLimit();
unsigned int TargetGetInitial();
unsigned int TargetGetLength(unsigned int nBits);
bool TargetSetLength(unsigned int nLength, unsigned int& nBits);
unsigned int TargetGetFractional(unsigned int nBits);
uint64 TargetGetFractionalDifficulty(unsigned int nBits);
bool TargetSetFractionalDifficulty(uint64 nFractionalDifficulty, unsigned int& nBits);
std::string TargetToString(unsigned int nBits);
unsigned int TargetFromInt(unsigned int nLength);

// Mine probable prime chain of form: n = h * p# +/- 1
bool MineProbablePrimeChain(Reap_CPU_param* state, Work& tempwork, mpz_class& mpzFixedMultiplier, bool& fNewBlock, unsigned int& nTriedMultiplier, unsigned int& nProbableChainLength, unsigned int& nTests, unsigned int& nPrimesHit, unsigned int& nChainsHit, mpz_class& mpzHash, unsigned int nPrimorialMultiplier);

// Check prime proof-of-work
enum // prime chain type
{
    PRIME_CHAIN_CUNNINGHAM1 = 1u,
    PRIME_CHAIN_CUNNINGHAM2 = 2u,
    PRIME_CHAIN_BI_TWIN     = 3u,
};

// prime chain type and length value
std::string GetPrimeChainName(unsigned int nChainType, unsigned int nChainLength);

#if defined(__i386__) || defined(_M_IX86) || defined(_X86_) || defined(__x86_64__) || defined(_M_X64)
#define USE_ROTATE
#endif

#include "CSieveOfEratosthenes.h"

// Prime Table
std::vector<unsigned int> vPrimes;
std::vector<unsigned int> vTwoInverses;
unsigned int nSieveSize = nDefaultSieveSize;
unsigned int nSievePercentage = nDefaultSievePercentage;

void GeneratePrimeTable()
{
    //nSievePercentage = (unsigned int)GetArg("-sievepercentage", nDefaultSievePercentage);
	nSievePercentage = globalconfs.coin.config.GetValue<uint>("sievepercentage");
    nSievePercentage = std::max(std::min(nSievePercentage, nMaxSievePercentage), nMinSievePercentage);

    //nSieveSize = (unsigned int)GetArg("-sievesize", nDefaultSieveSize);
	nSieveSize = globalconfs.coin.config.GetValue<uint>("sievesize");
    nSieveSize = std::max(std::min(nSieveSize, nMaxSieveSize), nMinSieveSize);

    printf("GeneratePrimeTable() : setting nSievePercentage = %u, nSieveSize = %u\n", nSievePercentage, nSieveSize);
    const unsigned nPrimeTableLimit = nSieveSize;
    vPrimes.clear();
    // Generate prime table using sieve of Eratosthenes
    std::vector<bool> vfComposite (nPrimeTableLimit, false);
    for (unsigned int nFactor = 2; nFactor * nFactor < nPrimeTableLimit; nFactor++)
    {
        if (vfComposite[nFactor])
            continue;
        for (unsigned int nComposite = nFactor * nFactor; nComposite < nPrimeTableLimit; nComposite += nFactor)
            vfComposite[nComposite] = true;
    }
    for (unsigned int n = 2; n < nPrimeTableLimit; n++)
        if (!vfComposite[n])
            vPrimes.push_back(n);
    printf("GeneratePrimeTable() : prime table [1, %d] generated with %lu primes", nPrimeTableLimit, vPrimes.size());
    //BOOST_FOREACH(unsigned int nPrime, vPrimes)
    //    printf(" %u", nPrime);
    printf("\n");
    
    const unsigned int nPrimes = vPrimes.size();
    vTwoInverses = std::vector<unsigned int> (nPrimes, 0);
    for (unsigned int nPrimeSeq = 1; nPrimeSeq < nPrimes; nPrimeSeq++)
    {
        vTwoInverses[nPrimeSeq] = int_invert(2, vPrimes[nPrimeSeq]);
    }
}

// Get next prime number of p
bool PrimeTableGetNextPrime(unsigned int& p)
{
    //BOOST_FOREACH(unsigned int nPrime, vPrimes)
	for(uint i=0; i<vPrimes.size(); ++i)
    {
		unsigned int nPrime = vPrimes[i];
        if (nPrime > p)
        {
            p = nPrime;
            return true;
        }
    }
    return false;
}

// Get previous prime number of p
bool PrimeTableGetPreviousPrime(unsigned int& p)
{
    unsigned int nPrevPrime = 0;
	//BOOST_FOREACH(unsigned int nPrime, vPrimes)
	for(uint i=0; i<vPrimes.size(); ++i)
    {
		unsigned int nPrime = vPrimes[i];
        if (nPrime >= p)
            break;
        nPrevPrime = nPrime;
    }
    if (nPrevPrime)
    {
        p = nPrevPrime;
        return true;
    }
    return false;
}

// Compute Primorial number p#
void Primorial(unsigned int p, mpz_class& mpzPrimorial)
{
    unsigned long nPrimorial = 1;
    unsigned int i;
    if (sizeof(unsigned long) >= 8)
    {
        // Fast 64-bit loop if possible
        for (i = 0; i < 15; i++)
        {
            unsigned int nPrime = vPrimes[i];
            if (nPrime > p) break;
            nPrimorial *= nPrime;
        }
    }
    else
    {
        // Fast 32-bit loop first
        for (i = 0; i < 9; i++)
        {
            unsigned int nPrime = vPrimes[i];
            if (nPrime > p) break;
            nPrimorial *= nPrime;
        }
    }

    mpzPrimorial = nPrimorial;
    for (; i < vPrimes.size(); i++)
    {
        unsigned int nPrime = vPrimes[i];
        if (nPrime > p) break;
        mpzPrimorial *= nPrime;
    }
}

// Compute Primorial number p#
// Fast 32-bit version assuming that p <= 23
unsigned int PrimorialFast(unsigned int p)
{
    unsigned int nPrimorial = 1;
	//BOOST_FOREACH(unsigned int nPrime, vPrimes)
	for(uint i=0; i<vPrimes.size(); ++i)
    {
		unsigned int nPrime = vPrimes[i];
        if (nPrime > p) break;
        nPrimorial *= nPrime;
    }
    return nPrimorial;
}

// Compute first primorial number greater than or equal to pn
void PrimorialAt(mpz_class& bn, mpz_class& mpzPrimorial)
{
	//BOOST_FOREACH(unsigned int nPrime, vPrimes)
	for(uint i=0; i<vPrimes.size(); ++i)
    {
		unsigned int nPrime = vPrimes[i];
        mpzPrimorial *= nPrime;
        if (mpzPrimorial >= bn)
            return;
    }
}

// Proof-of-work Target (prime chain target):
//   format - 32 bit, 8 length bits, 24 fractional length bits

unsigned int nTargetInitialLength = 7; // initial chain length target
unsigned int nTargetMinLength = 6;     // minimum chain length target

unsigned int TargetGetLimit()
{
    return (nTargetMinLength << nFractionalBits);
}

unsigned int TargetGetInitial()
{
    return (nTargetInitialLength << nFractionalBits);
}

unsigned int TargetGetLength(unsigned int nBits)
{
    return ((nBits & TARGET_LENGTH_MASK) >> nFractionalBits);
}

bool TargetSetLength(unsigned int nLength, unsigned int& nBits)
{
    if (nLength >= 0xff)
	{
		throw string("TargetSetLength() : invalid length=")+ToString(nLength);
        //return error("TargetSetLength() : invalid length=%u", nLength);
	}
    nBits &= TARGET_FRACTIONAL_MASK;
    nBits |= (nLength << nFractionalBits);
    return true;
}

static void TargetIncrementLength(unsigned int& nBits)
{
    nBits += (1 << nFractionalBits);
}

static void TargetDecrementLength(unsigned int& nBits)
{
    if (TargetGetLength(nBits) > nTargetMinLength)
        nBits -= (1 << nFractionalBits);
}

unsigned int TargetGetFractional(unsigned int nBits)
{
    return (nBits & TARGET_FRACTIONAL_MASK);
}

uint64 TargetGetFractionalDifficulty(unsigned int nBits)
{
    return (nFractionalDifficultyMax / (uint64) ((1llu<<nFractionalBits) - TargetGetFractional(nBits)));
}

bool TargetSetFractionalDifficulty(uint64 nFractionalDifficulty, unsigned int& nBits)
{
    if (nFractionalDifficulty < nFractionalDifficultyMin)
	{
		cout << "TargetSetFractionalDifficulty() : difficulty below min" << endl;
	}
    uint64 nFractional = nFractionalDifficultyMax / nFractionalDifficulty;
    if (nFractional > (1u<<nFractionalBits))
	{
		cout << "TargetSetFractionalDifficulty() : fractional overflow: nFractionalDifficulty=" << nFractionalDifficulty << endl;
	}
    nFractional = (1u<<nFractionalBits) - nFractional;
    nBits &= TARGET_LENGTH_MASK;
    nBits |= (unsigned int)nFractional;
    return true;
}

std::string TargetToString(unsigned int nBits)
{
    //return strprintf("%02x.%06x", TargetGetLength(nBits), TargetGetFractional(nBits));
	stringstream ss;
	ss << std::hex << std::setw(2) << std::setfill('0') << TargetGetLength(nBits) << "." << std::setw(6) << std::setfill('0') << TargetGetFractional(nBits);
	return ss.str();
}

unsigned int TargetFromInt(unsigned int nLength)
{
    return (nLength << nFractionalBits);
}

// Number of primes to test with fast divisibility testing
static const unsigned int nFastDivPrimes = 40;

class CPrimalityTestParams
{
public:
    // GMP variables
    mpz_t mpzN;
    mpz_t mpzE;
    mpz_t mpzR;
    mpz_t mpzRplusOne;
    
    // GMP C++ variables
    mpz_class mpzOriginMinusOne;
    mpz_class mpzOriginPlusOne;
    mpz_class N;

    // Values specific to a round
    unsigned int nBits;
    unsigned int nPrimorialSeq;

    // This is currently always false when mining
    static const bool fFermatTest = false;

    // Results
    unsigned int nChainLengthCunningham1;
    unsigned int nChainLengthCunningham2;
    unsigned int nChainLengthBiTwin;

    CPrimalityTestParams(unsigned int nBits, unsigned int nPrimorialSeq)
    {
        this->nBits = nBits;
        this->nPrimorialSeq = nPrimorialSeq;
        nChainLengthCunningham1 = 0;
        nChainLengthCunningham2 = 0;
        nChainLengthBiTwin = 0;
        mpz_init(mpzN);
        mpz_init(mpzE);
        mpz_init(mpzR);
        mpz_init(mpzRplusOne);
    }

    ~CPrimalityTestParams()
    {
        mpz_clear(mpzN);
        mpz_clear(mpzE);
        mpz_clear(mpzR);
        mpz_clear(mpzRplusOne);
    }
};

// Check Fermat probable primality test (2-PRP): 2 ** (n-1) = 1 (mod n)
// true: n is probable prime
// false: n is composite; set fractional length in the nLength output
static bool FermatProbablePrimalityTestFast(const mpz_class& n, unsigned int& nLength, CPrimalityTestParams& testParams, bool fFastDiv = false)
{
    // Faster GMP version
    mpz_t& mpzN = testParams.mpzN;
    mpz_t& mpzE = testParams.mpzE;
    mpz_t& mpzR = testParams.mpzR;
    const unsigned int nPrimorialSeq = testParams.nPrimorialSeq;

    mpz_set(mpzN, n.get_mpz_t());
    if (fFastDiv)
    {
        // Fast divisibility tests
        // Starting from the first prime not included in the round primorial
        const unsigned int nBeginSeq = nPrimorialSeq + 1;
        const unsigned int nEndSeq = nBeginSeq + nFastDivPrimes;
        for (unsigned int nPrimeSeq = nBeginSeq; nPrimeSeq < nEndSeq; nPrimeSeq++) {
            if (mpz_divisible_ui_p(mpzN, vPrimes[nPrimeSeq])) {
               return false;
            }
        }
    }

	++fermats;
    mpz_sub_ui(mpzE, mpzN, 1);
    mpz_powm(mpzR, mpzTwo.get_mpz_t(), mpzE, mpzN);
    if (mpz_cmp_ui(mpzR, 1) == 0)
    {
        return true;
    }
    // Failed Fermat test, calculate fractional length
    mpz_sub(mpzE, mpzN, mpzR);
    mpz_mul_2exp(mpzR, mpzE, nFractionalBits);
    mpz_tdiv_q(mpzE, mpzR, mpzN);
    unsigned int nFractionalLength = mpz_get_ui(mpzE);

    if (nFractionalLength >= (1 << nFractionalBits))
	{
		cout << "FermatProbablePrimalityTest() : fractional assert" << endl;
        return false;
	}
    nLength = (nLength & TARGET_LENGTH_MASK) | nFractionalLength;
    return false;
}

// Test probable primality of n = 2p +/- 1 based on Euler, Lagrange and Lifchitz
// fSophieGermain:
//   true:  n = 2p+1, p prime, aka Cunningham Chain of first kind
//   false: n = 2p-1, p prime, aka Cunningham Chain of second kind
// Return values
//   true: n is probable prime
//   false: n is composite; set fractional length in the nLength output
static bool EulerLagrangeLifchitzPrimalityTestFast(const mpz_class& n, bool fSophieGermain, unsigned int& nLength, CPrimalityTestParams& testParams, bool fFastDiv = false)
{
    // Faster GMP version
    mpz_t& mpzN = testParams.mpzN;
    mpz_t& mpzE = testParams.mpzE;
    mpz_t& mpzR = testParams.mpzR;
    mpz_t& mpzRplusOne = testParams.mpzRplusOne;
    const unsigned int nPrimorialSeq = testParams.nPrimorialSeq;

    mpz_set(mpzN, n.get_mpz_t());
    if (fFastDiv)
    {
        // Fast divisibility tests
        // Starting from the first prime not included in the round primorial
        const unsigned int nBeginSeq = nPrimorialSeq + 1;
        const unsigned int nEndSeq = nBeginSeq + nFastDivPrimes;
        for (unsigned int nPrimeSeq = nBeginSeq; nPrimeSeq < nEndSeq; nPrimeSeq++) {
            if (mpz_divisible_ui_p(mpzN, vPrimes[nPrimeSeq])) {
                return false;
            }
        }
    }
	++gandalfs;
    mpz_sub_ui(mpzE, mpzN, 1);
    mpz_tdiv_q_2exp(mpzE, mpzE, 1);
    mpz_powm(mpzR, mpzTwo.get_mpz_t(), mpzE, mpzN);
    unsigned int nMod8 = mpz_tdiv_ui(mpzN, 8);
    bool fPassedTest = false;
    if (fSophieGermain && (nMod8 == 7)) // Euler & Lagrange
        fPassedTest = !mpz_cmp_ui(mpzR, 1);
    else if (fSophieGermain && (nMod8 == 3)) // Lifchitz
    {
        mpz_add_ui(mpzRplusOne, mpzR, 1);
        fPassedTest = !mpz_cmp(mpzRplusOne, mpzN);
    }
    else if ((!fSophieGermain) && (nMod8 == 5)) // Lifchitz
    {
        mpz_add_ui(mpzRplusOne, mpzR, 1);
        fPassedTest = !mpz_cmp(mpzRplusOne, mpzN);
    }
    else if ((!fSophieGermain) && (nMod8 == 1)) // LifChitz
        fPassedTest = !mpz_cmp_ui(mpzR, 1);
    else
	{
		cout << "EulerLagrangeLifchitzPrimalityTest() : invalid n %% 8 = " << nMod8 << ", " << (fSophieGermain? "first kind" : "second kind") << endl;
		return false;
        //return error("EulerLagrangeLifchitzPrimalityTest() : invalid n %% 8 = %d, %s", nMod8, (fSophieGermain? "first kind" : "second kind"));
	}
    
    if (fPassedTest)
    {
        return true;
    }
    
    // Failed test, calculate fractional length
    mpz_mul(mpzE, mpzR, mpzR);
    mpz_tdiv_r(mpzR, mpzE, mpzN); // derive Fermat test remainder

    mpz_sub(mpzE, mpzN, mpzR);
    mpz_mul_2exp(mpzR, mpzE, nFractionalBits);
    mpz_tdiv_q(mpzE, mpzR, mpzN);
    unsigned int nFractionalLength = mpz_get_ui(mpzE);
    
    if (nFractionalLength >= (1 << nFractionalBits))
	{
		cout << "EulerLagrangeLifchitzPrimalityTest() : fractional assert" << endl;
        return false;
	}
    nLength = (nLength & TARGET_LENGTH_MASK) | nFractionalLength;
    return false;
}

// Test Probable Cunningham Chain for: n
// fSophieGermain:
//   true - Test for Cunningham Chain of first kind (n, 2n+1, 4n+3, ...)
//   false - Test for Cunningham Chain of second kind (n, 2n-1, 4n-3, ...)
// Return value:
//   true - Probable Cunningham Chain found (length at least 2)
//   false - Not Cunningham Chain
static bool ProbableCunninghamChainTestFast(const mpz_class& n, bool fSophieGermain, bool fFermatTest, unsigned int& nProbableChainLength, CPrimalityTestParams& testParams, bool use_gpu_fermat_test)
{
    nProbableChainLength = 0;
    mpz_class &N = testParams.N;
    N = n;

	if (!use_gpu_fermat_test && !FermatProbablePrimalityTestFast(N, nProbableChainLength, testParams, true))
	{
		return false;
	}

    // Euler-Lagrange-Lifchitz test for the following numbers in chain
    while (true)
    {
        TargetIncrementLength(nProbableChainLength);
        N = N + N + (fSophieGermain? 1 : (-1));
        if (fFermatTest)
        {
            if (!FermatProbablePrimalityTestFast(N, nProbableChainLength, testParams, true))
                break;
        }
        else
        {
            if (!EulerLagrangeLifchitzPrimalityTestFast(N, fSophieGermain, nProbableChainLength, testParams, true))
                break;
        }
    }
    return (TargetGetLength(nProbableChainLength) >= 2);
}

// Test probable prime chain for: nOrigin
// Return value:
//   true - Probable prime chain found (one of nChainLength meeting target)
//   false - prime chain too short (none of nChainLength meeting target)
static bool ProbablePrimeChainTestFast(const mpz_class& mpzPrimeChainOrigin, CPrimalityTestParams& testParams, uint sievenumber, bool use_gpu_fermat_test)
{
    const unsigned int nBits = testParams.nBits;
	const unsigned int nBits_masked = nBits&TARGET_LENGTH_MASK;
    unsigned int& nChainLengthCunningham1 = testParams.nChainLengthCunningham1;
    unsigned int& nChainLengthCunningham2 = testParams.nChainLengthCunningham2;
    unsigned int& nChainLengthBiTwin = testParams.nChainLengthBiTwin;
    const bool fFermatTest = testParams.fFermatTest;
    mpz_class& mpzOriginMinusOne = testParams.mpzOriginMinusOne;
    mpz_class& mpzOriginPlusOne = testParams.mpzOriginPlusOne;
    nChainLengthCunningham1 = 0;
    nChainLengthCunningham2 = 0;
    nChainLengthBiTwin = 0;
	
    // Test for Cunningham Chain of first kind
	if (sievenumber&1)
	{
		mpzOriginMinusOne = mpzPrimeChainOrigin - 1;
		ProbableCunninghamChainTestFast(mpzOriginMinusOne, true, fFermatTest, nChainLengthCunningham1, testParams, use_gpu_fermat_test);
		if ((nChainLengthCunningham1&TARGET_FRACTIONAL_MASK) == 0)
			nChainLengthCunningham1 |= TARGET_FRACTIONAL_MASK;
	}
    // Test for Cunningham Chain of second kind
	if (sievenumber&2)
	{
		mpzOriginPlusOne = mpzPrimeChainOrigin + 1;
		ProbableCunninghamChainTestFast(mpzOriginPlusOne, false, fFermatTest, nChainLengthCunningham2, testParams, use_gpu_fermat_test);
		if ((nChainLengthCunningham2&TARGET_FRACTIONAL_MASK) == 0)
			nChainLengthCunningham2 |= TARGET_FRACTIONAL_MASK;
	}
	// Figure out BiTwin Chain length
    // BiTwin Chain allows a single prime at the end for odd length chain
    nChainLengthBiTwin =
        (TargetGetLength(nChainLengthCunningham1) > TargetGetLength(nChainLengthCunningham2))?
            (nChainLengthCunningham2 + TargetFromInt(TargetGetLength(nChainLengthCunningham2)+1)) :
            (nChainLengthCunningham1 + TargetFromInt(TargetGetLength(nChainLengthCunningham1)));
			
	uint c1 = TargetGetLength(nChainLengthCunningham1);
	uint c2 = TargetGetLength(nChainLengthCunningham2);
	uint tw = TargetGetLength(nChainLengthBiTwin);
	
	if (c1 >= 6)
	{
		cout << "C1 " << nChainLengthCunningham1 << " --> " << TargetToString(nChainLengthCunningham1) << " found!" << endl;
	}
	if (c2 >= 6)
	{
		cout << "C2 " << nChainLengthCunningham2 << " --> " << TargetToString(nChainLengthCunningham2) << " found!" << endl;
	}
	if (tw >= 6)
	{
		cout << "TW " << nChainLengthBiTwin << " --> " << TargetToString(nChainLengthBiTwin) << " found!" << endl;
	}
	
	++chainspersec[c1];
	++chainspersec[c2];
	++chainspersec[tw];
	++totalpersec;
	
    return (nChainLengthCunningham1 >= nBits || nChainLengthCunningham2 >= nBits || nChainLengthBiTwin >= nBits);
}

// Mine probable prime chain of form: n = h * p# +/- 1
bool MineProbablePrimeChain(Reap_CPU_param* state, Work& tempwork, CSieveOfEratosthenes& psieve, mpz_class& mpzFixedMultiplier, bool& fNewBlock, unsigned int& nTriedMultiplier, unsigned int& nProbableChainLength, unsigned int& nTests, unsigned int& nPrimesHit, unsigned int& nChainsHit, mpz_class& mpzHash, unsigned int nPrimorialMultiplier)
{
    nProbableChainLength = 0;
    nPrimesHit = 0;
    nChainsHit = 0;
    //const unsigned int nBits = block.nBits;
	const unsigned int nBits = *(uint*)&tempwork.data[72];
	
	bool use_gpu_fermat_test = globalconfs.coin.config.GetValue<bool>("use_gpu_fermat_test");

    if (fNewBlock && psieve.inited)
    {
        // Must rebuild the sieve
		psieve.Deinit();
    }
    fNewBlock = false;

    int64 nStart; // microsecond timer
    if (!psieve.inited)
    {
        // Build sieve
        nStart = ticker()*1000;
		psieve.InitAndWeave(state, nSieveSize, nBits, mpzHash, mpzFixedMultiplier);
        if (globalconfs.coin.config.GetValue<bool>("debug"))
            printf("MineProbablePrimeChain() : new sieve (%u/%u) ready in %uus\n", psieve.CandidateList.size(), nSieveSize, (unsigned int) (ticker()*1000 - nStart));
    }

    mpz_class mpzHashMultiplier = mpzHash * mpzFixedMultiplier;
    mpz_class mpzChainOrigin;

    // Determine the sequence number of the round primorial
    unsigned int nPrimorialSeq = 0;
    while (vPrimes[nPrimorialSeq + 1] <= nPrimorialMultiplier)
        nPrimorialSeq++;

    // Allocate GMP variables for primality tests
    CPrimalityTestParams testParams(nBits, nPrimorialSeq);
    nStart = ticker()*1000;

    // References to counters;
    unsigned int& nChainLengthCunningham1 = testParams.nChainLengthCunningham1;
    unsigned int& nChainLengthCunningham2 = testParams.nChainLengthCunningham2;
    unsigned int& nChainLengthBiTwin = testParams.nChainLengthBiTwin;
	
	//cout << "PSIEVIOSIE" << psieve.CandidateList.size() << endl;
	
	for(uint i=0; i<psieve.CandidateList.size(); ++i)
    {
		nTriedMultiplier = psieve.CandidateList[i]&0x3FFFFFFFU;
		uint sievenumber = psieve.CandidateList[i]>>30;
		if (sievenumber == 0)
			sievenumber=3;
		if (nTriedMultiplier == 0) //will crash otherwise
			continue;
		++nTests;
		if (tempwork.time != current_work.time)
		{
			//cout << "Tempwork.time != curnetopqi" << tempwork.time << " " << current_work.time << endl;
			break;
		}
        mpzChainOrigin = mpzHashMultiplier * (nTriedMultiplier&0x3FFFFFFFU);
        nChainLengthCunningham1 = 0;
        nChainLengthCunningham2 = 0;
        nChainLengthBiTwin = 0;
        if (ProbablePrimeChainTestFast(mpzChainOrigin, testParams, sievenumber, use_gpu_fermat_test))
        {
			mpz_t mpzPrimeChainMultiplier; mpz_init(mpzPrimeChainMultiplier);
			mpz_mul_ui(mpzPrimeChainMultiplier,mpzFixedMultiplier.get_mpz_t(),nTriedMultiplier);
			{
				//gmp_printf("Found chain! Mult: %Zx\n",mpzPrimeChainMultiplier);
				vector<uchar> auxdata = XPM_create_auxdata(&mpzPrimeChainMultiplier);
				CPU_Got_share(state,tempwork,auxdata);
			}
			mpz_clear(mpzPrimeChainMultiplier);

            nProbableChainLength = std::max(std::max(nChainLengthCunningham1, nChainLengthCunningham2), nChainLengthBiTwin);
            return true;
        }
        nProbableChainLength = std::max(std::max(nChainLengthCunningham1, nChainLengthCunningham2), nChainLengthBiTwin);
        if(TargetGetLength(nProbableChainLength) >= 1)
            nPrimesHit++;
        if(TargetGetLength(nProbableChainLength) >= nStatsChainLength)
            nChainsHit++;
    }

	// power tests completed for the sieve
	//if (fDebug && GetBoolArg("-printmining"))
		//printf("MineProbablePrimeChain() : %u tests (%u primes and %u %d-chains) in %uus\n", nTests, nPrimesHit, nChainsHit, nStatsChainLength, (unsigned int) (GetTimeMicros() - nStart));
	psieve.Deinit();
	fNewBlock = true; // notify caller to change nonce
	return false; // stop as new block arrived
}

// prime chain type and length value
std::string GetPrimeChainName(unsigned int nChainType, unsigned int nChainLength)
{
	return (nChainType==PRIME_CHAIN_CUNNINGHAM1)? "1CC" : ((nChainType==PRIME_CHAIN_CUNNINGHAM2)? "2CC" : "TWN") + TargetToString(nChainLength);
}

unsigned int int_invert(unsigned int a, unsigned int nPrime)
{
    // Extended Euclidean algorithm to calculate the inverse of a in finite field defined by nPrime
    int rem0 = nPrime, rem1 = a % nPrime, rem2;
    int aux0 = 0, aux1 = 1, aux2;
    int quotient, inverse;

    while (1)
    {
        if (rem1 <= 1)
        {
            inverse = aux1;
            break;
        }

        rem2 = rem0 % rem1;
        quotient = rem0 / rem1;
        aux2 = -quotient * aux1 + aux0;

        if (rem2 <= 1)
        {
            inverse = aux2;
            break;
        }

        rem0 = rem1 % rem2;
        quotient = rem1 / rem2;
        aux0 = -quotient * aux2 + aux1;

        if (rem0 <= 1)
        {
            inverse = aux0;
            break;
        }

        rem1 = rem2 % rem0;
        quotient = rem2 / rem0;
        aux1 = -quotient * aux0 + aux2;
    }

    return (inverse + nPrime) % nPrime;
}

bool MinePrime_hp(Reap_CPU_param* state, Work& tempwork)
{
	uchar* tempdata = &tempwork.data[0];
	uchar hash[32];
	mysha256(hash,tempdata,80);
	mysha256(hash,hash,32);
	
	uint nBits = *(uint*)&tempdata[72];
	
	if (!(hash[31] & 0x80))
		return false; //hash is too small, abort

	Mpz_w hashnum;
	
	set_mpz_to_hash(&hashnum.n, hash);
	
	bool found = false;
	
	//START HP7 CODE
	CSieveOfEratosthenes psieve;
	static const unsigned int nPrimorialHashFactor = 7;
    unsigned int nPrimorialMultiplier = nPrimorialHashFactor;	
	unsigned int nHashFactor = PrimorialFast(nPrimorialHashFactor);
    int64 nTimeExpected = 0;   // time expected to prime chain (micro-second)
    int64 nTimeExpectedPrev = 0; // time expected to prime chain last time
    bool fIncrementPrimorial = true; // increase or decrease primorial factor

	while(true)
	{
		if (tempwork.time != current_work.time)
		{
			//cout << "OUTER" << endl;
			break;
		}
        int64 nStart = ticker()/1000;
        bool fNewBlock = true;
        unsigned int nTriedMultiplier = 0;

        // Primecoin: try to find hash divisible by primorial
        unsigned int nHashFactor = PrimorialFast(nPrimorialHashFactor);

        // Based on mustyoshi's patch from https://bitcointalk.org/index.php?topic=251850.msg2689981#msg2689981
        mpz_class mpzHash;
		mpzHash = mpz_class(hashnum.n);
        if (!mpz_divisible_ui_p(mpzHash.get_mpz_t(), nHashFactor))
			return false;

        // Use the hash that passed the tests
        // Primecoin: primorial fixed multiplier
        mpz_class mpzPrimorial;
        unsigned int nRoundTests = 0;
        unsigned int nRoundPrimesHit = 0;
        int64 nPrimeTimerStart = ticker()*1000;
        Primorial(nPrimorialMultiplier, mpzPrimorial);

        while(true)
        {
			if (tempwork.time != current_work.time)
			{
				//cout << "INNER" << endl;
				break;
			}

            unsigned int nTests = 0;
            unsigned int nPrimesHit = 0;
            unsigned int nChainsHit = 0;

            // Primecoin: adjust round primorial so that the generated prime candidates meet the minimum
            mpz_class mpzFixedMultiplier;
			
            mpzFixedMultiplier = 11*13*17*19*23*29*31;
			mpzFixedMultiplier *= 37*41*43*47*53;
			mpzFixedMultiplier *= 59;
			mpzFixedMultiplier *= 61*67*71*73*79;

            // Primecoin: mine for prime chain
            unsigned int nProbableChainLength;
            if (MineProbablePrimeChain(state, tempwork, psieve, mpzFixedMultiplier, fNewBlock, nTriedMultiplier, nProbableChainLength, nTests, nPrimesHit, nChainsHit, mpzHash, nPrimorialMultiplier))
				return true;

			nRoundTests += nTests;
            nRoundPrimesHit += nPrimesHit;

            if (fNewBlock)
				return false;
        }
	}
	
	//END HP7 CODE
	return found;
}

void* Reap_CPU_XPM_hp7(void* param)
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

			*(uint*)&tempwork.data[76] = state->thread_id<<28;
			current_server_id = tempwork.server_id;
		}

		bool result = MinePrime_hp(state,tempwork);
		if (result) 
		{
			pthread_mutex_lock(&current_work_mutex);
			current_work.old = true;
			pthread_mutex_unlock(&current_work_mutex);
		}
		++*(uint*)&tempwork.data[76];
	}
	pthread_exit(NULL);
	return NULL;
}
