#ifndef CSIEVEOFERATOSTHENES_H
#define CSIEVEOFERATOSTHENES_H

#include "gmpxx.h"
#include "CPUAlgos_global.h"

// Sieve of Eratosthenes for proof-of-work mining
class CSieveOfEratosthenes
{
    unsigned int nSieveSize; // size of the sieve
    unsigned int nBits; // target of the prime chain to search for
    mpz_class mpzFixedFactor; // fixed factor to derive the chain

    // final set of candidates for probable primality checking
    unsigned int *vfCandidates;
	
	unsigned int* vCunninghamMultipliers;
	unsigned int* vfCompositeCunningham;
    
    static const unsigned int nWordBits = 8 * sizeof(unsigned int);
    unsigned int nCandidatesWords;
    unsigned int nCandidatesBytes;

    unsigned int nPrimeSeq; // prime sequence number currently being processed
    unsigned int nCandidateCount; // cached total count of candidates
    unsigned int nCandidateMultiplier; // current candidate for power test
    
    unsigned int nChainLength;
    unsigned int nHalfChainLength;
    
    unsigned int GetWordNum(unsigned int nBitNum) {
        return nBitNum / nWordBits;
    }
    
    unsigned int GetBitMask(unsigned int nBitNum) {
        return 1UL << (nBitNum % nWordBits);
    }
    
    void AddMultiplier(unsigned int *vMultipliers, const unsigned int nSolvedMultiplier);

    void ProcessMultiplier(unsigned int *vfComposites, const unsigned int nPrime, unsigned int *vMultipliers)
    {
#ifdef USE_ROTATE
        const unsigned int nRotateBits = nPrime % nWordBits;
        for (unsigned int i = 0; i < nHalfChainLength; i++)
        {
            unsigned int nVariableMultiplier = vMultipliers[i];
            if (nVariableMultiplier == 0xFFFFFFFF) break;
            unsigned int lBitMask = GetBitMask(nVariableMultiplier);
            for (; nVariableMultiplier < nSieveSize; nVariableMultiplier += nPrime)
            {
                vfComposites[GetWordNum(nVariableMultiplier)] |= lBitMask;
                lBitMask = (lBitMask << nRotateBits) | (lBitMask >> (nWordBits - nRotateBits));
            }
            //vMultipliers[i] = nVariableMultiplier;
        }
#else
        for (unsigned int i = 0; i < nHalfChainLength; i++)
        {
            unsigned int nVariableMultiplier = vMultipliers[i];
            if (nVariableMultiplier == 0xFFFFFFFF) break;
            for (; nVariableMultiplier < nSieveSize; nVariableMultiplier += nPrime)
            {
                vfComposites[GetWordNum(nVariableMultiplier)] |= GetBitMask(nVariableMultiplier);
            }
            //vMultipliers[i] = nVariableMultiplier;
        }
#endif
    }

public:
	vector<uint> CandidateList;

	bool inited;
	CSieveOfEratosthenes()
	{
		inited = false;
	}
    ~CSieveOfEratosthenes()
    {
		Deinit();
    }

	void InitAndWeave(Reap_CPU_param* state, unsigned int nSieveSize, unsigned int nBits, mpz_class& mpzHash, mpz_class& mpzFixedMultiplier);
	void Deinit()
	{
		if (!inited)
			return;
		free(vfCandidates);
		vfCandidates = NULL;
		inited = false;
	}

    // Scan for the next candidate multiplier (variable part)
    // Return values:
    //   True - found next candidate; nVariableMultiplier has the candidate
    //   False - scan complete, no more candidate and reset scan
    bool GetNextCandidateMultiplier(unsigned int& nVariableMultiplier)
    {
        unsigned int lBits = vfCandidates[GetWordNum(nCandidateMultiplier)];
        while(true)
        {
            nCandidateMultiplier++;
            if (nCandidateMultiplier >= nSieveSize)
            {
                nCandidateMultiplier = 0;
                return false;
            }
            if (nCandidateMultiplier % nWordBits == 0)
            {
                lBits = vfCandidates[GetWordNum(nCandidateMultiplier)];
                if (lBits == 0)
                {
                    // Skip an entire word
                    nCandidateMultiplier += nWordBits - 1;
                    continue;
                }
            }
            if (lBits & GetBitMask(nCandidateMultiplier))
            {
                nVariableMultiplier = nCandidateMultiplier;
                return true;
            }
        }
    }

    void Weave();
	void GPUWeave(Reap_CPU_param* state);
};



#endif
