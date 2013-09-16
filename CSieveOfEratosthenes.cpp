#include "Global.h"
#include "CSieveOfEratosthenes.h"
#include "CPUAlgos_global.h"
#include "CPUAlgos.h"

void CSieveOfEratosthenes::AddMultiplier(unsigned int *vMultipliers, const unsigned int nSolvedMultiplier)
{
    // Eliminate duplicates
    for (unsigned int i = 0; i < nHalfChainLength; i++)
    {
        unsigned int& n = vMultipliers[i];
        if (n == 0xFFFFFFFF || n == nSolvedMultiplier)
        {
            n = nSolvedMultiplier;
            break;
        }
    }
}

void CSieveOfEratosthenes::InitAndWeave(Reap_CPU_param* state, unsigned int nSieveSize, unsigned int nBits, mpz_class& mpzHash, mpz_class& mpzFixedMultiplier)
{
	if (inited)
		return;

	this->nSieveSize = nSieveSize;
	if (nSieveSize%nWordBits != 0)
		cout << "Sieve size must be a multiple of " << nWordBits << endl;
	this->nBits = nBits;
	this->mpzFixedFactor = mpzFixedMultiplier * mpzHash;
	nPrimeSeq = 0;
	nCandidateCount = 0;
	nCandidateMultiplier = 0;
	nCandidatesWords = nSieveSize / nWordBits;
	nCandidatesBytes = nCandidatesWords * sizeof(unsigned int);
	vfCandidates = (unsigned int *)malloc(nCandidatesBytes);
	memset(vfCandidates, 0, nCandidatesBytes);
	
	inited = true;
	// Faster GMP version
	this->nChainLength = TargetGetLength(nBits);
	this->nHalfChainLength = (this->nChainLength + 1) / 2;
	
    const unsigned int nPrimes = (uint64)vPrimes.size() * nSievePercentage / 100;
	
	const unsigned int zSIZE = 4;
	const unsigned int ySIZE = nPrimes;
	const unsigned int xSIZE = nHalfChainLength;
	
	vCunninghamMultipliers = new unsigned int[zSIZE*ySIZE*xSIZE];//order: 1A, 2A, 1B, 2B
    memset(vCunninghamMultipliers, 0xFF, zSIZE*ySIZE*xSIZE*sizeof(unsigned int));

	vfCompositeCunningham = new unsigned int[4*nCandidatesWords];
    memset(vfCompositeCunningham, 0, 4*nCandidatesBytes);
	
	if (globalconfs.coin.config.GetValue<bool>("use_gpu"))
		GPUWeave(state);
	else
		Weave();

	delete [] vfCompositeCunningham;
	delete [] vCunninghamMultipliers;
}

void CSieveOfEratosthenes::GPUWeave(Reap_CPU_param* state)
{
	while(!app.opencl.inited);

	//uint devicenum=rand()%GPUstates.size(); //FUCK YOU ALL
	uint devicenum=state->thread_id%GPUstates.size();
	pthread_mutex_lock(&GPUstates[devicenum].device_mutex);
	uint starttime=ticker();
	
	const unsigned int nPrimes = (uint64)vPrimes.size() * nSievePercentage / 100;
	const unsigned int zSIZE = 4;
	const unsigned int ySIZE = nPrimes;
	const unsigned int xSIZE = nHalfChainLength;
	
	//this parameter hardcoded for AMD GPUs
	const uint localmemsize=4096;
	
	map<string,string> defines;
	defines["NCHAINLENGTH"] = ToString(nChainLength);
	defines["LOCAL_MEM_USED"] = ToString(localmemsize);
	app.opencl.RebuildKernels(GPUstates[devicenum],defines);
	
	app.opencl.WriteBufferPattern(devicenum,"vfCompositeCunningham", 4*nCandidatesBytes, vfCompositeCunningham,16);
	app.opencl.WriteBufferPattern(devicenum,"vfCunninghamMultipliers", zSIZE*ySIZE*xSIZE*sizeof(unsigned int), vCunninghamMultipliers,16);
	
	uint mpzFixedFactor_uint16[16] = {};
	size_t factorsize=0;
	mpz_export(mpzFixedFactor_uint16, &factorsize, -1, 4, -1, 0, mpzFixedFactor.get_mpz_t());
	if (globalconfs.coin.config.GetValue<bool>("opencldebug"))
		cout << "Exported size " << factorsize << " words" << endl;
	if (factorsize > 11)
		cout << "Factor size is too big! Kernel cannot handle yet :(" << endl;
	
	//__kernel void CalculateMultipliers(__global const uint*restrict vPrimes, const uint nPrimes, const uint16 nFixedFactor, const uint nCandidatesWords, __global uint* vfCunninghamMultipliers)	
	vector<std::pair<string,string> > args;
	args.push_back(std::pair<string,string>("buffer","vPrimes"));
	args.push_back(std::pair<string,string>("uint",string((char*)&nPrimes,sizeof(uint))));
	args.push_back(std::pair<string,string>("uint16",string((char*)&mpzFixedFactor_uint16,16*sizeof(uint))));
	args.push_back(std::pair<string,string>("uint",string((char*)&nCandidatesWords,sizeof(uint))));
	args.push_back(std::pair<string,string>("buffer","vfCunninghamMultipliers"));
	//nPrimes-1 because the first prime isnt used. i hate you, only even prime.
	app.opencl.RunKernel(devicenum, "CalculateMultipliers",nPrimes-1,args);

	args.clear();
	//__kernel void Sieve(__global uint* vfCompositeCunningham, __global const uint* vfCunninghamMultipliers, __global const uint* vPrimes, const uint nPrimes, const uint nCandidatesWords)
	args.push_back(std::pair<string,string>("buffer","vfCompositeCunningham"));
	args.push_back(std::pair<string,string>("buffer","vfCunninghamMultipliers"));
	args.push_back(std::pair<string,string>("buffer","vPrimes"));
	args.push_back(std::pair<string,string>("uint",string((char*)&nPrimes,sizeof(uint))));
	args.push_back(std::pair<string,string>("uint",string((char*)&nCandidatesWords,sizeof(uint))));
	{
		//these parameters hardcoded for AMD GPUs
		uint localsize=globalconfs.coin.config.GetValue<uint>("worksize");
		uint sieveperblock=localmemsize*8;
		uint blocks = (nSieveSize+sieveperblock-1)/sieveperblock;
		uint threads = blocks*localsize;
		app.opencl.RunKernel(devicenum, "Sieve",threads,args);
	}

	uint counter=0;
	app.opencl.WriteBuffer(devicenum,"Counter",&counter,4);

	//__kernel void Combine(__global const uint* vfCompositeCunningham, const uint nCandidatesWords, __global uint*restrict counter, __global uint*restrict CandidateList)
	args.clear();
	args.push_back(std::pair<string,string>("buffer","vfCompositeCunningham"));
	args.push_back(std::pair<string,string>("uint",string((char*)&nCandidatesWords,sizeof(uint))));
	args.push_back(std::pair<string,string>("buffer","Counter"));
	args.push_back(std::pair<string,string>("buffer","CandidateList"));
	app.opencl.RunKernel(devicenum, "Combine",nCandidatesWords,args);

	app.opencl.ReadBuffer(devicenum,"Counter",&counter, 4);
	CandidateList.assign(counter,0);
	app.opencl.ReadBuffer(devicenum,"CandidateList",&CandidateList[0],counter*sizeof(uint));

	
	if (globalconfs.coin.config.GetValue<bool>("use_gpu_fermat_test"))
	{
		//__kernel void Fermat(__global const uint* FermatNumbers, __global uint* FermatOutput, uint16 FixedMultiplier)
		args.clear();
		args.push_back(std::pair<string,string>("buffer","CandidateList"));
		args.push_back(std::pair<string,string>("buffer","FermatOutput"));
		args.push_back(std::pair<string,string>("uint16",string((char*)&mpzFixedFactor_uint16,16*sizeof(uint))));
		
		app.opencl.RunKernel(devicenum,"Fermat",counter*2,args);

		vector<uint> FermatOutput;
		FermatOutput.assign(counter*2,0);
		app.opencl.ReadBuffer(devicenum,"FermatOutput",&FermatOutput[0],counter*2*sizeof(uint));

		vector<uint> TrimmedCandidateList;
		for(uint i=0; i<counter-counter%32; ++i)
		{
			uint candidatenum = 0;
			if (FermatOutput[2*i+0] == 1)
				candidatenum |= CandidateList[i]|0x40000000;
			if (FermatOutput[2*i+1] == 1)
				candidatenum |= CandidateList[i]|0x80000000;
				
			if (candidatenum)
				TrimmedCandidateList.push_back(candidatenum);
		}
		CandidateList.swap(TrimmedCandidateList);
		
		fermats += (counter-counter%32)*2;
	}
	GPUstates[devicenum].in_use=false;
	uint stoptime=ticker();
	pthread_mutex_unlock(&GPUstates[devicenum].device_mutex);
	//cout << stoptime-starttime << "ms for GPU stuff! " << CandidateList.size() << " candidates found!" << endl;
}

void CSieveOfEratosthenes::Weave()
{
	uint starttime=ticker();
    // Process only a set percentage of the primes
    // Most composites are still found
    const unsigned int nPrimes = (uint64)vPrimes.size() * nSievePercentage / 100;
	const unsigned int zSIZE = 4;
	const unsigned int ySIZE = nPrimes;
	const unsigned int xSIZE = nHalfChainLength;

    for (unsigned int nPrimeSeq = 1; nPrimeSeq < nPrimes; nPrimeSeq++)
    {
        unsigned int nPrime = vPrimes[nPrimeSeq];
        unsigned int nFixedFactorMod = mpz_tdiv_ui(mpzFixedFactor.get_mpz_t(), nPrime);
		unsigned int nFixedInverse,nTwoInverse;
        if (nFixedFactorMod == 0)
			continue;
		nFixedInverse = int_invert(nFixedFactorMod, nPrime);
		nTwoInverse = vTwoInverses[nPrimeSeq];

		// Weave the sieve for the prime
		for (unsigned int nBiTwinSeq = 0; nBiTwinSeq < 2 * nChainLength; nBiTwinSeq+=2)
		{
			unsigned int arrnum;
			unsigned int nSolvedMultiplier;
			// Find the first number that's divisible by this prime
			arrnum = 2*(nBiTwinSeq+1 >= nChainLength);
			nSolvedMultiplier = (nPrime-nFixedInverse)%nPrime;
			AddMultiplier(&vCunninghamMultipliers[(arrnum+1)*ySIZE*xSIZE+nPrimeSeq*xSIZE], nSolvedMultiplier);

			arrnum = 2*(nBiTwinSeq >= nChainLength);
			nSolvedMultiplier = nFixedInverse;
			AddMultiplier(&vCunninghamMultipliers[arrnum*ySIZE*xSIZE+nPrimeSeq*xSIZE], nSolvedMultiplier);
			
			nFixedInverse = (uint64)nFixedInverse * nTwoInverse % nPrime;
		}
		for(uint i=0; i<4; ++i)
			ProcessMultiplier(vfCompositeCunningham+i*nCandidatesWords,nPrime,vCunninghamMultipliers+i*ySIZE*xSIZE+nPrimeSeq*xSIZE);
    }

    // Fast version
	const unsigned int nLongs = nSieveSize/8 / sizeof(unsigned int);
	for (unsigned int i = 0; i < nLongs; i++)
	{
		vfCandidates[i] = ~((vfCompositeCunningham[0*nCandidatesWords+i] | vfCompositeCunningham[2*nCandidatesWords+i]) &
						(vfCompositeCunningham[1*nCandidatesWords+i] | vfCompositeCunningham[3*nCandidatesWords+i]) &
						(vfCompositeCunningham[0*nCandidatesWords+i] | vfCompositeCunningham[1*nCandidatesWords+i]));
	}
	CandidateList.clear();
	for(uint i=0; i<nLongs; ++i)
	{
		for(uint j=0; j<sizeof(unsigned int)*8; ++j)
		{
			if (vfCandidates[i]&(1<<j))
				CandidateList.push_back(i*sizeof(unsigned int)*8+j);
		}
	}
	
	uint stoptime=ticker();
	cout << stoptime-starttime << "ms for CPU sieve! " << CandidateList.size() << " candidates found!" << endl;
}
