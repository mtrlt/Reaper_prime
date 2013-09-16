//max size for the FixedFactor is 32*10 bits
//max chain length is 16
//max sieve size 2^24 = ~16.7M

#define NHALFCHAINLENGTH ((NCHAINLENGTH+1)>>1)

#define ORIGIN_LENGTH 11

#if (ORIGIN_LENGTH != 11)
#error ORIGIN_LENGTH must be 11, other sizes are not supported yet.
#endif

uint GetWordNum(uint nBitNum) {
	return nBitNum >> 5; //32 bits in a uint
}

uint GetBitMask(uint nBitNum) {
	return 1UL << (nBitNum & 31);
}

//vMultipliers is thread-local
void AddMultiplier(__local uint* vMultipliers, const uint nSolvedMultiplier)
{
	uint broken=0;
    // Eliminate duplicates
#pragma unroll
    for (uint i = 0; i < NHALFCHAINLENGTH; i++)
    {
		uint temp = vMultipliers[WORKSIZE*i];
		bool result = !broken && (temp == 0xFFFFFFFF || temp == nSolvedMultiplier);
		vMultipliers[WORKSIZE*i] = select(temp,nSolvedMultiplier,(uint)result);
		broken = select(broken,1U,(uint)result);
    }
}

uint mod(uint a, uint b)
{
	return a%b;
}
uint mulmod(uint a, uint m)
{
	return (((ulong)a)<<32)%m;
}

uint modulo(uint16 n, uint m)
{
	uint result=0;
	result=n.sa%m;

	result=mulmod(result,m);
	result += mod(n.s9,m);
	result = select(result,result-m,result>=m);
	
	result=mulmod(result,m);
	result += mod(n.s8,m);
	result = select(result,result-m,result>=m);
	
	result=mulmod(result,m);
	result += mod(n.s7,m);
	result = select(result,result-m,result>=m);
	
	result=mulmod(result,m);
	result += mod(n.s6,m);
	result = select(result,result-m,result>=m);

	result=mulmod(result,m);
	result += mod(n.s5,m);
	result = select(result,result-m,result>=m);

	result=mulmod(result,m);
	result += mod(n.s4,m);
	result = select(result,result-m,result>=m);

	result=mulmod(result,m);
	result += mod(n.s3,m);
	result = select(result,result-m,result>=m);

	result=mulmod(result,m);
	result += mod(n.s2,m);
	result = select(result,result-m,result>=m);

	result=mulmod(result,m);
	result += mod(n.s1,m);
	result = select(result,result-m,result>=m);

	result=mulmod(result,m);
	result += mod(n.s0,m);
	result = select(result,result-m,result>=m);

	return result;
}

uint int_invert(uint a, uint nPrime)
{
    int rem0 = nPrime, rem1 = a;
    int aux0 = 0, aux1 = 1;
    int quotient, inverse=0;

	bool go_on = true;
    while (go_on)
    {
		bool test = (rem1<=1);
		inverse = select(inverse,aux1,(int)test);
		go_on = !test;

        quotient = rem0 / rem1;
		rem0 -= quotient * rem1;
        aux0 -= quotient * aux1;

		test = go_on && (rem0 <= 1);
		inverse = select(inverse,aux0,(int)test);
		go_on = select((int)go_on,false,(int)test);
		
        quotient = rem1 / rem0;
		rem1 -= quotient * rem0;
        aux1 -= quotient * aux0;
    }
	return select((uint)inverse,inverse+nPrime,inverse<0);
}

__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void CalculateMultipliers(__global const uint*restrict vPrimes, const uint nPrimes, const uint16 nFixedFactor, const uint nCandidatesWords, __global uint* vfCunninghamMultipliers)
{
	__local uint localmem[4*NHALFCHAINLENGTH*WORKSIZE];
	
#pragma unroll
	for(uint i=0; i<4*NHALFCHAINLENGTH*WORKSIZE; i+=WORKSIZE)
		localmem[i+get_local_id(0)] = 0xFFFFFFFF;
	
	__local uint* localCunninghamMultipliers = localmem+get_local_id(0);

	uint nPrimeSeq = get_global_id(0)+1;
    {
        uint nPrime = vPrimes[nPrimeSeq];
        uint nFixedFactorMod = modulo(nFixedFactor, nPrime);
        if (nFixedFactorMod != 0)
		{
			uint nFixedInverse;
			nFixedInverse = int_invert(nFixedFactorMod, nPrime);

			// Weave the sieve for the prime
#pragma unroll
			for (uint nBiTwinSeq = 0; nBiTwinSeq < 2 * NCHAINLENGTH; nBiTwinSeq+=2)
			{
				uint arrnum;
				uint nSolvedMultiplier;
				
				arrnum = select(1,3,nBiTwinSeq+1 >= NCHAINLENGTH);
				nSolvedMultiplier = select(nPrime-nFixedInverse,0U,nFixedInverse==0);
				AddMultiplier(localCunninghamMultipliers+WORKSIZE*arrnum*NHALFCHAINLENGTH, nSolvedMultiplier);

				arrnum = select(0,2,nBiTwinSeq >= NCHAINLENGTH);
				nSolvedMultiplier = nFixedInverse;
				AddMultiplier(localCunninghamMultipliers+WORKSIZE*arrnum*NHALFCHAINLENGTH, nSolvedMultiplier);
				
				//modular divide by 2
				nFixedInverse = select(nFixedInverse>>1,(nFixedInverse+nPrime)>>1,nFixedInverse&1);
			}
			
			//vCunninghamMultipliers+i*ySIZE*xSIZE+nPrimeSeq*xSIZE
			
#pragma unroll
			for(uint h=0; h<4; ++h)
			{
				uint global_base = h*nPrimes*NHALFCHAINLENGTH+nPrimeSeq*NHALFCHAINLENGTH;
				uint local_base = WORKSIZE*h*NHALFCHAINLENGTH;
#pragma unroll
				for(uint i=0; i<NHALFCHAINLENGTH; ++i)
					vfCunninghamMultipliers[global_base+i] = localCunninghamMultipliers[local_base+WORKSIZE*i];
			}
		}
    }
}

#define SIEVE_PER_GROUP (LOCAL_MEM_USED*8)

__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void Sieve(__global uint* vfCompositeCunningham, __global const uint* vfCunninghamMultipliers, __global const uint* vPrimes, const uint nPrimes, const uint nCandidatesWords)
{
	__local uint sievepart[SIEVE_PER_GROUP/32];
	
#pragma unroll
	for(uint g=0; g<4; ++g)
	{
#pragma unroll
		for(uint i=0; i<SIEVE_PER_GROUP/32; i+=WORKSIZE)
			sievepart[i+get_local_id(0)] = 0;
		
		uint sievepartition_start=get_group_id(0)*SIEVE_PER_GROUP;
		uint sievepartition_end=get_group_id(0)*SIEVE_PER_GROUP+SIEVE_PER_GROUP;
		sievepartition_end = min(sievepartition_end,(uint)NSIEVESIZE);

		for(uint h=get_local_id(0); h<nPrimes-1; h+=WORKSIZE)
		{
			uint global_base = g*nPrimes*NHALFCHAINLENGTH+(h+1)*NHALFCHAINLENGTH;
			uint nPrime = vPrimes[h+1];
			const uint nRotateBits = nPrime & 31;
			uint nVariableMultiplierStart = sievepartition_start-sievepartition_start%nPrime;
			
			uint broken_mask=0xFFFFFFFF;
#pragma unroll
			for (uint i = 0; i < NHALFCHAINLENGTH; i++)
			{
				uint array_value = vfCunninghamMultipliers[global_base+i];
				broken_mask = select(broken_mask,0U,array_value==0xFFFFFFFF);
				uint nVariableMultiplier = nVariableMultiplierStart+array_value;
				nVariableMultiplier = select(nVariableMultiplier,nVariableMultiplier+nPrime,nVariableMultiplier<sievepartition_start);
				uint lBitMask = GetBitMask(nVariableMultiplier);
				for (; broken_mask!=0 && nVariableMultiplier < sievepartition_end;)
				{
					uint inner_mask=0xFFFFFFFFU;
					atomic_or(sievepart+GetWordNum(nVariableMultiplier-sievepartition_start),lBitMask&broken_mask&inner_mask);
					lBitMask = rotate(lBitMask,nRotateBits);
					nVariableMultiplier += nPrime;
					inner_mask = select(0U,0xFFFFFFFFU,nVariableMultiplier<sievepartition_end);
					
					atomic_or(sievepart+GetWordNum(nVariableMultiplier-sievepartition_start),lBitMask&broken_mask&inner_mask);
					lBitMask = rotate(lBitMask,nRotateBits);
					nVariableMultiplier += nPrime;
					inner_mask = select(0U,0xFFFFFFFFU,nVariableMultiplier<sievepartition_end);

					atomic_or(sievepart+GetWordNum(nVariableMultiplier-sievepartition_start),lBitMask&broken_mask&inner_mask);
					lBitMask = rotate(lBitMask,nRotateBits);
					nVariableMultiplier += nPrime;
					inner_mask = select(0U,0xFFFFFFFFU,nVariableMultiplier<sievepartition_end);

					atomic_or(sievepart+GetWordNum(nVariableMultiplier-sievepartition_start),lBitMask&broken_mask&inner_mask);
					lBitMask = rotate(lBitMask,nRotateBits);
					nVariableMultiplier += nPrime;
				}
			}
		}
#pragma unroll
		for(uint i=0; i<(sievepartition_end-sievepartition_start)/32; i+=WORKSIZE)
			vfCompositeCunningham[g*nCandidatesWords+sievepartition_start/32+i+get_local_id(0)] = sievepart[i+get_local_id(0)];
	}
}

__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void Combine(__global const uint* vfCompositeCunningham, const uint nCandidatesWords, __global uint*restrict counter, __global uint*restrict CandidateList)
{
	uint i=get_global_id(0);
	uint n0=vfCompositeCunningham[0*nCandidatesWords+i];
	uint n1=vfCompositeCunningham[1*nCandidatesWords+i];
	uint n2=vfCompositeCunningham[2*nCandidatesWords+i];
	uint n3=vfCompositeCunningham[3*nCandidatesWords+i];
	uint bits = ~((n0|n2) & (n1|n3) & (n0|n1));

	//get rid of the candidate 0. the ALU is way underused in this kernel, so doing this doesn't slow anything down. i love tricks like this :)
	bits &= select(~0,~1,i==0);

#pragma unroll
	for(uint j=0; j<32; ++j)
	{
		if (bits&(1<<j))
		{
			uint counter_value = atomic_inc(counter);
			CandidateList[counter_value] = i*32+j;
		}
	}
}

//Fermat test code:
#define LoopDown(x) x(f); x(e); x(d); x(c); x(b); x(a); x(9); x(8); x(7); x(6); x(5); x(4); x(3); x(2); x(1); x(0); 
#define LoopUp(x) x(0); x(1); x(2); x(3); x(4); x(5); x(6); x(7); x(8); x(9); x(a); x(b); x(c); x(d); x(e); x(f); 

#define CompareIter(n) \
	highertmp |= select(0U,1U<<0x ## n,a.s ## n>b.s ## n);\
	lowertmp |= select(0U,1U<<0x ## n,a.s ## n<b.s ## n);


bool Compare(uint16 a, uint16 b, const bool lower, const bool equal, const bool higher)
{
	uint lowertmp=0,highertmp=0;
	
	LoopDown(CompareIter);
	bool result = select(select((uint)lower,(uint)higher,highertmp>lowertmp),(uint)equal,lowertmp==highertmp);
	
	return result;
}

bool a_LT_b(uint16 a, uint16 b)
{
	return Compare(a,b,true,false,false);
}
bool a_EQ_b(uint16 a, uint16 b)
{
	return Compare(a,b,false,true,false);
}
bool a_GT_b(uint16 a, uint16 b)
{
	return Compare(a,b,false,false,true);
}

uint GetLog2(uint v)
{
	uint r; // result of log2(v) will go here
	uint shift;
	r = select(0,1<<4,v > 0xFFFF); v >>= r;
	shift = select(0,1<<3,v > 0xFF); v >>= shift; r |= shift;
	shift = select(0,1<<2,v > 0xF ); v >>= shift; r |= shift;
	shift = select(0,1<<1,v > 0x3 ); v >>= shift; r |= shift;
	r |= (v >> 1);
	return r;
}

uint GetLog2_16bit(uint v)
{
	uint r; // result of log2(v) will go here
	uint shift;
	shift = select(0,1<<3,v > 0xFF); v >>= shift; r = shift;
	shift = select(0,1<<2,v > 0xF ); v >>= shift; r += shift;
	shift = select(0,1<<1,v > 0x3 ); v >>= shift; r += shift;
	r |= (v >> 1);
	return r;
}


#define GetLengthIter(n) \
	found = select((uint)found,(uint)true,(uint)(a.s ## n != 0));\
	length = select(length-1U,length,(uint)found);

uint GetLength(uint16 a)
{
	uint4 mask=0;

	mask.x |= select(0,1<< 0,a.s0!=0);
	mask.y |= select(0,1<< 1,a.s1!=0);
	mask.z |= select(0,1<< 2,a.s2!=0);
	mask.w |= select(0,1<< 3,a.s3!=0);

	mask.x |= select(0,1<< 4,a.s4!=0);
	mask.y |= select(0,1<< 5,a.s5!=0);
	mask.z |= select(0,1<< 6,a.s6!=0);
	mask.w |= select(0,1<< 7,a.s7!=0);

	mask.x |= select(0,1<< 8,a.s8!=0);
	mask.y |= select(0,1<< 9,a.s9!=0);
	mask.z |= select(0,1<<10,a.sa!=0);
	mask.w |= select(0,1<<11,a.sb!=0);

	mask.x |= select(0,1<<12,a.sc!=0);
	mask.y |= select(0,1<<13,a.sd!=0);
	mask.z |= select(0,1<<14,a.se!=0);
	mask.w |= select(0,1<<15,a.sf!=0);

	return GetLog2_16bit(mask.x|mask.y|mask.z|mask.w)+1;
}

#define GetBitLengthIter(n) \
	bits = select(bits,a.s ## n,(uint)(a.s ## n != 0 && !found));\
	length = select(length-32U,length,(uint)found);\
	found = select((uint)found,(uint)true,(uint)(a.s ## n != 0));

uint GetBitLength(uint16 a)
{
	uint length=512;
	uint bits=0;
	
	bool found=false;
	
	LoopDown(GetBitLengthIter);
	return length+GetLog2(bits)+1;
}

uint16 a_MINUS1_fast(uint16 a)
{
	a.s0 -= 1;
	return a;
}

uint GetIndex(const uint16 a, const uint i)
{
	uint8 b = select(a.s01234567,a.s89abcdef,(uint8)(select(0U,0xFFFFFFFFU,i&8)));
	uint4 c = select(b.s0123    ,b.s4567    ,(uint4)(select(0U,0xFFFFFFFFU,i&4)));
	uint2 d = select(c.s01      ,c.s23      ,(uint2)(select(0U,0xFFFFFFFFU,i&2)));
	uint  e = select(d.s0       ,d.s1       ,       (select(0U,0xFFFFFFFFU,i&1)));
	return e;
}

uint16 SetIndex(const uint a, const uint i)
{
	uint16 ret = (uint16)a;
	
	uint16 mask8 = (uint16)((uint8)0,(uint8)0xFFFFFFFF);
	uint16 mask4 = (uint16)((uint4)0,(uint4)0xFFFFFFFF,(uint4)0,(uint4)0xFFFFFFFF);
	uint16 mask2 = (uint16)((uint2)0,(uint2)0xFFFFFFFF,(uint2)0,(uint2)0xFFFFFFFF,(uint2)0,(uint2)0xFFFFFFFF,(uint2)0,(uint2)0xFFFFFFFF);
	uint16 mask1 = (uint16)(0,0xFFFFFFFF,0,0xFFFFFFFF,0,0xFFFFFFFF,0,0xFFFFFFFF,0,0xFFFFFFFF,0,0xFFFFFFFF,0,0xFFFFFFFF,0,0xFFFFFFFF);
	
	ret = select((uint16)0,ret,mask8^select(0xFFFFFFFFU,0U,i&8));
	ret = select((uint16)0,ret,mask4^select(0xFFFFFFFFU,0U,i&4));
	ret = select((uint16)0,ret,mask2^select(0xFFFFFFFFU,0U,i&2));
	ret = select((uint16)0,ret,mask1^select(0xFFFFFFFFU,0U,i&1));
	
	return ret;
}

uint lessthan(ulong a, ulong b)
{
	bool lt_hi = (a>>32)<(b>>32);
	bool lt_lo = (a&0xFFFFFFFF)<(b&0xFFFFFFFF);
	bool eq_hi = ((a>>32)==(b>>32));
	bool ret = select((uint)lt_hi,(uint)lt_lo,(uint)eq_hi);
	return (uint)ret;
}

#define MulloIter(n) \
	{\
		uint tmp_higher=0;\
		const int minj = max(0,0x ## n+1-blength);\
		const int maxj = min(0x ## n+1,alength);\
		for(int j=minj; j<maxj; ++j)\
		{\
			ulong oldtmp = tmp;\
			tmp += ((ulong)(GetIndex(b,0x ## n-j)))*GetIndex(a,j);\
			tmp_higher += select(0U,1U,lessthan(tmp,oldtmp));\
		}\
		(*lo).s ## n = (uint)tmp;\
		tmp = upsample(tmp_higher,tmp>>32);\
	}
	
#define MulhiIter(n) \
	{\
		uint tmp_higher=0;\
		const int minj = max(0,0x1 ## n+1-blength);\
		const int maxj = min(0x1 ## n+1,alength);\
		for(int j=minj; j<maxj; ++j)\
		{\
			ulong oldtmp = tmp;\
			tmp += ((ulong)(GetIndex(b,0x ## n-j)))*GetIndex(a,j);\
			tmp_higher += select(0U,1U,lessthan(tmp,oldtmp));\
		}\
		(*hi).s ## n =(uint)tmp;\
		tmp = upsample(tmp_higher,tmp>>32);\
	}

void Mul(uint16* hi, uint16* lo, uint16 a, uint16 b)
{
	ulong tmp=0;
	
	const int alength=ORIGIN_LENGTH+1;
	const int blength=ORIGIN_LENGTH+1;
	
	const int combined_length = alength+blength;
	LoopUp(MulloIter);
#if (ORIGIN_LENGTH >= 8)
	{
		uint tmp_higher=0;
		const int minj = 16-ORIGIN_LENGTH;
		const int maxj = ORIGIN_LENGTH+1;
#pragma unroll
		for(int j=minj; j<maxj; ++j)
		{
			ulong oldtmp = tmp;
			tmp += ((ulong)(GetIndex(b,0x0-j)))*GetIndex(a,j);
			tmp_higher += select(0U,1U,lessthan(tmp,oldtmp));
		}
		(*hi).s0 =(uint)tmp;
		tmp = upsample(tmp_higher,tmp>>32);
	}



	//MulhiIter(0);
	MulhiIter(1);
#endif
#if (ORIGIN_LENGTH >= 9)
	MulhiIter(2);
	MulhiIter(3);
#endif
#if (ORIGIN_LENGTH >= 10)
	MulhiIter(4);
	MulhiIter(5);
#endif
#if (ORIGIN_LENGTH >= 11)
	MulhiIter(6);
	MulhiIter(7);
#endif
#if (ORIGIN_LENGTH >= 12)
	MulhiIter(8);
	MulhiIter(9);
#endif
#if (ORIGIN_LENGTH >= 13)
	MulhiIter(a);
	MulhiIter(b);
#endif
#if (ORIGIN_LENGTH >= 14)
	MulhiIter(c);
	MulhiIter(d);
#endif
#if (ORIGIN_LENGTH >= 15)
	MulhiIter(e);
	MulhiIter(f);
#endif
}
void Mul_lower(uint16* lo, uint16 a, uint16 b)
{
	ulong tmp=0;
	
#pragma unroll
	for(int i=0; i<ORIGIN_LENGTH+1; ++i)
	{
		uint tmp_higher=0;
#pragma unroll
		for(int j=0; j<i+1; ++j)
		{
			ulong oldtmp = tmp;
			tmp += ((ulong)(GetIndex(b,i-j)))*GetIndex(a,j);
			tmp_higher += select(0U,1U,lessthan(tmp,oldtmp));
		}
		*lo |= SetIndex((uint)tmp,i);
		tmp = upsample(tmp_higher,tmp>>32);
	}
}

uint16 Mul_512b_32b(const uint16 a, uint b)
{
	ulong overflow=0;
	uint a_size = GetLength(a);
	uint16 ret=0;
	for(uint i=0; i<a_size; ++i)
	{
		overflow += ((ulong)(GetIndex(a,i)))*b;
		ret |= SetIndex((uint)overflow,i);
		overflow >>= 32;
	}
	ret |= SetIndex((uint)overflow,a_size);
	return ret;
}
bool Mul_512b_32b_compare(uint16 a, uint b, const uint16 compared)
{
	a = Mul_512b_32b(a,b);
	bool result = !a_LT_b(a,compared);
	return result;
}

uint16 Sub(uint16 a, uint16 b)
{
	uint a_size = GetLength(a);

	uint overflow=0, overflow_old=0;
	
	uint16 ret=0;
	
	for(uint i=0; i<a_size; ++i)
	{
		uint a_i = GetIndex(a,i);
		uint b_i = GetIndex(b,i);
		overflow_old = overflow;
		overflow = select(0U,1U,a_i<b_i+overflow_old);
		a_i -= b_i+overflow_old;
		ret |= SetIndex(a_i,i);
	}
	return ret;
}

uint16 Calculate_Mu(uint16 m, uint bits)
{
	uint onebit = bits+2; //the number of the bit that should be 1
	uint16 k = 0;
	k.s0 = 1<<(onebit&31);
	uint16 mu = 0;

	for(int i=onebit; i>=0; i-=32)
	{
		if (!a_LT_b(k,m))
		{
			uint x=0;
			//now k >= m
			//find highest x such that m*x<=k
			for(uint j=0; j<32; ++j)
			{
				uint coef = 0x80000000U>>j;
				bool gt_eq = Mul_512b_32b_compare(m,x|coef,k);
				x = select(x|coef,x,(uint)gt_eq);
			}
			if ((i>>5) < 16)
				mu |= SetIndex(x,i>>5);
			k = Sub(k,Mul_512b_32b(m,x));
		}
		k = k.sf0123456789abcde;
	}
	return mu;
}

uint16 SHR_big(uint16 hi, uint16 lo, uint n)
{
	while(n>=32)
	{
		lo = lo.s123456789abcdef0;
		hi = hi.s123456789abcdef0;
		lo.sf = hi.sf;
		n-=32;
	}
	
	if (n>0)
	{
		uint losf_tmp = lo.sf;
		
		lo = (lo>>n) | (lo.s123456789abcdef0 << 32-n);
		lo.sf = (losf_tmp >> n) | (hi.s0 << 32-n);
	}
	return lo;
}

uint GetFirstBitsMask(const uint bits, const uint wordnum)
{
	const uint word_start_bit = wordnum*32;
	const uint word_end_bit = wordnum*32+31;
	
	uint belowmask = 0;
	uint abovemask = 0xFFFFFFFFU;
	uint betweenmask = 0xFFFFFFFFU>>(32-(bits-word_start_bit));
	
	uint mask = betweenmask;
	mask = select(mask,belowmask,bits <= word_start_bit);
	mask = select(mask,abovemask,bits > word_end_bit);
	return mask;
}

#define FirstBitsIter(n) \
	a.s ## n &= GetFirstBitsMask(bits, 0x0 ## n);

uint16 FirstBits(uint16 a, uint bits)
{
	LoopUp(FirstBitsIter);
	return a;
}

uint16 Mod(uint16 hi, uint16 lo, uint16 modulus, uint modulus_length, uint16 mu)
{
	uint k = modulus_length;
	
	uint16 q = SHR_big(hi,lo,k);
	
	uint16 hi2=0,lo2=0;
	Mul(&hi2,&lo2,q,mu);
	q = SHR_big(hi2,lo2,k+2);
	uint16 tmp1 = FirstBits(lo,k+1);
	
	lo2=0;
	Mul_lower(&lo2,q,modulus);
	uint16 tmp2 = FirstBits(lo2,k+1);

	uint16 x = 0;
	
	bool altb = a_LT_b(tmp1,tmp2);
	if (altb)
	{
		tmp2 = Sub(tmp2,tmp1);
		x = SetIndex(1<<((k+1)&31),(k+1)>>5);
	}
	uint16 selected = select(tmp1,x,(uint16)(select(0U,0xFFFFFFFFU,(uint)altb)));
	x = Sub(selected,tmp2);
	x = select(Sub(x,modulus),x,(uint16)select(0U,0xFFFFFFFFU,(uint)a_LT_b(x,modulus)));
	x = select(Sub(x,modulus),x,(uint16)select(0U,0xFFFFFFFFU,(uint)a_LT_b(x,modulus)));
	return x;
}

uint16 MulMod(uint16 a, uint16 b, uint16 m, uint modulus_length, uint16 mu)
{
	uint16 hi=0,lo=0;
	Mul(&hi, &lo, a,b);
	return Mod(hi,lo,m,modulus_length,mu);
}
uint16 SquareMod(uint16 a, uint16 m, uint modulus_length, uint16 mu)
{
	uint16 hi=0,lo=0;
	Mul(&hi, &lo, a,a);
	return Mod(hi,lo,m,modulus_length,mu);
}

uint FermatTest(uint16 modulus)
{
	uint16 result = (uint16)(1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0);
	
	uint16 base = (uint16)(2,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0);
	uint16 exponent = a_MINUS1_fast(modulus);
	
	uint modulus_length = GetBitLength(modulus);
	uint16 mu = Calculate_Mu(modulus,modulus_length*2);

	uint exponent_size = ORIGIN_LENGTH+1;
	
	for(uint i=0; i<exponent_size; ++i)
	{
		uint current_exponent = GetIndex(exponent,i);
		for(uint b=0; b<32; ++b)
		{
			uint16 tmp = MulMod(result,base,modulus,modulus_length,mu);
			result = select(result,tmp,(uint16)(select(0U,0xFFFFFFFFU,(current_exponent>>b)&1)));
			base = SquareMod(base,modulus,modulus_length,mu);
		}
	}
	return result.s0; //good enough for our purposes, results in extremely rare false positives when the Fermat Test remainder is actually of the form X*2^32+1 where X is a positive integer.
}

__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void Fermat(__global const uint*restrict FermatNumbers, __global uint*restrict FermatOutput, const uint16 FixedMultiplier)
{
	uint16 a = Mul_512b_32b(FixedMultiplier,FermatNumbers[get_global_id(0)>>1]);

	//at this point, a has the chain origin. even threads subtract one, odd threads add one. this way primes on both sides of the origin are tested.
	uint added=select(1U,0xFFFFFFFFU,(get_global_id(0)&1) == 0);
	uint compared=select(0U,0xFFFFFFFFU,(get_global_id(0)&1) == 0);
	a.s0 += added;
	a.s1 += select(0U,added,a.s0==compared);
	a.s2 += select(0U,added,a.s1==compared);
	//could add more digits, I suppose. won't really change anything.

	uint ret = FermatTest(a);
	FermatOutput[get_global_id(0)] = ret;
}
