#ifndef SIEVE_H
#define SIEVE_H

#include "Primes.h"
#include <cmath>

class Sieve
{
public:
	static vector<vector<uint> > sieves; //sieves[prime_ordinal][number_to_check];
	
	//TODO: move this to a misc math util file
	static uint ModDivBy2(uint n, uint m) //modulus must be a prime!
	{
		return (n+m*(n&1))>>1;
	}
	
	static uint Get(uint prime_index, uint number)
	{
		//TODO, could check that number<prime_index. and prime_index<sieves.size().
		return sieves[prime_index][number];
	}
	
	enum KINDS
	{
		FIRST_KIND = 0x01,
		SECOND_KIND = 0x02,
		TWIN_KIND = 0x04
	};
	
	static vector<uint> GenerateSieve(const uint prime, uint chain_length_min)
	{
		vector<uint> sieve = vector<uint>(prime,0);
		vector<uint> first_kind(prime,2000000000);

		uint number=0,current_value=0;
		do
		{
			first_kind[number] = current_value;
			++current_value;
			number = ModDivBy2((number+prime-1)%prime,prime);
		}while(number!=0);

		for(uint i=0; i<prime; ++i)
		{
			uint first = first_kind[(prime+i-1)%prime];
			uint second = first_kind[(prime-i-1)%prime];
			if (first >= chain_length_min)
				sieve[i] |= FIRST_KIND;
			if (second >= chain_length_min)
				sieve[i] |= SECOND_KIND;
			if (min(first,second)*2+1 >= chain_length_min)
				sieve[i] |= FIRST_KIND|SECOND_KIND;
		}
		return sieve;
	}

	static void Generate(uint amount_of_sieves, uint chain_length_min)
	{
		if (Primes::v.size() == 0)
			cout << "Primes must be generated before generating the sieve!" << endl;
		
		for(uint i=0; i<amount_of_sieves; ++i)
		{
			sieves.push_back(GenerateSieve(Primes::v[i],chain_length_min));
		}
	}
};

#endif
