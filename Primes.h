#ifndef PRIMES_H
#define PRIMES_H

class Primes
{
public:
	static vector<uint> v;
	static void Generate(uint max_value)
	{
		cout << "Generating primes up to " << max_value << "." << endl;
		v.clear();
		vector<bool> eratosthenes(max_value,false);
		for(uint i=2; i<max_value; ++i)
			for(uint j=i*i; j<max_value; j+=i)
				eratosthenes[j] = true;
		for(uint i=2; i<max_value; ++i)
			if (!eratosthenes[i])
				v.push_back(i);
	}
};

#endif
