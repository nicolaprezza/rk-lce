/*
 *  This file is part of rk-lce.
 *  Copyright (c) by
 *  Nicola Prezza <nicolapr@gmail.com>
 *
 *   rk-lce is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   rk-lce is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 */

#include <includes.hpp>
#include "internal/rk_lce_bin.hpp"
#include "internal/bitv.hpp"

using namespace std;
using namespace rk_lce;

string output_file = string();
string input_text = string();
string input_pos = string();

int main(int argc, char** argv){

	srand(time(NULL));
	uint64_t n = 10000*127;

	vector<bool> B(n);

	//probability of a 1
	double p = 0.99;

	for(int i=0;i<n;++i){

		if(i<10)
			B[i] = 0;
		else
			B[i] = double(rand())/RAND_MAX < p;

	}

	auto lce = rk_lce_bin(B);

	cout << "Size of the input: " << n/8 << " bytes" << endl;
	cout << "Size of the structure: " << lce.bit_size()/8 << " bytes" << endl;

	for(int i=0;i<B.size();++i){

		assert(B[i]==lce[i]);

	}

}
