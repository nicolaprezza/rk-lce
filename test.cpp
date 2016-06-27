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

	vector<bool> BB {0,1,1,1,1,1,1,1,1,0,0,1,0,0,1,1,1,1,1,0};
	bitv bv(BB);

	for(int i=0;i<bv.size();++i) cout << bv[i];cout << endl;
	for(int i=0;i<bv.size();++i) cout << bv.rank(i) << " ";cout << endl;

	for(int i=0;i<bv.size();++i) cout << bv.predecessor_0(i) << " ";cout << endl;









	exit(0);




	srand(time(NULL));
	uint64_t n = 100000*127;

	vector<bool> B(n);

	//probability of a 1
	double p = 0.95;

	for(int i=0;i<n;++i){

		B[i] = double(rand())/RAND_MAX < p;

	}

	auto lce = rk_lce_bin(B);

	for(int i=0;i<B.size();++i){

		assert(B[i]==lce[i]);

	}

}
