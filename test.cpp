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

using namespace std;
using namespace rk_lce;

string output_file = string();
string input_text = string();
string input_pos = string();

int main(int argc, char** argv){

	srand(time(NULL));
	uint64_t n = (rand()%1000000)*127;

	vector<bool> B(n);

	for(int i=0;i<n;++i){

		B[i] = rand()%2;

	}

	auto lce = rk_lce_bin(B);

	for(int i=0;i<B.size();++i){

		assert(B[i]==lce[i]);

	}

}
