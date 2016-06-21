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

/*
 * rk_lce.hpp
 *
 *  Created on: Jun 4, 2016
 *      Author: Nicola Prezza
 *
 *  Encodes a text with suffix-inclusive Rabin-Karp hash function. This structure supports access to the text and
 *  fast computation of the RK function of any text suffix
 *
 *  Space: n*ceil(log_2 sigma) bits
 *  Supports:
 *  	- access to the text in in O(1) time per character
 *  	- LCP between any two text suffixes in O(log n) time
 *
 *  The class moreover can return a function to lexicographically-compare in O(log n) time any two
 *  text suffixes (useful for suffix-sorting in-place any subset of text positions)
 *
 *  Note: though there is a probability of getting a wrong LCP result due to hash collisions,
 *  this probability is less than 2^-120 for any reasonable-sized text (n<2^128): there is a
 *  (by far) higher chance of getting a wrong result by a cosmic ray flipping a bit in RAM during computation!
 *
 */

#ifndef INTERNAL_RK_LCE_HPP_
#define INTERNAL_RK_LCE_HPP_

#include <includes.hpp>

namespace rk_lce{

template<uint16_t w>
class rk_lce{

public:

	static constexpr uint128 q = (uint128(1)<<w)-1;

	/*
	 * Build RK-LCP structure over the text stored at this path
	 */
	rk_lce(string filename){

	    char_to_uint = vector<uint8_t>(256);
	    uint_to_char = vector<char>(256);

	    vector<bool> mapped(256,false);

	    // DETECT ALPHABET

	    //alphabet size. To be rounded to the next power of 2

	    n = 0;
	    int sigma_temp = 0;
	    {

			ifstream ifs(filename,std::ios::binary);

			uint8_t c;

			while(ifs >> c){

				if(not mapped[c]){

					mapped[c] = true;

					char_to_uint[c] = sigma_temp;
					uint_to_char[sigma_temp] = char(c);

					sigma_temp++;

				}

				n++;

			}

	    }

	    assert(sigma_temp>0);

	    sigma = 1;
	    log2_sigma = 0;

	    while(sigma < sigma_temp){

	    	sigma *= 2;
	    	log2_sigma++;

	    }

	    //if sigma_temp == 0, log2_sigma must be at least 1 bit
	    log2_sigma = log2_sigma == 0 ? 1 : log2_sigma;

	    //find largest B such that sigma^B <= q
	    B = 1;
	    uint128 temp = sigma;

	    while(q/temp >= sigma){

	    	temp *= sigma;
	    	B++;

	    }

	    assert(temp<q);

	    //we add a left padding
	    pad = n%B==0 ? 0 : B - (n%B);
	    N = n+pad;

	    assert(N%B == 0);

	}

	/*
	 * access i-th character of the text
	 *
	 * complexity: O(1)
	 *
	 */
	char operator[](uint64_t i){

		return 0;

	}

	/*
	 * LCE between i-th and j-th suffixes
	 *
	 * complexity:
	 *
	 * - O(1) if the LCE is shorter than T.block_size()
	 *
	 * - O(log n) otherwise
	 *
	 */
	uint64_t LCE(uint64_t i, uint64_t j){

		return 0;

	}

	/*
	 * returns a std::function F to lexicographically compare suffixes: F(i,j) = true iif i-th suffix is < than j-th suffix.
	 * Suffix here are enumerated from left (full text is suffix 0)
	 *
	 * Time: O(log n)
	 *
	 */
	std::function<bool (uint64_t, uint64_t)> lex_less_than() {

	    return [&](uint64_t i, uint64_t j) {

	    	if(i==j) return false;

	    	//rightmost suffix is the shortest
	    	uint64_t min = i<j ? j : i;

	    	auto lce = LCP(i,j);

	    	//in this case shortest suffix is the smallest
	    	if(lce == n-min) return i==min;

	    	assert(i+lce<n);
	    	assert(j+lce<n);

	    	//get characters following LCE
	    	auto ic = operator[](i+lce);
	    	auto jc = operator[](j+lce);

	    	return ic < jc;

	    };

	}

	uint64_t number_of_blocks(){ return 0; }

	uint64_t block_size(){ return B; }
	uint64_t length(){ return n; }

	uint64_t padding(){return pad;}
	uint64_t padded_size(){ return N; }
	uint64_t size(){ return n; }

	uint16_t alphabet_size(){return sigma;}
	uint16_t log_alphabet_size(){return log2_sigma;}

private:




	vector<uint8_t> char_to_uint;
	vector<char> uint_to_char;

	//text length
	uint64_t n;

	//padding at the left of the text to reach a size multiple of B
	//size of the stored text is n+pad
	uint64_t pad;

	//length of padded text = n+pad
	uint64_t N;

	//block size
	uint64_t B;

	//power of 2 immediately greater than or equal to alphabet size
	uint16_t sigma;

	//log_2(sigma)
	uint16_t log2_sigma;

};

}

#endif /* INTERNAL_RK_LCE_HPP_ */
