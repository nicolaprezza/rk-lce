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
 * rk_lce_bin.hpp
 *
 *  Created on: Jun 4, 2016
 *      Author: Nicola Prezza
 *
 *  Encodes a binary text with suffix-inclusive Rabin-Karp hash function. This structure supports access to the text and
 *  fast computation of the RK function of any text suffix
 *
 *  Space: n bits (n = size of original text)
 *  Supports:
 *  	- access to the text in in O(1) time per character
 *  	- LCE between any two text suffixes in O(log n) time
 *
 */

#ifndef INTERNAL_RK_LCE_BIN_HPP_
#define INTERNAL_RK_LCE_BIN_HPP_

#include <includes.hpp>
#include <bitv.hpp>

namespace rk_lce{

class rk_lce_bin{

public:

	//For efficiency reasons, the Mersenne number is fixed: 2^127-1
	static constexpr uint16_t w = 127;

	//modulo
	static constexpr uint128 q = (uint128(1)<<w)-1;

	/*
	 * Build RK-LCE structure over the bitvector. Note that the size of
	 * the bitvector must be a multiple of w!
	 */
	rk_lce_bin(vector<bool> & input_bitvector){

		n = input_bitvector.size();

		assert(n%w==0);

		//number of 127-bits blocks
		auto n_bl = n/w;

		assert(n_bl>0);

		//array P_vec: prefix sums of blocks different than q
		vector<uint128> P_vec;

		{

	    	//pack bits in blocks of w bits: array B

			auto B_vec = vector<uint128>(n_bl,0);

			uint64_t i = 0;
			for(auto b:input_bitvector){

				B_vec[i/w] |= (uint128(b) << (w-(i%w+1)) );
				i++;

			}

			//first block must be different than q
			assert(B_vec[0]!=q);

		    //detect full blocks

		    uint64_t number_full_blocks = 0;
		    for(auto bl : B_vec) number_full_blocks += bl==q;

		    //Build bitvector Q1

		    {

				vector<bool> Q_vec;

				i = 0;
				for(auto bl : B_vec) Q_vec.push_back(bl==q);

				Q1 = bitv(Q_vec);

				assert(Q1.rank(Q1.size()) == number_full_blocks);

		    }

			cout << Q1.rank(Q1.size()) << " full blocks"<<endl;
			cout << Q1.rank(Q1.size(),0) << " non-full blocks"<<endl;

		    // NOW COMPUTE PREFIX SUMS (array P' in the paper. here we overwrite B_vec with P')

		    for(uint64_t i=1;i<B_vec.size();++i){

		    	B_vec[i] = (B_vec[i] + mul_pow2<w>(B_vec[i-1],w))%q;

		    }

		    //P_vec is never empty because we require B[0] != q
		    P_vec = vector<uint128>(Q1.rank(Q1.size(),0));

		    uint64_t i_P_vec = 0;
		    for(uint64_t j = 0; j< Q1.size();++j){

		    	if(not Q1[j]){

			    	assert(i_P_vec < P_vec.size());
		    		P_vec[i_P_vec++] = B_vec[j];

		    	}

		    }

		}//destroy B_vec

		P = packed_vector_127(P_vec);

	}

	uint64_t bit_size(){

		return P.bit_size() + Q1.bit_size() + sizeof(this)*8;

	}


	/*
	 * access i-th bit
	 *
	 * complexity: O(1)
	 *
	 */
	bool operator[](uint64_t i){

		assert(i<n);

		//block containing position i
		auto ib = i/w;

		return (B(ib) >> (w-(i%w+1))) & uint128(1);

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

		return LCE_naive(i,j);

		/*
		assert(i<n);
		assert(j<n);

		//same suffix
		if(i==j) return n-i;

		//one of the two suffixes are empty
		if(i == n or j == n) return 0;

		//length of the suffixes
		auto i_len = n-i;
		auto j_len = n-j;

		auto i_block = get_block_from_position(i+pad);
		auto j_block = get_block_from_position(j+pad);

		if(i_block != j_block or i_len <= B or j_len <= B){

			auto leading_zeros = clz_u128(i_block ^ j_block);
			uint64_t lce = leading_zeros / log2_sigma;

			auto min = std::min(i_len,j_len);

			lce = lce > min ? min : lce;

			assert(lce == LCE_naive(i,j));

			return lce;

		}

		auto lce = LCE_binary(i,j);

		if(lce != LCE_naive(i,j)){

			cout << endl << i << ", " << j <<endl;
			cout << "LCE = " << lce << endl;
			cout << "LCE trivial = " << LCE_naive(i,j) << endl;

			for(uint64_t k = i; k<i+100 and k<n;++k) cout << operator[](k); cout << endl;
			for(uint64_t k = j; k<j+100 and k<n;++k) cout << operator[](k); cout << endl;


		}

		assert(lce == LCE_naive(i,j));

		return lce;*/

	}

	/*
	 * true iif (with low probability of failure) T[i, ..., i+l-1] == T[j, ..., j+l-1]
	 *
	 * Time: O(1)
	 *
	 * the fingerprints i_fp and j_fp of T[i,...,n-1] and T[j,...,n-1]
	 * can be provided for more efficiency
	 *
	 */
	/*bool equals(uint64_t i, uint64_t j, uint64_t l, uint128 i_fp = q, uint128 j_fp = q){

		assert(i+l-1<n);
		assert(j+l-1<n);

		i += pad;
		j += pad;

		if(i==j) return true;

		//reorder i and j such that i<j

		if(i>j){

			auto t = i;
			i = j;
			j = t;

		}

		uint128 Xi = get_pref_suf_fp(i,i+l-1, i_fp);
		uint128 Xj = get_pref_suf_fp(j,j+l-1, j_fp);

		uint64_t exp = (j-i)*log2_sigma;

		Xj = mul_pow2<w>(Xj,exp);

		return Xj == Xi;

	}*/


	/*
	 * O(n)-time implementation of LCE
	 */
	uint64_t LCE_naive(uint64_t i, uint64_t j){

		if(i==j) return n-i;

		uint64_t lce = 0;

		while( i+lce < n and j+lce<n and operator[](i+lce) == operator[](j+lce)) lce++;

		return lce;

	}

	uint64_t number_of_blocks(){ return Q1.size(); }

	uint64_t block_size(){ return w; }

	uint64_t length(){ return n; }
	uint64_t size(){ return n; }

private:

	/*
	 * get fingerprint of prefix staring at i-th block (included)
	 * this is array P' of the paper
	 *
	 * complexity: O(log m)
	 *
	 */
	uint128 P1(uint64_t i){

		//if there are no full blocks, speed up computation of P'[i]

		if(Q1.rank(Q1.size()) == 0) return P[i];

		auto t = Q1.predecessor_0(i);

		assert(Q1.rank(t,0) < P.size());
		assert(t<=i);

		return mul_pow2<w>(P[Q1.rank(t,0)],w*(i-t));

	}

	/*
	 * get i-th block
	 * this is array B of the paper
	 *
	 * complexity: O(log m)
	 *
	 */
	uint128 B(uint64_t i){

		assert(i<Q1.size());

		//first block must not be equal to q
		assert(not Q1[0]);

		//if the above is satisfied, P[0] contains the first block
		if(i==0) return P[0];

		//deal separately with blocks equal to q
		if(Q1[i]) return q;

		//if the block is not equal to q, we can compute it with
		//modular operations using array P'
		return sub<w>( P1(i), mul_pow2<w>(P1(i-1),w) );

	}

	/*
	 * LCE between i-th and j-th suffixes
	 *
	 * complexity: O(log n) (a binary search)
	 *
	 */
	/*uint64_t LCE_binary(uint64_t i, uint64_t j){

		//length of the suffixes
		auto i_len = n-i;
		auto j_len = n-j;

		assert(i_len > w);
		assert(j_len > w);

    	auto sc = suffix_comparator(this, i, j);

    	cout << endl;
    	for(int i=0;i< 100;++i) cout << sc[i];cout << endl;

    	//lce = (position of first 1 in bitvector sc) - 1
    	auto lce = uint64_t(std::upper_bound(sc.begin(), sc.end(), false)) - 1;

    	return lce;

	}


	inline int clz_u128 (uint128 u) {
	  uint64_t hi = u>>64;
	  uint64_t lo = u;
	  int retval[3]={
	    __builtin_clzll(hi),
	    __builtin_clzll(lo)+64,
	    128
	  };
	  int idx = !hi + ((!lo)&(!hi));
	  return retval[idx];
	}*/


	/*
	 * returns the B characters following position i (included).
	 * The value is left-aligned: first character begins at the leftmost bit.
	 * A padding of 0 is
	 * added to the right if there are less than B characters on the right of
	 * position i
	 *
	 * Time: O(1)
	 *
	 */
	/*uint128 get_block_from_position(uint64_t i){

		assert(i<N);

		uint16_t l_shift = 128 - B*log2_sigma;
		uint64_t ib = i/B;

		if(i%B==0) return get_block(ib)<<l_shift;

		//resp. number of charachters to take from
		//right and left block
		uint16_t r_c = i%B;
		uint16_t l_c = B - r_c;

		uint64_t nb = blocks.size();

		auto l_block = get_block(ib);//block containing position i
		auto r_block = ib==nb-1 ? 0 : get_block(ib+1);//next block, if it exists

		uint128 X = l_block << log2_sigma*r_c;
		X |= r_block >> log2_sigma*l_c;

		return X<<l_shift;

	}*/

	/*
	 * returns the Rabin-Karp fingerprint of T[i,...,j]000...000 (inclusive),
	 * where there are n - (j+1) zeros on the right
	 *
	 * Time: O(1)
	 *
	 * the fingerprint i_fp of T[i,...,n-1] can be provided for more efficiency
	 *
	 */
	/*uint128 get_pref_suf_fp(uint64_t i, uint64_t j, uint128 i_fp = q){

		return sub<w>(  i_fp == q ? get_suf_fp(i) : i_fp , get_suf_fp(j+1) );

	}*/

	/* get fingerprint of i-th suffix (from the leftmost: full text is suffix 0.
	 * n-th suffix corresponds to empty string and has fingerprint = 0
	 *
	 * O(1) time
	 *
	 */
	/*uint128 get_suf_fp(uint64_t i){

		assert(i<=N);

		if(i==N) return 0;

		uint64_t ib = i/B;
		uint64_t nb = blocks.size();

		//content of the block following the one containing position i (0 if it does not exist)
		uint128 next_block = ib == nb-1 ? 0 : blocks[ib+1];
		uint128 this_block = blocks[ib];

		uint128 R = i%B == 0 ? this_block : next_block;

		uint128 MASK = i%B == 0 ? 0 : ( uint128(1) << (log2_sigma*(B-i%B)) ) - 1;

		//number of blocks on the right of this block
		uint64_t rb = nb - (ib+1);

		uint128 L = i%B == 0 ? 0 : mul_pow2<w>(get_block(ib) & MASK,rb*B*log2_sigma);

		return (L+R)%q;

	}*/

	/*uint128 get_suf_fp1(uint64_t i){

		assert(i<=N);

		if(i==N) return 0;

		i -= pad;

		uint64_t e = 0;

		uint128 X = 0;

		for(int64_t j=n-1;j>=int64_t(i);--j){


			X = (X + mul_pow2<w>(char_to_uint[uint16_t(operator[](j))],e))%q;
			e += log2_sigma;

		}

		return X;

	}*/

	/*
	 * this class is built upon two suffixes i != j. Let i<j (other way is symmetric). The class offers an abstraction
	 * over a virtual vector of bool of length (n - j) + 1. The vector is of the form 000...0011...111. The 0s are in
	 * correspondence of equal prefixes of the two suffixes, while a 1 means that the two prefixes are different. A virtual
	 * terminator symbol (smaller than all characters) is appended at the end of the text to guarantee that there is at least
	 * one '1' at the end of this vector (this is why the vector has size (n - j) + 1 instead of n - j)
	 *
	 * The class implements an iterator so that std binary search can be used on it.
	 *
	 */
	/*class suffix_comparator{

	   class sc_iterator  : public std::iterator<random_access_iterator_tag,bool >{

		   	friend class suffix_comparator;

			suffix_comparator *_sci = nullptr;
			uint64_t _index = 0;

			sc_iterator(suffix_comparator *v, uint64_t index)
				: _sci(v), _index(index) { }

		public:

			sc_iterator() = default;
			sc_iterator(sc_iterator const&) = default;


			// iterator to index number
			operator uint64_t(){
				return _index;
			}

			sc_iterator &operator=(sc_iterator const&) = default;

			// Iterator
			bool operator*() const {
				return (*_sci)[_index];
			}

			sc_iterator &operator++() {
				++_index;
				return *this;
			}

			// EqualityComparable
			bool operator==(sc_iterator it) const {
				return _index == it._index;
			}

			// ForwardIterator
			bool operator!=(sc_iterator it) const {
				return _index != it._index;
			}

			sc_iterator operator++(int) {
				sc_iterator it(*this);
				++_index;
				return it;
			}

			// BidirectionalIterator
			sc_iterator &operator--() {
				--_index;
				return *this;
			}

			sc_iterator operator--(int) {
				sc_iterator it(*this);
				--_index;
				return it;
			}

			// RandomAccessIterator
			sc_iterator &operator+=(uint64_t n) {
				_index += n;
				return *this;
			}

			sc_iterator operator+(uint64_t n) const {
				sc_iterator it(*this);
				it += n;
				return it;
			}

			friend sc_iterator operator+(uint64_t n, sc_iterator it)
			{
				return it + n;
			}

			sc_iterator &operator-=(uint64_t n) {
				_index -= n;
				return *this;
			}

			sc_iterator operator-(uint64_t n) const {
				sc_iterator it(*this);
				it -= n;
				return it;
			}

			friend sc_iterator operator-(uint64_t n, sc_iterator it) {
				return it - n;
			}

			uint64_t operator-(sc_iterator it) {
				return uint64_t(_index) - uint64_t(it._index);
			}

			bool  operator[](uint64_t i) const {
				return (*_sci)[_index + i];
			}

			bool operator<(sc_iterator it) const {
				return _index < it._index;
			}

			bool operator<=(sc_iterator it) const {
				return _index <= it._index;
			}

			bool operator>(sc_iterator it) const {
				return _index > it._index;
			}

			bool operator>=(sc_iterator it) const {
				return _index >= it._index;
			}

		};

	public:

		suffix_comparator(rk_lce_bin<w> * T, uint64_t i, uint64_t j){

			//check these conditions outside this class
			assert(i!=j);
			assert(i<T->size());
			assert(j<T->size());

			if(j<i){

				auto t = i;
				i = j;
				j = t;

			}

			// now i<j holds

			this->T = T;
			this->i = i;
			this->j = j;

			i_fp = T->get_suf_fp(i+T->padding());
			j_fp = T->get_suf_fp(j+T->padding());

			/*
			 * T->size() - j is the length of suffix j. We add a virtual
			 * terminator, so the suffix has length (T->size() - j) + 1.
			 * => there are (T->size() - j) + 2 total possible lengths
			 * (because we include length 0)
			 */


	/*
	 *
			n = (T->size() - j) + 2;
			log2_sigma = T->log_alphabet_size();
			B = T->block_size();

		}

		uint64_t size(){return n;}

		bool operator[](uint64_t t){

			assert(t<n);

			//prefixes of length n-1 are different because of the terminator
			if(t==n-1) return true;
			//empty prefixes are equal
			if(t == 0) return false;

			return not T->equals(i,j,t);
			//return not T->equals(i,j,t,i_fp,j_fp);

		}

		sc_iterator begin() { return sc_iterator(this, 0); }
		sc_iterator end()   { return sc_iterator(this, n); }



	private:

		rk_lce_bin<w> * T;

		uint64_t n;

		uint64_t i;
		uint64_t j;

		uint16_t log2_sigma;
		uint16_t B;

		//fingerprints of i-th and j-th suffixes
		uint128 i_fp;
		uint128 j_fp;

	};*/

	//127-bits blocks. Array P in the paper
	packed_vector_127  P;

	//sparse bitvector marking with a 1 blocks of value q
	//array Q' in the paper
	bitv Q1;

	//bitvector length
	uint64_t n = 0;

};

}

#endif /* INTERNAL_RK_LCE_BIN_HPP_ */
