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
	 * the bitvector must be a multiple of w and first w bits must not be equal to 1111...111
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

		if(i==j) return n-i;

		auto lce = LCE_binary(i,j);

		//assert(lce == LCE_naive(i,j));

		return lce;

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
	 * the fingerprints i_fp and j_fp of T[0,...,i-1] and T[0,...,j-1]
	 * can be provided for more efficiency
	 *
	 */
	bool equals(uint64_t i, uint64_t j, uint64_t l, uint128 i_fp = q, uint128 j_fp = q){

		assert(i+l-1<n);
		assert(j+l-1<n);

		uint128 fpi = RK(i,i+l-1,i_fp);
		uint128 fpj = RK(j,j+l-1,j_fp);

		return fpi == fpj;

	}


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
	 * rabin-karp fingerprint of T[0,...,i]
	 *
	 */
	uint128 RK(uint64_t i){

		auto j = i/w;

		uint128 L = j==0 ? 0 : mul_pow2<w>(P1(j-1),i-j*w+1);
		uint128 R = B(j) >> (w-(i%w+1));

		return (L+R)%q;

	}

	/*
	 * rabin-karp fingerprint of T[i,...,j]
	 *
	 * for efficiency, rki=RK(i-1) can be specified as input
	 *
	 */
	uint128 RK(uint64_t i, uint64_t j, uint128 rki=q){

		assert(j>=i);

		rki = rki != q ? rki : ( i==0 ? 0 : RK(i-1) );

		rki = mul_pow2<w>(rki, j-i+1);

		uint128 rkj = RK(j);

		return sub<w>(rkj,rki);

	}


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
	uint64_t LCE_binary(uint64_t i, uint64_t j){

		assert(i!=j);

		//length of the suffixes
		auto i_len = n-i;
		auto j_len = n-j;

    	auto sc = suffix_comparator(this, i, j);

    	//cout << endl;
    	//for(int i=0;i< 100;++i) cout << sc[i];cout << endl;

    	//lce = (position of first 1 in bitvector sc) - 1
    	auto lce = uint64_t(std::upper_bound(sc.begin(), sc.end(), false)) - 1;

    	return lce;

	}


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
	class suffix_comparator{

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

		suffix_comparator(rk_lce_bin * T, uint64_t i, uint64_t j){

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

			i_fp = i==0 ? 0 : T->RK(i-1);
			j_fp = j==0 ? 0 : T->RK(j-1);

			/*
			 * T->size() - j is the length of suffix j. We add a virtual
			 * terminator, so the suffix has length (T->size() - j) + 1.
			 * => there are (T->size() - j) + 2 total possible lengths
			 * (because we include length 0)
			 */

			n = (T->size() - j) + 2;

		}

		uint64_t size(){return n;}

		bool operator[](uint64_t t){

			assert(t<n);

			//prefixes of length n-1 are different because of the terminator
			if(t==n-1) return true;
			//empty prefixes are equal
			if(t == 0) return false;

			//return not T->equals(i,j,t);
			return not T->equals(i,j,t,i_fp,j_fp);

		}

		sc_iterator begin() { return sc_iterator(this, 0); }
		sc_iterator end()   { return sc_iterator(this, n); }



	private:

		rk_lce_bin * T;

		uint64_t n;

		uint64_t i;
		uint64_t j;

		//fingerprints of (i-1)-th and (j-1)-th prefixes
		uint128 i_fp;
		uint128 j_fp;

	};

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
