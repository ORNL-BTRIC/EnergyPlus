// Tencent is pleased to support the open source community by making RapidJSON available.
//
// Copyright (C) 2015 THL A29 Limited, a Tencent company, and Milo Yip. All rights reserved.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
// associated documentation files (the "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// The following code is inspired by
// https://github.com/miloyip/rapidjson/blob/1a4cce275d7a855264f714789ed02cdc991dd62a/include/rapidjson/reader.h#L295
// however it has been generalized for use in EnergyPlus.

#ifndef count_until_hpp_INCLUDED
#define count_until_hpp_INCLUDED

#if defined(_MSC_VER)
#include <intrin.h>
#pragma intrinsic(_BitScanForward)
#endif

#include <emmintrin.h>

inline size_t simd_count_until( char const * p, char const c1, char const c2, char const c3, char const c4 ) {
	const __m128i w0 = _mm_set1_epi8( c1 );
	const __m128i w1 = _mm_set1_epi8( c2 );
	const __m128i w2 = _mm_set1_epi8( c3 );
	const __m128i w3 = _mm_set1_epi8( c4 );

	size_t index = 0;

	for (;; index += 16) {
		const __m128i s = _mm_loadu_si128( reinterpret_cast<const __m128i *>( p + index ) );
		__m128i x = _mm_cmpeq_epi8(s, w0);
		x = _mm_or_si128(x, _mm_cmpeq_epi8(s, w1));
		x = _mm_or_si128(x, _mm_cmpeq_epi8(s, w2));
		x = _mm_or_si128(x, _mm_cmpeq_epi8(s, w3));
		unsigned short r = static_cast<unsigned short>(~_mm_movemask_epi8(x));
		if (r != 0) {
#ifdef _MSC_VER
			unsigned long offset;
			_BitScanForward(&offset, r);
			index += offset;
			return index;
			// return p + offset;
#else
			index += __builtin_ffs(r) - 1;
			return index;
#endif
		}
	}
	return index;
}

inline size_t simd_count_until( char const * p, char const c1, char const c2 ) {
	const __m128i w0 = _mm_set1_epi8( c1 );
	const __m128i w1 = _mm_set1_epi8( c2 );

	size_t index = 0;

	for (;; index += 16) {
		const __m128i s = _mm_loadu_si128( reinterpret_cast<const __m128i *>( p + index ) );
		__m128i x = _mm_cmpeq_epi8(s, w0);
		x = _mm_or_si128(x, _mm_cmpeq_epi8(s, w1));
		unsigned short r = static_cast<unsigned short>(~_mm_movemask_epi8(x));
		if (r != 0) {
#ifdef _MSC_VER
			unsigned long offset;
			_BitScanForward(&offset, r);
			index += offset;
			return index;
			// return p + offset;
#else
			index += __builtin_ffs(r) - 1;
			return index;
#endif
		}
	}
	return index;
}

inline size_t simd_count_until( char const * p, char const c1 ) {
	const __m128i w0 = _mm_set1_epi8( c1 );

	size_t index = 0;

	for (;; index += 16) {
		const __m128i s = _mm_loadu_si128( reinterpret_cast<const __m128i *>( p + index ) );
		__m128i x = _mm_cmpeq_epi8(s, w0);
		unsigned short r = static_cast<unsigned short>(~_mm_movemask_epi8(x));
		if (r != 0) {
#ifdef _MSC_VER
			unsigned long offset;
			_BitScanForward(&offset, r);
			index += offset;
			return index;
			// return p + offset;
#else
			index += __builtin_ffs(r) - 1;
			return index;
#endif
		}
	}
	return index;
}

#endif // count_until_hpp_INCLUDED
