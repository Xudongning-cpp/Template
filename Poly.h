#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wdangling-else"
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#include <bits/stdc++.h>
using namespace std;
namespace FastIO {
#define USE_FastIO
#if ( defined(LOCAL) || defined(_WIN32) ) && !defined(DISABLE_MMAP)
#define DISABLE_MMAP
#endif
#ifdef LOCAL
	inline constexpr void _chk_i() {}
	inline char _gc_nochk() {
		return getchar();
	}
	inline char _gc() {
		return getchar();
	}
	inline void _chk_o() {}
	inline void _pc_nochk(char c) {
		putchar(c);
	}
	inline void _pc(char c) {
		putchar(c);
	}
	template < int n > inline void _pnc_nochk(const char* c) {
		for ( int i = 0 ; i < n ; i++ ) putchar(c[i]);
	}
#else
#ifdef DISABLE_MMAP
	inline constexpr int _READ_SIZE = 1 << 18;
	inline static char _read_buffer[_READ_SIZE + 40], *_read_ptr = nullptr, *_read_ptr_end = nullptr;
	static inline bool _eof = false;
	inline void _chk_i() {
		if ( __builtin_expect(!_eof, true) && __builtin_expect(_read_ptr_end - _read_ptr < 40, false) ) {
			int sz = _read_ptr_end - _read_ptr;
			if ( sz ) memcpy(_read_buffer, _read_ptr, sz);
			char* beg = _read_buffer + sz;
			_read_ptr = _read_buffer, _read_ptr_end = beg + fread(beg, 1, _READ_SIZE, stdin);
			if ( __builtin_expect(_read_ptr_end != beg + _READ_SIZE, false) ) _eof = true, *_read_ptr_end = EOF;
		}
	}
	inline char _gc_nochk() {
		return __builtin_expect(_eof && _read_ptr == _read_ptr_end, false) ? EOF : *_read_ptr++;
	}
	inline char _gc() {
		_chk_i();
		return _gc_nochk();
	}
#else
#include<sys/mman.h>
#include<sys/stat.h>
	inline static char* _read_ptr = (char*)mmap(nullptr, [] { struct stat s; return fstat(0, &s), s.st_size; } (), 1, 2, 0, 0);
	inline constexpr void _chk_i() {}
	inline char _gc_nochk() {
		return *_read_ptr++;
	}
	inline char _gc() {
		return *_read_ptr++;
	}
#endif
	inline constexpr int _WRITE_SIZE = 1 << 18;
	inline static char _write_buffer[_WRITE_SIZE + 40], *_write_ptr = _write_buffer;
	inline void _chk_o() {
		if ( __builtin_expect(_write_ptr - _write_buffer > _WRITE_SIZE, false) ) fwrite(_write_buffer, 1, _write_ptr - _write_buffer, stdout), _write_ptr = _write_buffer;
	}
	inline void _pc_nochk(char c) {
		*_write_ptr++ = c;
	}
	inline void _pc(char c) {
		*_write_ptr++ = c, _chk_o();
	}
	template < int n > inline void _pnc_nochk(const char* c) {
		memcpy(_write_ptr, c, n), _write_ptr += n;
	}
	inline struct _auto_flush {
		inline ~_auto_flush() {
			fwrite(_write_buffer, 1, _write_ptr - _write_buffer, stdout);
		}
	} _auto_flush;
#endif
#define println println_
	template < class T > inline constexpr bool _is_signed = numeric_limits < T >::is_signed;
	template < class T > inline constexpr bool _is_unsigned = numeric_limits < T >::is_integer && !_is_signed < T >;
#if __SIZEOF_LONG__ == 64
	template <> inline constexpr bool _is_signed < __int128 > = true;
	template <> inline constexpr bool _is_unsigned < __uint128_t > = true;
#endif
	inline bool _isgraph(char c) {
		return c >= 33;
	}
	inline bool _isdigit(char c) {
		return 48 <= c && c <= 57;
	}
	constexpr struct _table {
#ifndef LOCAL
		int i[65536];
#endif
		char o[40000];
		constexpr _table() :
#ifndef LOCAL
			i {},
#endif
			o {} {
#ifndef LOCAL
			for ( int x = 0 ; x < 65536 ; x++ ) i[x] = -1;
			for ( int x = 0 ; x <= 9 ; x++ ) for ( int y = 0 ; y <= 9 ; y++ ) i[x + y * 256 + 12336] = x * 10 + y;
#endif
			for ( int x = 0 ; x < 10000 ; x++ ) for ( int y = 3, z = x ; ~y ; y-- ) o[x * 4 + y] = z % 10 + 48, z /= 10;
		}
	} _table;
	template < class T, int digit > inline constexpr T _pw10 = 10 * _pw10 < T, digit - 1 >;
	template < class T > inline constexpr T _pw10 < T, 0 > = 1;
	inline void read(char& c) {
		do c = _gc();
		while ( !_isgraph(c) );
	}
	inline void read_cstr(char* s) {
		char c = _gc();
		while ( !_isgraph(c) ) c = _gc();
		while ( _isgraph(c) ) *s++ = c, c = _gc();
		*s = 0;
	}
	inline void read(string& s) {
		char c = _gc();
		s.clear();
		while ( !_isgraph(c) ) c = _gc();
		while ( _isgraph(c) ) s.push_back(c), c = _gc();
	}
	template < class T, bool neg >
#ifndef LOCAL
	__attribute__((no_sanitize("undefined")))
#endif
	inline void _read_int_suf(T& x) {
		_chk_i();
		char c;
		while
#ifndef LOCAL
		( ~_table.i[*reinterpret_cast < unsigned short*& >(_read_ptr)] ) if constexpr ( neg ) x = x * 100 - _table.i[*reinterpret_cast < unsigned short*& >(_read_ptr)++];
			else x = x * 100 + _table.i[*reinterpret_cast < unsigned short*& >(_read_ptr)++];
		if
#endif
		( _isdigit(c = _gc_nochk()) ) if constexpr ( neg ) x = x * 10 - ( c & 15 );
			else x = x * 10 + ( c & 15 );
	}
	template < class T, enable_if_t < _is_signed < T >, int > = 0 > inline void read(T& x) {
		char c;
		while ( !_isdigit(c = _gc()) ) if ( c == 45 ) {
				_read_int_suf < T, true >(x = -( _gc_nochk() & 15 ));
				return;
			}
		_read_int_suf < T, false >(x = c & 15);
	}
	template < class T, enable_if_t < _is_unsigned < T >, int > = 0 > inline void read(T& x) {
		char c;
		while ( !_isdigit(c = _gc()) );
		_read_int_suf < T, false >(x = c & 15);
	}
	inline void write(bool x) {
		_pc(x | 48);
	}
	inline void write(char c) {
		_pc(c);
	}
	inline void write_cstr(const char* s) {
		while ( *s ) _pc(*s++);
	}
	inline void write(const string& s) {
		for ( char c : s ) _pc(c);
	}
	template < class T, bool neg, int digit > inline void _write_int_suf(T x) {
		if constexpr ( digit == 4 ) _pnc_nochk < 4 >(_table.o + ( neg ? -x : x ) * 4);
		else _write_int_suf < T, neg, digit / 2 > (x / _pw10 < T, digit / 2 > ), _write_int_suf < T, neg, digit / 2 > (x % _pw10 < T, digit / 2 > );
	}
	template < class T, bool neg, int digit > inline void _write_int_pre(T x) {
		if constexpr ( digit <= 4 ) if ( digit >= 3 && ( neg ? x <= -100 : x >= 100 ) ) if ( digit >= 4 && ( neg ? x <= -1000 : x >= 1000 ) ) _pnc_nochk < 4 >(_table.o + ( neg ? -x : x ) * 4);
				else _pnc_nochk < 3 >(_table.o + ( neg ? -x : x ) * 4 + 1);
			else if ( digit >= 2 && ( neg ? x <= -10 : x >= 10 ) ) _pnc_nochk < 2 >(_table.o + ( neg ? -x : x ) * 4 + 2);
			else _pc_nochk(( neg ? -x : x ) | 48);
		else {
			constexpr int cur = 1 << __lg(digit - 1);
			if ( neg ? x <= -_pw10 < T, cur > : x >= _pw10 < T, cur > ) _write_int_pre < T, neg, digit - cur > (x / _pw10 < T, cur >), _write_int_suf < T, neg, cur >(x % _pw10 < T, cur >);
			else _write_int_pre < T, neg, cur >(x);
		}
	}
	template < class T, enable_if_t < _is_signed < T >, int > = 0 > inline void write(T x) {
		if ( x >= 0 ) _write_int_pre < T, false, numeric_limits < T >::digits10 + 1 > (x);
		else _pc_nochk(45), _write_int_pre < T, true, numeric_limits < T >::digits10 + 1 > (x);
		_chk_o();
	}
	template < class T, enable_if_t < _is_unsigned < T >, int > = 0 > inline void write(T x) {
		_write_int_pre < T, false, numeric_limits < T >::digits10 + 1 > (x), _chk_o();
	}
	template < size_t N, class ...T > inline void _read_tuple(tuple < T... >& x) {
		read(get < N >(x));
		if constexpr ( N + 1 != sizeof...(T) ) _read_tuple < N + 1, T... > (x);
	}
	template < size_t N, class ...T > inline void _write_tuple(const tuple < T... >& x) {
		write(get < N >(x));
		if constexpr ( N + 1 != sizeof...(T) ) _pc(32), _write_tuple < N + 1, T... > (x);
	}
	template < class ...T > inline void read(tuple < T... >& x) {
		_read_tuple < 0, T... >(x);
	}
	template < class ...T > inline void write(const tuple < T... >& x) {
		_write_tuple < 0, T... >(x);
	}
	template < class T1, class T2 > inline void read(pair < T1, T2 >& x) {
		read(x.first), read(x.second);
	}
	template < class T1, class T2 > inline void write(const pair < T1, T2 >& x) {
		write(x.first), _pc(32), write(x.second);
	}
	template < class T > inline auto read(T& x) -> decltype(x.read(), void()) {
		x.read();
	}
	template < class T > inline auto write(const T& x) -> decltype(x.write(), void()) {
		x.write();
	}
	template < class T1, class ...T2 > inline void read(T1& x, T2& ...y) {
		read(x), read(y...);
	}
	template < class ...T > inline void read_cstr(char* x, T* ...y) {
		read_cstr(x), read_cstr(y...);
	}
	template < class T1, class ...T2 > inline void write(const T1& x, const T2& ...y) {
		write(x), write(y...);
	}
	template < class ...T > inline void write_cstr(const char* x, const T* ...y) {
		write_cstr(x), write_cstr(y...);
	}
	template < class T > inline void print(const T& x) {
		write(x);
	}
	inline void print_cstr(const char* x) {
		write_cstr(x);
	}
	template < class T1, class ...T2 > inline void print(const T1& x, const T2& ...y) {
		write(x), _pc(32), print(y...);
	}
	template < class ...T > inline void print_cstr(const char* x, const T* ...y) {
		write_cstr(x), _pc(32), print_cstr(y...);
	}
	inline void println() {
		_pc(10);
	}
	inline void println_cstr() {
		_pc(10);
	}
	template < class ...T > inline void println(const T& ...x) {
		print(x...), _pc(10);
	}
	template < class ...T > inline void println_cstr(const T* ...x) {
		print_cstr(x...), _pc(10);
	}
}
using FastIO::read, FastIO::read_cstr, FastIO::write, FastIO::write_cstr, FastIO::println, FastIO::println_cstr;
#pragma GCC diagnostic pop
#include <bits/stdc++.h>
using namespace std;
#define POLY_ALLOW_AVX2
template <uint32_t mod> class Z_32 {
	private:
		static constexpr uint32_t get_r() {
			uint32_t ret = mod;
			for (int i = 0; i < 4; i++) ret *= 2 - mod * ret;
			return ret;
		}
		static constexpr uint32_t r = get_r();
		static constexpr uint32_t n2 = -uint64_t(mod) % mod;
		static_assert(r * mod == 1 && mod < (1 << 30)&& mod & 1);
		uint32_t v;
	public:
		constexpr Z_32(): v(0) {}
		template <class int_t> constexpr Z_32(const int_t x): v(reduce(uint64_t((sizeof(int_t) < sizeof(uint32_t) ? x : x % int_t(mod)) + mod) * n2)) {};
		static constexpr inline uint32_t reduce(const uint64_t x) {
			return (x + uint64_t(uint32_t(x) * (-r)) * mod) >> 32;
		}
		constexpr inline Z_32& operator += (const Z_32& r) {
			if (int32_t(v += r.v - 2 * mod) < 0) {
				v += 2 * mod;
			}
			return *this;
		}
		constexpr inline Z_32& operator -= (const Z_32& r) {
			if (int32_t(v -= r.v) < 0) {
				v += 2 * mod;
			}
			return *this;
		}
		constexpr inline Z_32& operator *= (const Z_32& r) {
			return v = reduce((uint64_t)v * r.v), *this;
		}
		constexpr inline Z_32& operator /= (const Z_32& r) {
			return *this *= r.inv();
		}
		constexpr inline friend Z_32 operator + (Z_32 l, const Z_32& r) {
			return l += r;
		}
		constexpr inline friend Z_32 operator - (Z_32 l, const Z_32& r) {
			return l -= r;
		}
		constexpr inline friend Z_32 operator * (Z_32 l, const Z_32& r) {
			return l *= r;
		}
		constexpr inline friend Z_32 operator / (Z_32 l, const Z_32& r) {
			return l /= r;
		}
		constexpr inline Z_32 operator - () const {
			return Z_32() - Z_32(*this);
		}
		template <class int_t> constexpr Z_32 pow(int_t r) const {
			Z_32 res(1), w(*this);
			for (; r; r >>= 1, w *= w) if (r & 1) res *= w;
			return res;
		}
		constexpr inline Z_32 inv() const {
			return pow(mod - 2);
		}
		constexpr inline uint32_t value() const {
			uint32_t res = reduce(v);
			return res >= mod ? res - mod : res;
		}
		constexpr inline Z_32 div_2() const {
			Z_32 res = *this;
			return res.v & 1 && (res.v += mod), res.v >>= 1, res;
		}
		constexpr inline Z_32 shrink() const {
			Z_32 res = *this;
			return res.v >= mod && (res.v -= mod), res;
		}
};
#ifdef POLY_ALLOW_AVX2
#pragma GCC target("avx2")
#include <immintrin.h>
template <uint32_t mod> union Z_8x32 {
private:
	static constexpr uint32_t get_r() {
		uint32_t ret = mod;
		for (int i = 0; i < 4; i++) ret *= 2 - mod * ret;
		return ret;
	}
	static constexpr uint32_t r = get_r();
	Z_32<mod> v[8];
	__m256i x;
	static inline __m256i reduce_4x64(const __m256i x) {
		__m256i v1 = _mm256_mul_epu32(x, _mm256_set1_epi32(- r));
		__m256i v2 = _mm256_mul_epu32(v1, _mm256_set1_epi32(mod));
		return _mm256_add_epi64(x, v2);
	}
	inline Z_8x32(__m256i x): x(x) {}
public:
	inline Z_8x32(Z_32<mod> v = 0): x(_mm256_set1_epi32(*(uint32_t*) & v)) {}
	inline Z_8x32(Z_32<mod> v0, Z_32<mod> v1): x(_mm256_set_epi32(*(uint32_t*) & v1, *(uint32_t*) & v0, *(uint32_t*) & v1, *(uint32_t*) & v0, *(uint32_t*) & v1, *(uint32_t*) & v0, *(uint32_t*) & v1, *(uint32_t*) & v0)) {}
	inline Z_8x32(Z_32<mod> v0, Z_32<mod> v1, Z_32<mod> v2, Z_32<mod> v3): x(_mm256_set_epi32(*(uint32_t*) & v3, *(uint32_t*) & v2, *(uint32_t*) & v1, *(uint32_t*) & v0, *(uint32_t*) & v3, *(uint32_t*) & v2, *(uint32_t*) & v1, *(uint32_t*) & v0)) {}
	inline Z_8x32(Z_32<mod> v0, Z_32<mod> v1, Z_32<mod> v2, Z_32<mod> v3, Z_32<mod> v4, Z_32<mod> v5, Z_32<mod> v6, Z_32<mod> v7): x(_mm256_set_epi32(*(uint32_t*) & v7, *(uint32_t*) & v6, *(uint32_t*) & v5, *(uint32_t*) & v4, *(uint32_t*) & v3, *(uint32_t*) & v2, *(uint32_t*) & v1, *(uint32_t*) & v0)) {}
	inline static Z_8x32 load(void* p) {
		return Z_8x32(_mm256_loadu_si256((__m256i*)p));
	}
	inline void store(void* p) const {
		_mm256_storeu_si256((__m256i*)p, x);
	}
	inline static void shuffle_2x2x4(Z_8x32& l, Z_8x32& r) {
		__m256i u = _mm256_permute2x128_si256(l.x, r.x, 0x20);
		__m256i v = _mm256_permute2x128_si256(l.x, r.x, 0x31);
		l = u, r = v;
	}
	inline static void shuffle_2x4x2(Z_8x32& l, Z_8x32& r) {
		l.x = _mm256_permute4x64_epi64(l.x, 0xD8);
		r.x = _mm256_permute4x64_epi64(r.x, 0xD8);
	}
	Z_32<mod>& operator [] (int idx) {
		return v[idx];
	}
	uint64_t& u64(int idx) {
		return idx[(uint64_t*)this];
	}
	inline Z_8x32 operator + (const Z_8x32& o) const {
		__m256i tmp = _mm256_add_epi32(x, o.x), m2x8 = _mm256_set1_epi32(mod * 2);
		return Z_8x32(_mm256_min_epu32(tmp, _mm256_sub_epi32(tmp, m2x8)));
	}
	inline Z_8x32 operator - (const Z_8x32& o) const {
		__m256i tmp = _mm256_sub_epi32(x, o.x), m2x8 = _mm256_set1_epi32(mod * 2);
		return Z_8x32(_mm256_min_epu32(tmp, _mm256_add_epi32(tmp, m2x8)));
	}
	inline Z_8x32 operator * (const Z_8x32& o) const {
		__m256i lo = _mm256_srli_epi64(reduce_4x64(_mm256_mul_epu32(x, o.x)), 32);
		__m256i hi = reduce_4x64(_mm256_mul_epu32(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(o.x, 32)));
		return Z_8x32(_mm256_blend_epi32(lo, hi, 0xAA));
	}
	inline Z_8x32 operator * (const Z_32<mod>& o) const {
		uint32_t w;
		memcpy(&w, (void*)&o, 4);
		__m256i v = _mm256_set1_epi32(w);
		__m256i lo = _mm256_srli_epi64(reduce_4x64(_mm256_mul_epu32(x, v)), 32);
		__m256i hi = reduce_4x64(_mm256_mul_epu32(_mm256_srli_epi64(x, 32), v));
		return Z_8x32(_mm256_blend_epi32(lo, hi, 0xAA));
	}
	inline Z_8x32 DFT_1() const {
		__m256i lo = _mm256_add_epi32(x, _mm256_srli_epi64(x, 32));
		__m256i hi = _mm256_sub_epi32(_mm256_slli_epi64(x, 32), x);
		__m256i u = _mm256_blend_epi32(lo, hi, 0xAA);
		__m256i v = _mm256_add_epi32(u, _mm256_set1_epi64x(uint64_t(mod * 2) << 32 | -mod * 2));
		return Z_8x32(_mm256_min_epu32(u, v));
	}
	inline Z_8x32 DIF_small(Z_8x32<mod> w2, Z_8x32<mod> w4) const {
		Z_8x32 v(_mm256_permutevar8x32_epi32(x, _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0)));
		v = v.DFT_1().mul_hilohi(w4);
		v.x = _mm256_permutevar8x32_epi32(v.x, _mm256_set_epi32(7, 3, 5, 1, 6, 2, 4, 0));
		v = v.DFT_1().mul_hilohi(w2);
		return v.x = _mm256_shuffle_epi32(v.x, 0xD8), v.DFT_1();
	}
	inline Z_8x32 DIT_small(Z_8x32<mod> w2, Z_8x32<mod> w4) const {
		Z_8x32 v = DFT_1();
		v.x = _mm256_shuffle_epi32(v.x, 0xD8);
		v = v.mul_hilohi(w2).DFT_1();
		v.x = _mm256_permutevar8x32_epi32(v.x, _mm256_set_epi32(7, 3, 5, 1, 6, 2, 4, 0));
		v = v.mul_hilohi(w4).DFT_1();
		return v.x = _mm256_permutevar8x32_epi32(v.x, _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0)), v;
	}
	inline Z_8x32 mul_hi(const Z_8x32& o) const {
		__m256i hi = reduce_4x64(_mm256_mul_epu32(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(o.x, 32)));
		return Z_8x32(_mm256_blend_epi32(x, hi, 0xAA));
	}
	inline Z_8x32 mul_hilohi(const Z_8x32& o) const {
		__m256i hi = reduce_4x64(_mm256_mul_epu32(_mm256_srli_epi64(x, 32), o.x));
		return Z_8x32(_mm256_blend_epi32(x, hi, 0xAA));
	}
	inline static Z_8x32 sbm(const Z_8x32& a, const Z_8x32& b, const Z_8x32& c) {
		__m256i tmp = _mm256_sub_epi32(a.x, b.x), m2x8 = _mm256_set1_epi32(mod * 2);
		return Z_8x32(_mm256_add_epi32(tmp, m2x8)) * c;
	}
private:
	void debug() const {
		for (size_t i = 0; i < 8; i++) cerr << v[i].value() << ' ';
		cerr << endl;
	}
};
#endif
template < uint32_t mod, uint32_t g, size_t LIM = 1 << __builtin_ctz(mod - 1) > class Poly_32 {
private:
	static constexpr size_t UNROLL_LOOP_LEVEL = 16;
	union SIMD {
		Z_32<mod> v[LIM];
#ifdef POLY_ALLOW_AVX2
		Z_8x32<mod> v8[LIM / 8];
#endif
		SIMD() {}
		Z_32<mod>& operator [] (size_t idx) {
			return v[idx];
		}
		const Z_32<mod>& operator [] (size_t idx) const {
			return v[idx];
		}
	};
	static SIMD init_w0() {
		SIMD w0;
		w0[LIM >> 1] = 1;
		Z_32<mod> w = omega(LIM)      ;
		for (size_t i = 1; i < (LIM >> 1); i++) w0[(LIM >> 1) + i] = w0[(LIM >> 1) + i - 1] * w;
		for (size_t i = (LIM >> 1) - 1; i; i--) w0[i] = w0[i << 1];
		for (size_t i = 0; i < (LIM >> 1); i++) w0[i] = w0[i].shrink();
		return w0;
	}
	static SIMD init_w1() {
		SIMD w1;
		w1[LIM >> 1] = 1;
		Z_32<mod> w = omega(LIM).inv();
		for (size_t i = 1; i < (LIM >> 1); i++) w1[(LIM >> 1) + i] = w1[(LIM >> 1) + i - 1] * w;
		for (size_t i = (LIM >> 1) - 1; i; i--) w1[i] = w1[i << 1];
		for (size_t i = 0; i < (LIM >> 1); i++) w1[i] = w1[i].shrink();
		return w1;
	}
	static inline SIMD butterfly(const SIMD& w) {
		SIMD w_;
		size_t id[LIM >> 1];
		id[0] = 0;
		for (size_t i = 0; i < (LIM >> 1); i++) id[i] = (id[i >> 1] | (i & 1) * LIM) >> 1;
		for (size_t k = 1; k < LIM; k <<= 1) for (size_t i = 0; i < k; i++) w_[k + i] = w[k + id[i] / (LIM / k)];
		return w_;
	}
	static inline SIMD w0 = init_w0(), w1 = init_w1();
	static inline SIMD w0_ = butterfly(w0), w1_ = butterfly(w1);
	static inline Z_32<mod> fac[LIM], inv_fac[LIM], inv[LIM];
	static inline int _ = []() {
		fac[0] = 1;
		for (size_t i = 1; i < LIM; i++) fac[i] = fac[i - 1] * i;
		inv_fac[LIM - 1] = fac[LIM - 1].inv();
		for (size_t i = LIM - 1; i; i--) inv_fac[i - 1] = inv_fac[i] * i;
		for (size_t i = 1; i < LIM; i++) inv[i] = inv_fac[i] * fac[i - 1];
		return 0;
	}();
	static constexpr Z_32<mod> omega(size_t n) {
		return (mod - 1) % n ? throw domain_error("") : Z_32<mod>(g).pow((mod - 1) / n);
	}
	template <size_t j, auto f, SIMD& w, class It> static inline void DFT_small_0(It a, size_t i, size_t k, Z_32<mod>& u, Z_32<mod>& v) {
		if constexpr(j) DFT_small_0 < j - 1, f, w > (a, i, k, u, v), f(a[i + j - 1], a[i + j + k - 1], u, v, w[j + k - 1]);
	}
	template <size_t j, auto f, SIMD& w, class It> static inline void DFT_small(It a, size_t i, size_t j0, size_t k, Z_32<mod>& u, Z_32<mod>& v) {
		if constexpr(j) DFT_small < j - 1, f, w > (a, i, j0, k, u, v), f(a[i + j0 + j - 1], a[i + j0 + j + k - 1], u, v, w[j + j0 + k - 1]);
	}
	template <size_t k, auto f, SIMD& w, class It> static inline void DFT_small_loop(It a, size_t n) {
		for (size_t i = 0; i + (k << 1) <= n; i += k << 1) {
			Z_32<mod> u, v;
			DFT_small_0<k, f, w>(a, i, k, u, v);
		}
	}
	template <size_t k, auto f, class It> static inline void DIF_small_loops(It a, size_t n) {
		if constexpr(k) DFT_small_loop<k, f, w0>(a, n), DIF_small_loops < (k >> 1), f > (a, n);
	}
	template <size_t k, auto f, class It> static inline void DIT_small_loops(It a, size_t n) {
		if constexpr(k) DIT_small_loops < (k >> 1), f > (a, n), DFT_small_loop<k, f, w1>(a, n);
	}
public:
	template <class It> static void DIF(It a, size_t n) {
		auto f = [](Z_32<mod>& x, Z_32<mod>& y, Z_32<mod>& u, Z_32<mod>& v, Z_32<mod> w) {
			u = x, v = y, x = (u + v), y = (u - v) * w;
		};
#ifdef POLY_ALLOW_AVX2
		if (n < 8) return DIF_small_loops<4, f>(a, n);
		for (size_t k = n >> 1; k > 4; k >>= 1) {
			for (size_t i = 0; i < n; i += k << 1) for (size_t j = 0; j < k; j += 8) {
					Z_8x32<mod> u = Z_8x32<mod>::load(a + i + j    );
					Z_8x32<mod> v = Z_8x32<mod>::load(a + i + j + k);
					Z_8x32<mod> x = u + v, y = Z_8x32<mod>::sbm(u, v, w0.v8[(j + k) >> 3]);
					x.store(a + i + j    );
					y.store(a + i + j + k);
				}
		}
		Z_8x32<mod> w2(1, 0, w0[3], 0, 1, 0, w0[3], 0), w4(1, 0, w0[5], 0, w0[6], 0, w0[7], 0);
		for (size_t i = 0; i < n; i += 8) Z_8x32<mod>::load(a + i).DIF_small(w2, w4).store(a + i);
#else
		for (size_t k = n >> 1; k > UNROLL_LOOP_LEVEL; k >>= 1) for (size_t i = 0; i < n; i += k << 1) {
				Z_32<mod> u, v;
				for (size_t j = 0; j < k; j += UNROLL_LOOP_LEVEL)
					DFT_small<UNROLL_LOOP_LEVEL, f, w0>(a, i, j, k, u, v);
			}
		DIF_small_loops<UNROLL_LOOP_LEVEL, f>(a, n);
#endif
	}
	template <class It> static void DIT(It a, size_t n) {
		auto f = [](Z_32<mod>& x, Z_32<mod>& y, Z_32<mod>& u, Z_32<mod>& v, Z_32<mod> w) {
			u = x, v = y * w, x = (u + v), y = (u - v);
		};
#ifdef POLY_ALLOW_AVX2
		if (n < 8) return DIT_small_loops<4, f>(a, n), mul_c_n(a, -Z_32<mod>((mod - 1) >> __lg(n)), n, a);
		Z_8x32<mod> w2(1, 0, w1[3], 0, 1, 0, w1[3], 0), w4(1, 0, w1[5], 0, w1[6], 0, w1[7], 0);
		for (size_t i = 0; i < n; i += 8) Z_8x32<mod>::load(a + i).DIT_small(w2, w4).store(a + i);
		for (size_t k = 4 << 1; k < n; k <<= 1) {
			for (size_t i = 0; i < n; i += k << 1) for (size_t j = 0; j < k; j += 8) {
					Z_8x32<mod> u = Z_8x32<mod>::load(a + i + j    );
					Z_8x32<mod> v = Z_8x32<mod>::load(a + i + j + k);
					v = v * w1.v8[(j + k) >> 3];
					Z_8x32<mod> x = u + v, y = u - v;
					x.store(a + i + j    );
					y.store(a + i + j + k);
				}
		}
#else
		DIT_small_loops<UNROLL_LOOP_LEVEL, f>(a, n);
		for (size_t k = UNROLL_LOOP_LEVEL << 1; k < n; k <<= 1) for (size_t i = 0; i < n; i += k << 1) {
				Z_32<mod> u, v;
				for (size_t j = 0; j < k; j += UNROLL_LOOP_LEVEL)
					DFT_small<UNROLL_LOOP_LEVEL, f, w1>(a, i, j, k, u, v);
			}
#endif
		mul_c_n(a, -Z_32<mod>((mod - 1) >> __lg(n)), n, a);
	}
	template <class It0, class It1> static inline void neg_n(It0 a, size_t n, It1 b) {
		for (size_t i = 0; i < n; i++) b[i] = -a[i];
	}
	template <class It0, class It1, class It2> static inline void add_n(It0 a, It1 b, size_t n, It2 c) {
		for (size_t i = 0; i < n; i++) c[i] = a[i] + b[i];
	}
	template <class It0, class It1, class It2> static inline void sub_n(It0 a, It1 b, size_t n, It2 c) {
		for (size_t i = 0; i < n; i++) c[i] = a[i] - b[i];
	}
	template <class It0, class It1, class It2> static inline void dot_n(It0 a, It1 b, size_t n, It2 c) {
#ifdef POLY_ALLOW_AVX2
		for (size_t i = 0; i + 7 < n; i += 8) (Z_8x32<mod>::load(a + i) * Z_8x32<mod>::load(b + i)).store(c + i);
		for (size_t i = n / 8 * 8; i < n; i++) c[i] = a[i] * b[i];
#else
		for (size_t i = 0; i < n; i++) c[i] = a[i] * b[i];
#endif
	}
	template <class It0, class It1> static inline void mul_c_n(It0 a, Z_32<mod> b, size_t n, It1 c) {
#ifdef POLY_ALLOW_AVX2
		for (size_t i = 0; i + 7 < n; i += 8) (Z_8x32<mod>::load(a + i) * b).store(c + i);
		for (size_t i = n / 8 * 8; i < n; i++) c[i] = a[i] * b;
#else
		for (size_t i = 0; i < n; i++) c[i] = a[i] * b;
#endif
	}
	template <class It0, class It1> static inline void div_2_n(It0 a, size_t n, It1 b) {
		for (size_t i = 0; i < n; i++) b[i] = a[i].div_2();
	}
	template <class It0, class It1> static inline void comp_ax_n(It0 a, Z_32<mod> b, size_t n, It1 c) {
#ifdef POLY_ALLOW_AVX2
		Z_8x32<mod> w;
		w[0] = 1;
		for (size_t i = 1; i < 8; i++) w[i] = w[i - 1] * b;
		Z_8x32<mod> v(w[7] * b);
		for (size_t i = 0; i + 7 < n; i += 8) (Z_8x32<mod>::load(a + i) * w).store(c + i), w = w * v;
		for (size_t i = n / 8 * 8; i < n; i++) c[i] = a[i] * w[0], w[0] *= b;
#else
		Z_32<mod> w = 1;
		for (size_t i = 0; i < n; i++) c[i] = a[i] * w, w *= b;
#endif
	}
	template <class It0, class It1> static inline void comp_x2_2n(It0 a, size_t n, It1 b) {
		for (size_t i = n; i; i--) b[(i - 1) << 1] = a[i - 1], b[(i - 1) << 1 | 1] = 0;
	}
	template <class It0, class It1> static inline void der_implace_n(It0 a, size_t n, It1 b) {
#ifdef POLY_ALLOW_AVX2
		Z_8x32<mod> w, v(8);
		for (size_t i = 0; i < 8; i++) w[i] = i;
		for (size_t i = 0; i + 7 < n; i += 8) (Z_8x32<mod>::load(a + i) * w).store(b + i), w = w + v;
		for (size_t i = n / 8 * 8; i < n; i++) b[i] = a[i] * i;
#else
		for (size_t i = 0; i < n; i++) b[i] = a[i] * i;
#endif
	}
	template <class It0, class It1> static inline void int_implace_n(It0 a, size_t n, It1 b) {
#ifdef POLY_ALLOW_AVX2
		for (size_t i = _; i + 7 < n; i += 8) (Z_8x32<mod>::load(a + i) * Z_8x32<mod>::load(inv + i)).store(b + i);
		for (size_t i = n / 8 * 8; i < n; i++) b[i] = a[i] * inv[i];
#else
		for (size_t i = _; i < n; i++) b[i] = a[i] * inv[i];
#endif
	}
	template <class It0, class It1> static inline void der_n(It0 a, size_t n, It1 b) {
		for (size_t i = 1; i < n; i++) b[i - 1] = a[i] * i;
		b[n - 1] = 0;
	}
	template <class It0, class It1> static inline void int_n(It0 a, size_t n, It1 b) {
		for (size_t i = n - 1; i; i--) b[i] = a[i - 1] * inv[i];
		b[0] = _;
	}
	template <class It> static inline void DFT_high_n(It a, size_t n) {
		for (size_t i = 0; i < n; i++) a[i] *= w0[n + i];
		DIF(a, n);
	}
	template <class It> static inline void DFT_keep_0_halfn(It a, size_t n) {
		for (size_t i = 0; i < (n >> 1); i++) a[i] = a[i << 1] + a[i << 1 | 1];
		div_2_n(a, n >> 1, a);
	}
	template <class It> static inline void DFT_keep_1_halfn(It a, size_t n) {
		for (size_t i = 0; i < (n >> 1); i++) a[i] = (a[i << 1] - a[i << 1 | 1]) * w1_[(n >> 1) + i];
		div_2_n(a, n >> 1, a);
	}
	template <class It> static inline void DFT_mul_x_n(It a, size_t n) {
		for (size_t i = 0; i < n; i += 2) {
			Z_32<mod> w = w0_[(n + i) >> 1];
			a[i] *= w, a[i | 1] *= -w;
		}
	}
	template <class It> static inline void DFT_comp_x2_2n(It a, size_t n) {
		comp_x2_2n(a, n, a);
		for (size_t i = 0; i < n; i++) a[i << 1 | 1] = a[i << 1];
	}
	template <class It0, class It1, class It2, class It3 = Z_32<mod>*, class It4 = Z_32<mod>*> static inline void conv_n(It0 a, It1 b, size_t n, It2 c, It3 a_ = (Z_32<mod>*)nullptr, It4 b_ = (Z_32<mod>*)nullptr) {
		alignas(32) Z_32<mod> x[n];
		if (a_ == nullptr) copy_n(a, n, x), DIF(x, n);
		else copy_n(a_, n, x);
		if (b_ == nullptr) copy_n(b, n, c), DIF(c, n);
		else copy_n(b_, n, c);
		dot_n(x, c, n, c), DIT(c, n);
	}
private:
	template <class It0, class It1, class It2, class It3 = Z_32<mod>*, class It4 = Z_32<mod>*> static inline void inv_work_n(It0 a, It1 b, size_t n, It2 c, It3 a_ = (Z_32<mod>*)nullptr, It4 b_ = (Z_32<mod>*)nullptr) {
		alignas(32) Z_32<mod> x[n], y[n];
		if (a_ == nullptr) copy_n(a, n, x), DIF(x, n);
		else copy_n(a_, n, x);
		if (b_ == nullptr) copy_n(b, n, y), DIF(y, n);
		else copy_n(b_, n, y);
		dot_n(x, y, n, c), dot_n(c, y, n, c), sub_n(c, y, n, c), DIT(c, n);
	}
	template <class It0, class It1, class It2 = Z_32<mod>*, class It3 = Z_32<mod>*> static inline void inv_step_2n(It0 a, size_t n, It1 b, It2 a_ = (Z_32<mod>*)nullptr, It3 b0_ = (Z_32<mod>*)nullptr) {
		try {
			alignas(32) Z_32<mod> x[n << 1], y[n], w = omega(n << 2), i = omega(4);
			comp_ax_n(a, w, n << 1, x), comp_ax_n(b, w, n, y), add_n(x, x + n, n, x), inv_work_n(x, y, n, y), comp_ax_n(y, w.inv(), n, y);
			fill_n(b + n, n, 0), inv_work_n(a, b, n << 1, x, a_, b0_), add_n(x, y, n, y), mul_c_n(y, -i, n, y);
			add_n(x + n, y, n, x + n), div_2_n(x + n, n, b + n), neg_n(b + n, n, b + n);
		} catch (domain_error&) {
			alignas(32) Z_32<mod> x[n << 1], y[n << 1];
			if (b0_ == nullptr) copy_n(b, n, y), fill_n(y + n, n, 0), DIF(y, n << 1);
			else copy_n(b0_, n << 1, y);
			conv_n(b, a, n << 1, x, y, a_), fill_n(x, n, 0), conv_n(b, x, n << 1, x, y), neg_n(x + n, n, b + n);
		}
	}
public:
	template <class It0, class It1, class It2> static inline void mul_n(It0 a, It1 b, size_t n, It2 c) {
		alignas(32) Z_32<mod> x[n << 1];
		copy_n(a, n, x), fill_n(x + n, n, 0);
		alignas(32) Z_32<mod> y[n << 1];
		copy_n(b, n, y), fill_n(y + n, n, 0);
		conv_n(x, y, n << 1, x), copy_n(x, n, c);
	}
	template <class It0, class It1, class It2 = Z_32<mod>*> static inline void inv_n(It0 a, size_t n, It1 b, It2 a_ = (Z_32<mod>*)nullptr) {
		alignas(32) Z_32<mod> x[n];
		x[0] = a[0].value() == 1 ? 1 : a[0].inv();
		for (size_t k = 1; k < n; k <<= 1) inv_step_2n(a, k, x, (k << 1) == n ? a_ : nullptr);
		copy_n(x, n, b);
	}
	template <class It0, class It1> static inline void sqrt_n(It0 a, size_t n, It1 b) {
		alignas(32) Z_32<mod> x[n], xi[n], y[n], z[n], u[n], v[n];
		x[0] = xi[0] = 1, u[0] = x[0], u[1] = 0, DIF(u, 2);
		for (size_t k = 1; k < n; k <<= 1) {
			copy_n(xi, k, v), fill_n(v + k, k, 0), DIF(v, k << 1), copy_n(a, k << 1, z), DIF(z, k << 1);
			dot_n(u, u, k << 1, y), sub_n(y, z, k << 1, y), conv_n(xi, (Z_32<mod>*)nullptr, k << 1, y, v, y);
			div_2_n(y + k, k, y + k), neg_n(y + k, k, x + k);
			if ((k << 1) < n) {
				copy_n(x, k << 1, u), fill_n(u + (k << 1), k << 1, 0), DIF(u, k << 2);
				inv_step_2n(x, k, xi, u, v);
			}
		}
		copy_n(x, n, b);
	}
	template <class It0, class It1> static inline void log_n(It0 a, size_t n, It1 b) {
		alignas(32) Z_32<mod> x[n];
		der_n(a, n, x), inv_n(a, n, b), mul_n(x, b, n, b), int_n(b, n, b);
	}
	template <class It0, class It1, class It2 = Z_32<mod>*> static inline void exp_n(It0 a, size_t n, It1 b, It2 c = (Z_32<mod>*)nullptr) {
		alignas(32) Z_32<mod> x[n], xi[n], y[n], z[n], u[n], v[n];
		x[0] = xi[0] = 1, u[0] = x[0], u[1] = 0, DIF(u, 2);
		for (size_t k = 1; k < n; k <<= 1) {
			der_implace_n(a, k, y), fill_n(y + k, k, 0), DIF(y, k << 1), dot_n(y, u, k << 1, y);
			der_implace_n(x, k, z), fill_n(z + k, k, 0), DIF(z, k << 1), sub_n(z, y, k << 1, y);
			copy_n(xi, k, v), fill_n(v + k, k, 0), DIF(v, k << 1);
			conv_n(xi, (Z_32<mod>*)nullptr, k << 1, y, v, y), fill_n(y, k, 0);
			der_implace_n(a, k, z), add_n(y, z, k, y), int_implace_n(y, k << 1, y), sub_n(y, a, k << 1, y);
			conv_n(x, y, k << 1, y, u), neg_n(y + k, k, x + k);
			if ((k << 1) < n || c) {
				copy_n(x, k << 1, u), fill_n(u + (k << 1), k << 1, 0), DIF(u, k << 2);
				inv_step_2n(x, k, xi, u, v);
			}
		}
		copy_n(x, n, b);
		if (c) copy_n(xi, n, c);
	}
	template <class It0, class Int, class It1> static inline void pow_n(It0 a, Int c, size_t n, It1 b) {
		if (!c) {
			fill_n(a, n, 0), a[0] = 1;
			return;
		}
		size_t ord = 0;
		while (ord < n && !a[ord].value()) ord++;
		if (ord == n || (ord && c > (n / ord))) {
			fill_n(a, n, 0);
			return;
		}
		Z_32<mod> low = a[ord];
		rotate(a, a + ord, a + n), mul_c_n(a, low.inv(), n, b);
		log_n(b, n, b), mul_c_n(b, c % mod, n, b), exp_n(b, n, b);
		mul_c_n(a, low.pow(c), n, b), ord *= c;
		for (size_t i = n; i > ord; i--) a[i - 1] = a[i - ord - 1];
		fill_n(a, ord, 0);
	}
	template <class It0, class It1> static inline void sin_n(It0 a, size_t n, It1 b) {
		alignas(32) Z_32<mod> x[n], w = omega(4);
		mul_c_n(a, w, n, x), exp_n(x, n, x, b);
		sub_n(x, b, n, b), mul_c_n(b, (w * 2).inv(), n, b);
	}
	template <class It0, class It1> static inline void cos_n(It0 a, size_t n, It1 b) {
		alignas(32) Z_32<mod> x[n], w = omega(4);
		mul_c_n(a, w, n, x), exp_n(x, n, x, b);
		add_n(x, b, n, b), div_2_n(b, n, b);
	}
	template <class It0, class It1> static inline void tan_n(It0 a, size_t n, It1 b) {
		alignas(32) Z_32<mod> x[n], y[n], w = omega(4);
		mul_c_n(a, w, n, x), exp_n(x, n, x, y);
		sub_n(x, y, n, b), add_n(x, y, n, x), inv_n(x, n, x), mul_n(x, b, n, b);
		mul_c_n(b, w.inv(), n, b);
	}
	template <class It0, class It1> static inline void asin_n(It0 a, size_t n, It1 b) {
		alignas(32) Z_32<mod> x[n];
		sqr_n(a, n, x), mul_c_n(x, -1, n, x), x[0] += 1, sqrt_n(x, n, x);
		der_n(a, n, b), inv_n(x, n, x), mul_n(x, b, n, b), int_n(b, n, b);
	}
	template <class It0, class It1> static inline void atan_n(It0 a, size_t n, It1 b) {
		alignas(32) Z_32<mod> x[n];
		sqr_n(a, n, x), x[0] += 1;
		der_n(a, n, b), inv_n(x, n, x), mul_n(x, b, n, b), int_n(b, n, b);
	}
private:
	template <class It0, class It1> static inline void comp_step_nxm(It0 a, It1 b, size_t n, size_t m) {
		alignas(32) Z_32<mod> x[n * m << 1], y[n * m << 1];
		for (size_t i = _; i < n; i++) copy_n(b + i * m, m, x + (i * m << 1)), fill_n(x + (i * m << 1) + m, m, 0);
		DIF(x, n * m << 1);
		for (size_t i = 0; i < (n * m << 1); i += 2) y[i] = x[i | 1], y[i | 1] = x[i];
		neg_n(y, n * m << 1, y);
		if (m > 2) {
			dot_n(x, y, n * m << 1, x), DFT_keep_0_halfn(x, n * m << 1), DIT(x, n * m);
			for (size_t i = 0; i < n; i++) copy_n(x + i * m + (m >> 1), m >> 1, x + (i * m >> 1));
			x[(m - 1) >> 1] -= 1;
			fill_n(x + (n * m >> 1), n * m >> 1, 0), x[(n * m + m - 1) >> 1] = 1;
			comp_step_nxm(a, x, n << 1, m >> 1), fill_n(x, n * m, 0);
		}
		for (size_t i = 0; i < n; i++) copy_n(a + (i * m >> 1), m >> 1, x + i * m);
		DIF(x, n * m), DFT_comp_x2_2n(x, n * m), DFT_mul_x_n(x, n * m << 1);
		dot_n(x, y, n * m << 1, y), DIT(y, n * m << 1);
		for (size_t i = 0; i < n; i++) copy_n(y + n * m + (i * m << 1) + m - 1, m, a + (i * m));
	}
public:
	template <class It0, class It1, class It2> static inline void comp_n(It0 a, It1 b, size_t n, It2 c) {
		alignas(32) Z_32<mod> x[n << 1], y[n << 1];
		fill_n(x + n, n, 0), fill_n(y + n, n, 0);
		for (size_t i = 0; i < n; i++) y[i] = fac[n + i - 1] * inv_fac[i] * inv_fac[n - 1];
		comp_ax_n(y, b[0], n, y), reverse_copy(a, a + n, x), mul_n(x, y, n, x), reverse(x, x + n);
		neg_n(b, n, y), reverse(y, y + n), y[(n << 1) - 1] = 1;
		comp_step_nxm(x, y, 2, n), reverse_copy(x, x + n, c);
	}
	template <class It0, class It1> static inline void comp_inv_n(It0 a, size_t n, It1 b) {
		alignas(32) Z_32<mod> v = a[1].inv(), x[n << 2], y[n << 2], z[n << 2];
		fill_n(x, n << 2, 0), fill_n(y, n << 2, 0), x[0] = y[0] = 1, mul_c_n(a, -v, n, y + n);
		for (size_t m = n, n = 2; m > 1; n <<= 1, m >>= 1) {
			for (size_t i = n; i; i--) {
				copy_n(x + (i - 1) * m, m, x + ((i - 1) << 1) * m);
				fill_n(x + ((i - 1) << 1 | 1) * m, m, 0);
				copy_n(y + (i - 1) * m, m, y + ((i - 1) << 1) * m);
				fill_n(y + ((i - 1) << 1 | 1) * m, m, 0);
			}
			DIF(x, n * m << 1), DIF(y, n * m << 1), copy_n(y, n * m << 1, z);
			for (size_t i = 0; i < n * m; i++) swap(z[i << 1], z[i << 1 | 1]);
			dot_n(x, z, n * m << 1, x), DFT_keep_1_halfn(x, n * m << 1), DIT(x, n * m);
			dot_n(y, z, n * m << 1, y), DFT_keep_0_halfn(y, n * m << 1), DIT(y, n * m);
			for (size_t i = 0; i < n; i++) {
				copy_n(x + i * m, (m >> 1), x + i * (m >> 1));
				copy_n(y + i * m, (m >> 1), y + i * (m >> 1));
			}
			fill_n(x + (n * m >> 1), n * m >> 1, 0);
			fill_n(y + (n * m >> 1), n * m >> 1, 0);
			copy(y + 1, y + (m >> 1), y + (n * m >> 1) + 1);
			fill(y + 1, y + (m >> 1), 0);
		}
		reverse_copy(x + 1, x + n, b);
		for (size_t i = 0; i < n - 1; i++) b[i] = b[i] * (n - 1) / (n - i - 1);
		comp_ax_n(b, v, n, b), pow_n(b, (-Z_32<mod>(n - 1).inv()).value(), n, b);
		mul_c_n(b, v, n, b), rotate(b, b + n - 1, b + n), b[0] = 0;
	}
private:
	template <class It0, class It1, class It2> static inline void multipoint_step_n(It0 a, It1 b, size_t n, size_t m, It2 b_) {
		alignas(32) Z_32<mod> v[n << 1];
		for (size_t i = 0; i < (n << 1); i += m << 1) {
			dot_n(b_ + i, b_ + i + m, m, v + i), copy_n(v + i, m, b + i), DIT(b + i, m);
			b[i] -= 1, fill_n(b + i + m, m, 0), b[i + m] = 1;
			if (m < n) copy_n(b + i, m, v + i + m), v[i + m] -= 1, DFT_high_n(v + i + m, m);
		}
		if (m < n) multipoint_step_n(a, b, n, m << 1, v);
		else {
			reverse(b, b + n + 1), inv_n(b, n, b), fill_n(b + n, n, 0), reverse(b, b + (n << 1));
			conv_n(a, b, n << 1, a), rotate(a, a + (n << 1) - 1, a + (n << 1));
		}
		for (size_t i = 0; i < (n << 1); i += m << 1) {
			DIF(a + i, m);
			dot_n(a + i, b_ + i, m, a + i + m), DIT(a + i + m, m);
			dot_n(a + i, b_ + i + m, m, a + i    ), DIT(a + i, m);
			copy_n(a + i + m + (m >> 1), m >> 1, a + i + m);
			copy_n(a + i     + (m >> 1), m >> 1, a + i    );
		}
	}
public:
	template <class It0, class It1, class It2> static inline void multipoint_n(It0 a, It1 b, size_t n, It2 c) {
		alignas(32) Z_32<mod> x[n << 1], y[n << 1], z[n << 1];
		for (size_t i = 0; i < n; i++) {
			y[i << 1] = -b[i], y[i << 1 | 1] = 1;
			copy_n(y + (i << 1), 2, z + (i << 1)), DIF(z + (i << 1), 2);
		}
		copy_n(a, n, x), fill_n(x + n, n, 0), multipoint_step_n(x, y, n, 2, z);
		for (size_t i = 0; i < n; i++) c[i] = x[i << 1];
	}
	template <class It0, class It1, class It2> static inline void interpolate_n(It0 a, It1 b, size_t n, size_t m, It2 c) {
		size_t lim = __lg(n);
		alignas(32) Z_32<mod> x[lim][n << 1], y[n << 1];
		for (size_t i = 0; i < n; i++) {
			if (i < m) x[0][i << 1] = -a[i], x[0][i << 1 | 1] = 1;
			else x[0][i << 1] = 1, x[0][i << 1 | 1] = 0;
			DIF(x[0] + (i << 1), 2);
		}
		for (size_t t = 1, k = 2; k <  n; t++, k <<= 1) for (size_t i = 0; i < (n << 1); i += k << 1) {
				dot_n(x[t - 1] + i, x[t - 1] + i + k, k, x[t] + i);
				copy_n(x[t] + i, k, x[t] + i + k), DIT(x[t] + i + k, k);
				if (i + (k << 1) <= (m << 1)) x[t][i + k] -= 2;
				DFT_high_n(x[t] + i + k, k);
			}
		dot_n(x[lim - 1], x[lim - 1] + n, n, y), DIT(y, n);
		if (n <= m) y[0] -= 1;
		der_n(y, n, y);
		if (n <= m) y[n - 1] = n;
		multipoint_n(y, a, n, y);
		for (size_t i = 0; i < n; i++) y[i] = b[i] * y[i].inv();
		DFT_comp_x2_2n(y, n);
		for (size_t t = 0, k = 2; k <= n; t++, k <<= 1) for (size_t i = 0; i < (n << 1); i += k << 1) {
				dot_n(y + i, x[t] + i + k, k, y + i), dot_n(y + i + k, x[t] + i, k, y + i + k), add_n(y + i, y + i + k, k, y + i);
				if (k < n) copy_n(y + i, k, y + i + k), DIT(y + i + k, k), DFT_high_n(y + i + k, k);
			}
		DIT(y, n), copy_n(y, n, c);
	}
};
