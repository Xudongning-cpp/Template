#pragma once
#define BIGINTEGER_H
#include <iomanip>
#include <cstring>
#include <climits>
#include <random>
#include <algorithm>
#include <functional>
#include <numeric>
#include <ranges>
class ZeroDivisionError : public std::runtime_error {
public:
    ZeroDivisionError() : std::runtime_error("Division by zero") {}
};
class FFTLimitExceededError : public std::runtime_error {
public:
    FFTLimitExceededError() : std::runtime_error("FFT size limit exceeded") {}
};
class NegativeRadicandError : public std::runtime_error {
public:
    NegativeRadicandError() : std::runtime_error("Negative radicand in root operation") {}
};
template<typename T>
concept IntegerLike = std::integral<T> || std::same_as<T, __int128>;
template<typename T>
concept SignedIntegerLike = std::signed_integral<T> || std::same_as<T, __int128>;
template<typename T>
concept UnsignedIntegerLike = std::unsigned_integral<T>;
class BigInteger {
protected:
    static constexpr int WIDTH = 8;
    static constexpr long long BASE = 1e8;
    static constexpr long long FFT_LIMIT = 32;
    static constexpr long long NEWTON_DIV_LIMIT = 64;
    static constexpr long long NEWTON_DIV_MIN_LEVEL = 16;
    static constexpr long long NEWTON_SQRT_LIMIT = 48;
    static constexpr long long NEWTON_SQRT_MIN_LEVEL = 5;
    static constexpr size_t POOL_CHUNK_SIZE = 2048;
    long long* digits;
    int capacity, size;
    bool flag;
    inline void push(const long long&);
    inline void pop();
    inline void resize(const int&);
    inline int compare(const BigInteger&) const;
    BigInteger& build_binary(const std::vector<bool>&);
    static inline BigInteger fft_mul(const BigInteger&, const BigInteger&);
    BigInteger newton_inv(int) const;
    std::pair<BigInteger, BigInteger> newton_div(const BigInteger&) const;
    BigInteger newton_invsqrt() const;
    BigInteger sqrt_normal() const;
    std::vector<long long*> memory_pool;
    inline long long* allocate_digits(size_t);
    inline void deallocate_digits(long long*, size_t);
public:
    void reserve(const int&);
    void clear();
    ~BigInteger() = default;
    int digit_length() const { return size; }
    bool is_positive() const { return flag; }
    bool is_negative() const { return !flag; }
    bool is_zero() const { return size == 1 && digits[1] == 0; }
    friend std::ostream& operator<<(std::ostream& out, const BigInteger& x) {
        if (!x.flag) out << '-';
        out << x.digits[x.size];
        for (int i = x.size - 1; i >= 1; i--)
            out << std::setw(WIDTH) << std::setfill('0') << x.digits[i];
        return out;
    }
    friend std::istream& operator>>(std::istream& in, BigInteger& x) {
        std::string s;
        in >> s;
        x = s;
        return in;
    }
    BigInteger() : digits(nullptr), capacity(0), size(0), flag(true) {}
    BigInteger(const BigInteger& x) : digits(nullptr) { *this = x; }
    template<IntegerLike T>
    BigInteger(T x) : digits(nullptr) { *this = x; }
    BigInteger(const std::string& s) : digits(nullptr) { *this = s; }
    BigInteger(const std::vector<bool>& b) : digits(nullptr) { *this = b; }
    template<std::input_iterator It>
    BigInteger(It begin, It end) : digits(nullptr) { *this = std::vector<bool>(begin, end); }
    BigInteger(const long long* in, int len) : digits(nullptr), flag(true) {
        while (len > 0 && in[len - 1] == 0) len--, in++;
        size = len;
        digits = allocate_digits(len + 1);
        for (int i = 0; i < len; i++) digits[i + 1] = in[i];
    }
    BigInteger& operator=(const BigInteger&);
    template<IntegerLike T>
    BigInteger& operator=(T x) {
        if constexpr (std::same_as<T, bool>) {
            *this = x ? 1LL : 0LL;
        } else {
            flag = (x >= 0);
            reserve(4);
            if (x == 0) return size = 1, digits[1] = 0, *this;
            if constexpr (std::same_as<T, __int128>) {
                if (x == __int128{0}) return size = 1, digits[1] = 0, *this;
                __int128 n = x >= 0 ? x : -x;
                do { push(n % BASE), n /= BASE; } while (n);
            } else {
                long long n = std::fabs(x);
                do { push(n % BASE), n /= BASE; } while (n);
            }
        }
        return *this;
    }
    BigInteger& operator=(const std::string&);
    BigInteger& operator=(const std::vector<bool>&);
    BigInteger& operator=(const __int128&);
    std::string to_string() const;
    long long to_long_long() const;
    std::vector<bool> to_binary() const;
    __int128 to_int128() const;
    std::string_view to_string_view() const noexcept {
        return std::string_view(to_string());
    }
    BigInteger operator-() const;
    BigInteger operator~() const;
    BigInteger abs() const;
    bool operator==(const BigInteger&) const;
    auto operator<=>(const BigInteger&) const;
    template<IntegerLike T>
    bool operator==(T x) const { return *this == BigInteger(x); }
    template<IntegerLike T>
    auto operator<=>(T x) const { return *this <=> BigInteger(x); }
    BigInteger div2() const;
    std::pair<BigInteger, BigInteger> divmod(const BigInteger&, bool = false) const;
    BigInteger operator+(const BigInteger&) const;
    BigInteger operator-(const BigInteger&) const;
    BigInteger operator*(const long long&) const;
    BigInteger operator*(const BigInteger&) const;
    BigInteger operator/(const long long&) const;
    BigInteger operator/(const BigInteger&) const;
    BigInteger operator%(const long long&) const;
    BigInteger operator%(const BigInteger&) const;
    BigInteger pow(const long long&) const;
    template <typename T>
    BigInteger pow(const long long&, const T&) const;
    BigInteger root(const long long& = 2) const;
    BigInteger sqrt() const;
    BigInteger gcd(const BigInteger&) const;
    BigInteger lcm(const BigInteger&) const;
    inline BigInteger _move_l(int) const;
    inline BigInteger _move_r(int) const;
    BigInteger& operator+=(const BigInteger&);
    BigInteger& operator-=(const BigInteger&);
    BigInteger& operator*=(long long);
    BigInteger& operator*=(const BigInteger&);
    BigInteger& operator/=(const long long&);
    BigInteger& operator/=(const BigInteger&);
    BigInteger& operator%=(const long long&);
    BigInteger& operator%=(const BigInteger&);
    BigInteger operator<<(const long long&);
    BigInteger operator>>(const long long&);
    BigInteger& operator<<=(const long long&);
    BigInteger& operator>>=(const long long&);
    BigInteger operator&(const BigInteger&);
    BigInteger operator|(const BigInteger&);
    BigInteger operator^(const BigInteger&);
    BigInteger& operator&=(const BigInteger&);
    BigInteger& operator|=(const BigInteger&);
    BigInteger& operator^=(const BigInteger&);
    BigInteger& operator++();
    BigInteger operator++(int);
    BigInteger& operator--();
    BigInteger operator--(int);
    static BigInteger random(const int&);
    auto as_span() const noexcept {
        return std::span<const long long>(digits + 1, size);
    }
    class iterator {
        const BigInteger* parent;
        int index;
    public:
        iterator(const BigInteger* p, int i) : parent(p), index(i) {}
        long long operator*() const { return parent->digits[index]; }
        iterator& operator++() { ++index; return *this; }
        iterator operator++(int) { auto tmp = *this; ++index; return tmp; }
        bool operator==(const iterator& other) const { return index == other.index; }
        bool operator!=(const iterator& other) const { return index != other.index; }
    };
    iterator begin() const { return iterator(this, 1); }
    iterator end() const { return iterator(this, size + 1); }
    auto keys() const {
        return std::views::iota(0, size) | std::views::transform([this](int i) { return digits[i+1]; });
    }
};
inline void BigInteger::push(const long long& val) {
    if (size == capacity) {
        int new_capacity = (capacity < 256) ? capacity * 2 : static_cast<int>(std::pow(capacity, 1.125));
        long long* new_digits = allocate_digits(new_capacity + 1);
        std::memcpy(new_digits, digits, sizeof(long long) * (capacity + 1));
        deallocate_digits(digits, capacity + 1);
        digits = new_digits;
        capacity = new_capacity;
    }
    digits[++size] = val;
}
inline void BigInteger::pop() {
    if (size > 0) digits[size--] = 0;
}
inline void BigInteger::resize(const int& sz) {
    reserve(sz), size = sz;
}
inline int BigInteger::compare(const BigInteger& x) const {
    if (flag && !x.flag) return 1;
    if (!flag && x.flag) return -1;
    int sgn = (flag && x.flag ? 1 : -1);
    if (size > x.size) return sgn;
    if (size < x.size) return -sgn;
    for (int i = size; i >= 1; i--) {
        if (digits[i] > x.digits[i]) return sgn;
        if (digits[i] < x.digits[i]) return -sgn;
    }
    return 0;
}
BigInteger& BigInteger::build_binary(const std::vector<bool>& b) {
    *this = 0ll;
    if (b.empty() || (b.size() == 1 && b[0] == 0)) return *this;
    BigInteger pw2 = 1;
    for (int i = b.size() - 1; i >= 0; i--, pw2 += pw2) if (b[i]) *this += pw2;
    return *this;
}
inline long long* BigInteger::allocate_digits(size_t size) {
    for (auto it = memory_pool.begin(); it != memory_pool.end(); ++it) {
        size_t* block_size = reinterpret_cast<size_t*>(*it) - 1;
        if (*block_size >= size) {
            long long* block = *it;
            memory_pool.erase(it);
            std::memset(block, 0, size * sizeof(long long));
            return block;
        }
    }
    size_t* mem = static_cast<size_t*>(::operator new(sizeof(size_t) + size * sizeof(long long)));
    *mem = size;
    long long* ptr = reinterpret_cast<long long*>(mem + 1);
    std::memset(ptr, 0, size * sizeof(long long));
    return ptr;
}
inline void BigInteger::deallocate_digits(long long* ptr, size_t size) {
    if (!ptr) return;
    size_t* block_size = reinterpret_cast<size_t*>(ptr) - 1;
    if (*block_size >= POOL_CHUNK_SIZE) memory_pool.push_back(ptr);
    else ::operator delete(block_size);
}
void BigInteger::reserve(const int& sz) {
    if (sz < 0) return;
    if (digits != nullptr) deallocate_digits(digits, capacity + 1);
    capacity = sz;
    size = 0;
    digits = allocate_digits(sz + 1);
}
void BigInteger::clear() {
    deallocate_digits(digits, capacity + 1);
    digits = nullptr;
    size = capacity = 0;
}
BigInteger& BigInteger::operator=(const BigInteger& x) {
    if (this != &x) {
        long long* new_digits = allocate_digits(x.size + 1);
        deallocate_digits(digits, capacity + 1);
        digits = new_digits;
        capacity = x.size;
        flag = x.flag;
        size = x.size;
        std::memcpy(digits, x.digits, sizeof(long long) * (x.size + 1));
    }
    return *this;
}
BigInteger& BigInteger::operator=(const std::string& s) {
    flag = true, reserve(s.size() / WIDTH + 1);
    if (s.empty() || s == "-") return *this = 0ll;
    int i = 0;
    if (s[0] == '-') flag = false, i++;
    for (int j = s.size() - 1; j >= i; j -= WIDTH) {
        int start = std::max(i, j - WIDTH + 1), len = j - start + 1;
        push(std::stoll(s.substr(start, len)));
    }
    while (size > 1 && digits[size] == 0) pop();
    return *this;
}
BigInteger& BigInteger::operator=(const std::vector<bool>& a) {
    if (a == std::vector<bool>{0ll}) return *this = 0ll;
    std::vector<bool> res;
    if (a[0] == 0) return build_binary(a);
    res.resize(a.size());
    for (int i = 0; i < static_cast<int>(a.size()); ++i) res[i] = !a[i];
    BigInteger x;
    x.build_binary(res);
    return *this = -x - 1;
}
std::string BigInteger::to_string() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}
std::vector<bool> BigInteger::to_binary() const {
    if (*this == 0) return {0, 0};
    std::vector<bool> res;
    if (flag) {
        for (BigInteger x = *this; x != 0; x = x.div2()) res.emplace_back(x.digits[1] & 1);
        res.emplace_back(0);
        std::reverse(res.begin(), res.end());
        return res;
    }
    else {
        for (BigInteger x = -*this - 1; x != 0; x = x.div2()) res.emplace_back(!(x.digits[1] & 1));
        res.emplace_back(1);
        std::reverse(res.begin(), res.end());
        return res;
    }
}
BigInteger BigInteger::operator-() const {
    if (*this == 0) return 0;
    BigInteger res = *this;
    res.flag = !flag;
    return res;
}
BigInteger BigInteger::operator~() const {
    return -(*this) - 1;
}
BigInteger BigInteger::abs() const {
    BigInteger res = *this;
    res.flag = true;
    return res;
}
bool BigInteger::operator==(const BigInteger& x) const { return compare(x) == 0; }
auto BigInteger::operator<=>(const BigInteger& x) const { return compare(x); }
BigInteger BigInteger::operator+(const BigInteger& x) const {
    if (!x.flag) return *this - x.abs();
    if (!flag) return x - abs();
    BigInteger res;
    res.flag = !(flag ^ x.flag);
    int n = std::max(size, x.size) + 1;
    res.reserve(n);
    long long carry = 0;
    for (int i = 1; i <= n; i++) {
        long long d1 = i <= size ? digits[i] : 0, d2 = i <= x.size ? x.digits[i] : 0;
        res.push(d1 + d2 + carry);
        carry = res.digits[i] / BASE;
        res.digits[i] %= BASE;
    }
    while (res.size > 1 && res.digits[res.size] == 0) res.pop();
    return res;
}
BigInteger BigInteger::operator-(const BigInteger& x) const {
    if (!x.flag) return *this + x.abs();
    if (!flag) return -(abs() + x);
    BigInteger res;
    if (*this < x) res.flag = false;
    long long carry = 0;
    int n = std::max(size, x.size);
    res.reserve(n);
    for (int i = 1; i <= n; i++) {
        long long d1 = i <= size ? digits[i] : 0, d2 = i <= x.size ? x.digits[i] : 0;
        if (res.flag) res.push(d1 - d2 - carry);
        else res.push(d2 - d1 - carry);
        if (res.digits[i] < 0) res.digits[i] += BASE, carry = 1;
        else carry = 0;
    }
    while (res.size > 1 && res.digits[res.size] == 0) res.pop();
    return res;
}
namespace FastFourierTransformation {
    constexpr long long FFT_BASE = 1e4;
    constexpr double PI2 = 6.283185307179586231995927;
    constexpr double PI6 = 18.84955592153875869598778;
    constexpr int RECALC_WIDTH = 10;
    constexpr int RECALC_BASE = (1 << RECALC_WIDTH) - 1;
    struct complex {
        double real, imag;
        complex(double x = 0.0, double y = 0.0) : real(x), imag(y) {}
        complex operator+(const complex& other) const { return complex(real + other.real, imag + other.imag); }
        complex operator-(const complex& other) const { return complex(real - other.real, imag - other.imag); }
        complex operator*(const complex& other) const { return complex(real * other.real - imag * other.imag, real * other.imag + other.real * imag); }
        complex& operator+=(const complex& other) { return real += other.real, imag += other.imag, *this; }
        complex& operator-=(const complex& other) { return real -= other.real, imag -= other.imag, *this; }
        complex& operator*=(const complex& other) { return *this = *this * other; }
    };
    inline thread_local complex* arr = nullptr;
    inline void init(int n) {
        delete[] arr;
        arr = new (std::nothrow) complex[n + 1];
        if (!arr) throw std::bad_alloc();
    }
    template <const int n> inline void fft(complex* a) {
        const int n2 = n >> 1, n4 = n >> 2;
        complex w(1.0, 0.0), w3(1.0, 0.0);
        const complex wn(std::cos(PI2 / n), std::sin(PI2 / n)), wn3(std::cos(PI6 / n), std::sin(PI6 / n));
        for (int i = 0; i < n4; i++, w *= wn, w3 *= wn3) {
            if (!(i & RECALC_BASE)) w = complex(std::cos(PI2 * i / n), std::sin(PI2 * i / n)), w3 = w * w * w;
            complex x = a[i] - a[i + n2], y = a[i + n4] - a[i + n2 + n4];
            y = complex(y.imag, -y.real);
            a[i] += a[i + n2], a[i + n4] += a[i + n2 + n4];
            a[i + n2] = (x - y) * w, a[i + n2 + n4] = (x + y) * w3;
        }
        fft<n2>(a), fft<n4>(a + n2), fft<n4>(a + n2 + n4);
    }
    template <> inline void fft<0>(complex*) {}
    template <> inline void fft<1>(complex*) {}
    template <> inline void fft<2>(complex* a) {
        complex x = a[0], y = a[1];
        a[0] += y, a[1] = x - y;
    }
    template <> inline void fft<4>(complex* a) {
        complex a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
        complex x = a0 - a2, y = a1 - a3;
        y = complex(y.imag, -y.real);
        a[0] += a2, a[1] += a3, a[2] = x - y, a[3] = x + y;
        fft<2>(a);
    }
    template <const int n> inline void ifft(complex* a) {
        const int n2 = n >> 1, n4 = n >> 2;
        ifft<n2>(a), ifft<n4>(a + n2), ifft<n4>(a + n2 + n4);
        complex w(1.0, 0.0), w3(1.0, 0.0);
        const complex wn(std::cos(PI2 / n), -std::sin(PI2 / n)), wn3(std::cos(PI6 / n), -std::sin(PI6 / n));
        for (int i = 0; i < n4; i++, w *= wn, w3 *= wn3) {
            if (!(i & RECALC_BASE)) w = complex(std::cos(PI2 * i / n), -std::sin(PI2 * i / n)), w3 = w * w * w;
            complex p = w * a[i + n2], q = w3 * a[i + n2 + n4];
            complex x = a[i], y = p + q, x1 = a[i + n4], y1 = p - q;
            y1 = complex(y1.imag, -y1.real);
            a[i] += y, a[i + n4] += y1, a[i + n2] = x - y, a[i + n2 + n4] = x1 - y1;
        }
    }
    template <> inline void ifft<0>(complex*) {}
    template <> inline void ifft<1>(complex*) {}
    template <> inline void ifft<2>(complex* a) {
        complex x = a[0], y = a[1];
        a[0] += y, a[1] = x - y;
    }
    template <> inline void ifft<4>(complex* a) {
        ifft<2>(a);
        complex p = a[2], q = a[3];
        complex x = a[0], y = p + q, x1 = a[1], y1 = p - q;
        y1 = complex(y1.imag, -y1.real);
        a[0] += y, a[1] += y1, a[2] = x - y, a[3] = x1 - y1;
    }
    inline void dft(complex* a, int n) {
        if (n <= 1) return;
        switch (n) {
            case 1 << 2: fft<1 << 2>(a); break;
            case 1 << 3: fft<1 << 3>(a); break;
            case 1 << 4: fft<1 << 4>(a); break;
            case 1 << 5: fft<1 << 5>(a); break;
            case 1 << 6: fft<1 << 6>(a); break;
            case 1 << 7: fft<1 << 7>(a); break;
            case 1 << 8: fft<1 << 8>(a); break;
            case 1 << 9: fft<1 << 9>(a); break;
            case 1 << 10: fft<1 << 10>(a); break;
            case 1 << 11: fft<1 << 11>(a); break;
            case 1 << 12: fft<1 << 12>(a); break;
            case 1 << 13: fft<1 << 13>(a); break;
            case 1 << 14: fft<1 << 14>(a); break;
            case 1 << 15: fft<1 << 15>(a); break;
            case 1 << 16: fft<1 << 16>(a); break;
            case 1 << 17: fft<1 << 17>(a); break;
            case 1 << 18: fft<1 << 18>(a); break;
            case 1 << 19: fft<1 << 19>(a); break;
            case 1 << 20: fft<1 << 20>(a); break;
            case 1 << 21: fft<1 << 21>(a); break;
            case 1 << 22: fft<1 << 22>(a); break;
            case 1 << 23: fft<1 << 23>(a); break;
            case 1 << 24: fft<1 << 24>(a); break;
            case 1 << 25: fft<1 << 25>(a); break;
            case 1 << 26: fft<1 << 26>(a); break;
            case 1 << 27: fft<1 << 27>(a); break;
            case 1 << 28: fft<1 << 28>(a); break;
            throw FFTLimitExceededError();
        }
    }
    inline void idft(complex* a, int n) {
        if (n <= 1) return;
        switch (n) {
            case 1 << 2: ifft<1 << 2>(a); break;
            case 1 << 3: ifft<1 << 3>(a); break;
            case 1 << 4: ifft<1 << 4>(a); break;
            case 1 << 5: ifft<1 << 5>(a); break;
            case 1 << 6: ifft<1 << 6>(a); break;
            case 1 << 7: ifft<1 << 7>(a); break;
            case 1 << 8: ifft<1 << 8>(a); break;
            case 1 << 9: ifft<1 << 9>(a); break;
            case 1 << 10: ifft<1 << 10>(a); break;
            case 1 << 11: ifft<1 << 11>(a); break;
            case 1 << 12: ifft<1 << 12>(a); break;
            case 1 << 13: ifft<1 << 13>(a); break;
            case 1 << 14: ifft<1 << 14>(a); break;
            case 1 << 15: ifft<1 << 15>(a); break;
            case 1 << 16: ifft<1 << 16>(a); break;
            case 1 << 17: ifft<1 << 17>(a); break;
            case 1 << 18: ifft<1 << 18>(a); break;
            case 1 << 19: ifft<1 << 19>(a); break;
            case 1 << 20: ifft<1 << 20>(a); break;
            case 1 << 21: ifft<1 << 21>(a); break;
            case 1 << 22: ifft<1 << 22>(a); break;
            case 1 << 23: ifft<1 << 23>(a); break;
            case 1 << 24: ifft<1 << 24>(a); break;
            case 1 << 25: ifft<1 << 25>(a); break;
            case 1 << 26: ifft<1 << 26>(a); break;
            case 1 << 27: ifft<1 << 27>(a); break;
            case 1 << 28: ifft<1 << 28>(a); break;
            throw FFTLimitExceededError();
        }
    }
}
BigInteger BigInteger::fft_mul(const BigInteger& a, const BigInteger& b) {
    int least = (a.size + b.size) << 1, lim = 1 << std::__lg(least);
    if (lim < least) lim <<= 1;
    FastFourierTransformation::init(lim);
    for (int i = 0; i < a.size; i++) {
        FastFourierTransformation::arr[i << 1].real = a.digits[i + 1] % FastFourierTransformation::FFT_BASE;
        FastFourierTransformation::arr[i << 1 | 1].real = a.digits[i + 1] / FastFourierTransformation::FFT_BASE % FastFourierTransformation::FFT_BASE;
    }
    for (int i = 0; i < b.size; i++) {
        FastFourierTransformation::arr[i << 1].imag = b.digits[i + 1] % FastFourierTransformation::FFT_BASE;
        FastFourierTransformation::arr[i << 1 | 1].imag = b.digits[i + 1] / FastFourierTransformation::FFT_BASE % FastFourierTransformation::FFT_BASE;
    }
    FastFourierTransformation::dft(FastFourierTransformation::arr, lim);
    for (int i = 0; i < lim; i++) FastFourierTransformation::arr[i] *= FastFourierTransformation::arr[i];
    FastFourierTransformation::idft(FastFourierTransformation::arr, lim);
    BigInteger res;
    res.resize(a.size + b.size + 1);
    long long carry = 0;
    double inv = 0.5 / lim;
    for (int i = 0; i <= a.size + b.size; i++) {
        carry += (long long)(FastFourierTransformation::arr[i << 1].imag * inv + 0.5);
        carry += (long long)(FastFourierTransformation::arr[i << 1 | 1].imag * inv + 0.5) * FastFourierTransformation::FFT_BASE;
        res.digits[i + 1] += carry % BASE, carry /= BASE;
    }
    while (res.size > 1 && res.digits[res.size] == 0) res.pop();
    return res;
}
BigInteger BigInteger::operator*(const BigInteger& x) const {
    BigInteger zero = 0;
    if (*this == zero || x == zero) return zero;
    int n = size, m = x.size;
    long long lim = 1LL * n * m;
    if (lim >= FFT_LIMIT) {
        BigInteger res = fft_mul(*this, x);
        res.flag = !(flag ^ x.flag);
        return res;
    }
    BigInteger res;
    res.flag = !(flag ^ x.flag);
    res.resize(n + m + 2);
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            res.digits[i + j - 1] += digits[i] * x.digits[j];
            res.digits[i + j] += res.digits[i + j - 1] / BASE;
            res.digits[i + j - 1] %= BASE;
        }
    }
    for (int i = 1; i <= n + m + 1; i++) {
        res.digits[i + 1] += res.digits[i] / BASE;
        res.digits[i] %= BASE;
    }
    while (res.size > 1 && res.digits[res.size] == 0) res.pop();
    return res;
}
BigInteger& BigInteger::operator*=(long long x) {
    if (x == 0 || *this == 0) return *this = 0ll;
    if (x < 0) flag = !flag, x = -x;
    long long carry = 0;
    for (int i = 1; i <= size || carry; i++) {
        if (i > size) push(0);
        long long cur = digits[i] * x + carry;
        carry = cur / BigInteger::BASE;
        digits[i] = cur % BigInteger::BASE;
    }
    while (size > 1 && digits[size] == 0) pop();
    return *this;
}
BigInteger BigInteger::operator*(const long long& x) const {
    BigInteger t = *this;
    return t *= x;
}
BigInteger BigInteger::div2() const {
    BigInteger res = *this;
    for (int i = size; i >= 1; i--) {
        if ((res.digits[i] & 1) && (i > 1)) res.digits[i - 1] += BASE;
        res.digits[i] >>= 1;
    }
    while (res.size > 1 && res.digits[res.size] == 0) res.pop();
    return res;
}
BigInteger BigInteger::operator/(const long long& x) const {
    if (x == 0) throw ZeroDivisionError();
    if (*this == 0) return 0;
    if (x == 2) return div2();
    if (x == -2) return -div2();
    BigInteger res;
    res.flag = !(flag ^ (x >= 0));
    long long cur = 0, div = std::abs(x);
    res.resize(size);
    for (int i = size; i >= 1; i--) {
        cur = cur * BASE + digits[i];
        res.digits[i] = res.flag ? (cur / div) : (-cur / -div);
        cur %= div;
    }
    while (res.size > 1 && res.digits[res.size] == 0) res.pop();
    return res;
}
inline BigInteger BigInteger::_move_r(int d) const {
    if (*this == 0 || d >= size) return 0;
    if (d == 0) return *this;
    BigInteger res;
    res.reserve(size - d + 1);
    for (int i = d + 1; i <= size; i++) res.push(digits[i]);
    return res;
}
inline BigInteger BigInteger::_move_l(int d) const {
    if (*this == 0) return 0;
    if (d == 0) return *this;
    BigInteger res;
    res.reserve(size + d + 1);
    for (int i = 1; i <= d; i++) res.push(0);
    for (int i = 1; i <= size; i++) res.push(digits[i]);
    return res;
}
BigInteger BigInteger::newton_inv(int n) const {
    if (*this == 0) throw ZeroDivisionError();
    if (std::min(size, n - size) <= NEWTON_DIV_MIN_LEVEL) {
        BigInteger a;
        a.resize(n + 1);
        std::memset(a.digits, 0, sizeof(long long) * a.size);
        a.digits[n + 1] = 1;
        return a.divmod(*this, true).first;
    }
    int k = (n - size + 2) >> 1, k2 = k > size ? 0 : size - k;
    BigInteger x = _move_r(k2);
    int n2 = k + x.size;
    BigInteger y = x.newton_inv(n2), a = y + y, b = (*this) * y * y;
    return a._move_l(n - n2 - k2) - b._move_r(2 * (n2 + k2) - n) - 1;
}
std::pair<BigInteger, BigInteger> BigInteger::newton_div(const BigInteger& x) const {
    int k = size - x.size + 2, k2 = k > x.size ? 0 : x.size - k;
    BigInteger x2 = x._move_r(k2);
    if (k2 != 0) x2 += 1;
    int n2 = k + x2.size;
    BigInteger u = (*this) * x2.newton_inv(n2);
    BigInteger q = u._move_r(n2 + k2), r = (*this) - q * x;
    while (r >= x) q += 1, r -= x;
    return {q, r};
}
std::pair<BigInteger, BigInteger> BigInteger::divmod(const BigInteger& x, bool dis_newton) const {
    static const int base = BigInteger::BASE;
    BigInteger a = abs(), b = x.abs();
    if (b == 0) throw ZeroDivisionError();
    if (a < b) return {0, flag ? a : -a};
    if (!dis_newton && size > NEWTON_DIV_LIMIT) return newton_div(x);
    int t = base / (x.digits[x.size] + 1);
    a *= t, b *= t;
    int n = a.size, m = b.size;
    BigInteger q = 0, r = 0;
    q.resize(n);
    for (int i = n; i >= 1; i--) {
        r *= base, r += a.digits[i];
        long long d1 = m < r.size ? r.digits[m + 1] : 0, d2 = m - 1 < r.size ? r.digits[m] : 0;
        int d = (d1 * base + d2) / b.digits[m];
        r -= b * d;
        while (!r.flag) r += b, d--;
        q.digits[i] = d;
    }
    q.flag = !(flag ^ x.flag), r.flag = flag;
    while (q.size > 1 && q.digits[q.size] == 0) q.pop();
    return {q, r / t};
}
BigInteger BigInteger::operator/(const BigInteger& x) const { return divmod(x).first; }
BigInteger BigInteger::operator%(const long long& x) const {
    if (x == 2 || x == 4 || x == 5) return digits[1] % x;
    return *this - (*this / x * x);
}
BigInteger BigInteger::operator%(const BigInteger& x) const { return divmod(x).second; }
BigInteger BigInteger::pow(const long long& x) const {
    BigInteger res = 1, a = *this;
    for (long long t = x; t != 0; t >>= 1) {
        if (t & 1) res = res * a;
        a = a * a;
    }
    return res;
}
template <typename T>
BigInteger BigInteger::pow(const long long& x, const T& p) const {
    BigInteger res = 1, a = *this % p;
    for (long long t = x; t != 0; t >>= 1) {
        if (t & 1) res = res * a % p;
        a = a * a % p;
    }
    return res;
}
BigInteger BigInteger::root(const long long& m) const {
    if (!flag && m % 2 == 0) throw NegativeRadicandError();
    if (m == 1 || *this == 0) return *this;
    if (m == 2) return sqrt();
    static constexpr long long base = BigInteger::BASE;
    if (size <= m) {
        long long l = 0, r = base - 1;
        while (l < r) {
            long long mid = (l + r + 1) >> 1;
            if (BigInteger(mid).pow(m) <= *this) l = mid;
            else r = mid - 1;
        }
        return l;
    }
    if (size <= m * 2) {
        BigInteger res; res.resize(2);
        long long l = 0, r = base - 1;
        while (l < r) {
            long long mid = (l + r + 1) >> 1;
            res.digits[2] = mid;
            if (res.pow(m) <= *this) l = mid;
            else r = mid - 1;
        }
        res.digits[2] = l, l = 0, r = base - 1;
        while (l < r) {
            long long mid = (l + r + 1) >> 1;
            res.digits[1] = mid;
            if (res.pow(m) <= *this) l = mid;
            else r = mid - 1;
        }
        return res.digits[1] = l, res;
    }
    int n = size, t = n / m / 2;
    BigInteger s = (_move_r(t * m).root(m) + 1)._move_l(t);
    BigInteger res = (s * (m - 1) + *this / s.pow(m - 1)) / m;
    long long l = std::max<long long>(res.digits[1] - 100, 0), r = std::min(res.digits[1] + 100, base - 1);
    while (l < r) {
        long long mid = (l + r + 1) >> 1;
        res.digits[1] = mid;
        if (res.pow(m) <= *this) l = mid;
        else r = mid - 1;
    }
    return res.digits[1] = l, res;
}
BigInteger BigInteger::newton_invsqrt() const {
    const BigInteger& x = *this;
    static constexpr long long base = BigInteger::BASE;
    if (x.size <= 2) {
        if (x.size == 0) return 0;
        uint64_t t = 0;
        if (x.size == 1) t = x.digits[1];
        else t = x.digits[1] + x.digits[2] * BASE;
        return uint64_t(base * base / std::sqrt(t));
    }
    if (x.size <= NEWTON_SQRT_MIN_LEVEL) {
        BigInteger b1 = BigInteger(base).pow(x.size / 2);
        return b1 / x.sqrt_normal();
    }
    int n2 = x.size % 2 == 0 ? x.size : x.size + 1, k2 = (n2 + 2) / 4 * 2;
    BigInteger x2k(x.digits + 1 + (n2 - k2), k2 + x.size - n2);
    BigInteger s = x2k.newton_invsqrt()._move_l((n2 - k2) / 2);
    BigInteger x2 = (s * 3).div2() - (s * s * s * x).div2()._move_r(2 * n2);
    BigInteger rx = BigInteger(1)._move_l(2 * n2) - x * x2 * x2, delta = 1;
    if (!rx.flag) {
        for (; !rx.flag; delta *= 2) {
            BigInteger t = (x2 * 2 - delta + delta * delta) * x;
            x2 -= delta, rx += t;
        }
    }
    else {
        while (true) {
            BigInteger t = (x2 * 2 + delta) * delta * x;
            if (t > rx) break;
            x2 += delta, rx -= t;
            delta *= 2;
        }
    }
    for (; delta > 0; delta = delta.div2()) {
        BigInteger t = (x2 * 2 + delta) * delta * x;
        if (t <= rx) x2 += delta, rx -= t;
    }
    return x2;
}
BigInteger BigInteger::sqrt_normal() const {
    static constexpr long long base = BigInteger::BASE;
    BigInteger n = *this, x0 = BigInteger(base)._move_l((n.size + 2) >> 1);
    BigInteger x = (x0 + n / x0).div2();
    while (x < x0) std::swap(x, x0), x = (x0 + n / x0).div2();
    return x0;
}
BigInteger BigInteger::sqrt() const {
    if (!flag) throw NegativeRadicandError();
    if (*this <= 1) return *this;
    if (size <= NEWTON_SQRT_LIMIT) return sqrt_normal();
    const BigInteger& x = *this;
    int n2 = x.size % 2 == 0 ? x.size : x.size + 1;
    BigInteger res = (x * x.newton_invsqrt())._move_r(n2);
    BigInteger r = x - res * res, delta = 1;
    while (true) {
        BigInteger dr = (res * 2 + delta) * delta;
        if (dr > r) break;
        r -= dr, res += delta, delta *= 2;
    }
    for (; delta > 0; delta = delta.div2()) {
        BigInteger dr = (res * 2 + delta) * delta;
        if (dr <= r) r -= dr, res += delta;
    }
    return res;
}
BigInteger BigInteger::gcd(const BigInteger& x) const {
	BigInteger a = *this, b = x;
	if (a < b) std::swap(a, b);
	if (b == 0) return a;
	int t = 0;
	while (a % 2 == 0 && b % 2 == 0) a = a.div2(), b = b.div2(), t++;
	while (b > 0) {
		if (a % 2 == 0) a = a.div2();
		else if (b % 2 == 0) b = b.div2();
		else a -= b;
		if (a < b) std::swap(a, b);
	}
	while (t--) a += a;
	return a;
}
BigInteger BigInteger::lcm(const BigInteger& x) const {
    return *this / gcd(x) * x;
}
BigInteger& BigInteger::operator+= (const BigInteger& x) {return *this = *this + x;}
BigInteger& BigInteger::operator-= (const BigInteger& x) {return *this = *this - x;}
BigInteger& BigInteger::operator*= (const BigInteger& x) {return *this = *this * x;}
BigInteger& BigInteger::operator/= (const long long& x) {return *this = *this / x;}
BigInteger& BigInteger::operator/= (const BigInteger& x) {return *this = *this / x;}
BigInteger& BigInteger::operator%= (const long long& x) {return *this = *this % x;}
BigInteger& BigInteger::operator%= (const BigInteger& x) {return *this = *this % x;}
BigInteger BigInteger::operator<< (const long long& x) {return *this * BigInteger(2).pow(x);}
BigInteger BigInteger::operator>> (const long long& x) {return *this / BigInteger(2).pow(x);}
BigInteger& BigInteger::operator<<= (const long long& x) {return *this = *this << x;}
BigInteger& BigInteger::operator>>= (const long long& x) {return *this = *this >> x;}
BigInteger BigInteger::operator& (const BigInteger& x) {
 std::vector<bool> a = to_binary(), b = x.to_binary();
	int n = a.size(), m = b.size(), lim = std::max(n, m);
	std::vector<bool> res(lim), temp(lim);
	if (m == lim) {
		a.resize(lim);
		for (int i = 0; i < n; ++i) temp[lim - n + i] = a[i];
		a = temp;
	}
	else {
		b.resize(lim);
		for (int i = 0; i < m; ++i) temp[lim - m + i] = b[i];
		b = temp;
	}
	for (int i = 0; i < lim; ++i) res[i] = a[i] & b[i];
	return res;
}
BigInteger BigInteger::operator| (const BigInteger& x) {
	std::vector<bool> a = to_binary(), b = x.to_binary();
	int n = a.size(), m = b.size(), lim = std::max(n, m);
	std::vector<bool> res(lim), temp(lim);
	if (m == lim) {
		a.resize(lim);
		for (int i = 0; i < n; ++i) temp[lim - n + i] = a[i];
		a = temp;
	}
	else {
		b.resize(lim);
		for (int i = 0; i < m; ++i) temp[lim - m + i] = b[i];
		b = temp;
	}
	for (int i = 0; i < lim; ++i) res[i] = a[i] | b[i];
	return res;
}
BigInteger BigInteger::operator^ (const BigInteger& x) {
	std::vector<bool> a = to_binary(), b = x.to_binary();
	int n = a.size(), m = b.size(), lim = std::max(n, m);
	std::vector<bool> res(lim), temp(lim);
	if (m == lim) {
		a.resize(lim);
		for (int i = 0; i < n; ++i) temp[lim - n + i] = a[i];
		a = temp;
	}
	else {
		b.resize(lim);
		for (int i = 0; i < m; ++i) temp[lim - m + i] = b[i];
		b = temp;
	}
	for (int i = 0; i < lim; ++i) res[i] = a[i] ^ b[i];
	return res;
}
BigInteger& BigInteger::operator&= (const BigInteger& x) {return *this = *this & x;}
BigInteger& BigInteger::operator|= (const BigInteger& x) {return *this = *this | x;}
BigInteger& BigInteger::operator^= (const BigInteger& x) {return *this = *this ^ x;}
BigInteger& BigInteger::operator++ () {return *this += 1;}
BigInteger BigInteger::operator++ (int) {BigInteger t = *this; return *this += 1, t;}
BigInteger& BigInteger::operator-- () {return *this -= 1;}
BigInteger BigInteger::operator-- (int) {BigInteger t = *this; return *this -= 1, t;}
BigInteger BigInteger::random(const int& num_digits){
    std::random_device rd;
    std::mt19937 gen(rd());
    BigInteger res;
    res.reserve(num_digits / 8 + 8);
	res.size = 0;
    std::uniform_int_distribution<long long> last_digits(1, std::pow(10, (num_digits - 1) % 8 + 1) - 1);
    if (num_digits % 8) res.push(last_digits(gen));
    std::uniform_int_distribution<long long> other_digits(0, 99999999);
    for (int i = 1; i <= num_digits / 8; ++i) res.push(other_digits(gen));
    return res;
}
