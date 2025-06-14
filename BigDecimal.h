#pragma once
#define BIGDECIMAL_H
#include "BigInteger.h"
class NaNError : public std::runtime_error {
public:
	NaNError() : std::runtime_error("Not a number") {}
};
enum __round_types {
	ROUND_CEIL, ROUND_FLOOR, ROUND_05, ROUND_04
};
int __global_precision = -65;
int __global_round_type = ROUND_04;
class BigDecimal {
private:
	BigInteger val;
	int e;
	BigDecimal(const BigInteger& v, int _e) : val(v), e(_e) {reset();}
	std::string to_string_base(std::function<std::string(std::string&)>) const;
	std::string round_045(int = INT_MAX, bool = true) const;
	inline int compare(const BigDecimal&) const;
public:
	BigDecimal() : val(0), e(0) {}
	BigDecimal(const BigInteger& x) {*this = x;}
	BigDecimal(const BigDecimal& x) {*this = x;}
	BigDecimal(const long double& x) {*this = x;}
	BigDecimal(const long double& x, int precision) {
		std::stringstream ss;
		ss << std::setprecision(precision) << x;
		*this = ss.str();
	}
	BigDecimal(const std::string& s) {*this = s;}
	BigDecimal& operator= (const BigInteger&);
	BigDecimal& operator= (const BigDecimal&);
	BigDecimal& operator= (const long double&);
	BigDecimal& operator= (const std::string&);
	static BigInteger move_l(const BigInteger& num, int x) {
		std::string s = num.to_string();
		for (int i = 1; i <= x; i++) s += '0';
		return BigInteger(s);
	}
	static BigInteger move_r(const BigInteger& num, int x) {
		std::string s = num.to_string();
		for (int i = 1; i <= x && !s.empty(); i++) s.pop_back();
		return BigInteger(s);
	}
	int len() const {return val.to_string().size();}
	int precision() const {return e;}
	void reset() {
		const int p = __global_precision;
		if (e < p) val = move_r(val, p - e);
		else if (e > p) val = move_l(val, e - p);
		e = p;
	}
	std::string to_string(int = INT_MAX) const;
	std::string floor(int = INT_MAX) const;
	std::string ceil(int = INT_MAX) const;
	std::string round_04(int = INT_MAX) const;
	std::string round_05(int = INT_MAX) const;
	BigInteger to_bigint() const;
	friend std::ostream& operator<< (std::ostream& out, const BigDecimal& x) {
		return out << x.to_string();
	}
	friend std::istream& operator>> (std::istream& in, BigDecimal& x) {
		std::string s;
		in >> s;
		x = s;
		return in;
	}
	bool operator== (const BigDecimal&) const;
	auto operator<=> (const BigDecimal&) const;
	BigDecimal operator+ (const BigDecimal&) const;
	BigDecimal operator- (const BigDecimal&) const;
	BigDecimal operator* (const BigDecimal&) const;
	BigDecimal operator/ (const BigDecimal&) const;
	BigDecimal operator% (const BigDecimal&) const;
	BigDecimal& operator+= (const BigDecimal&);
	BigDecimal& operator-= (const BigDecimal&);
	BigDecimal& operator*= (const BigDecimal&);
	BigDecimal& operator/= (const BigDecimal&);
	BigDecimal& operator%= (const BigDecimal&);
	bool is_zero() const;
	BigDecimal abs() const;
	BigDecimal operator- () const;
	BigDecimal exp() const;
	BigDecimal ln() const;
	BigDecimal sqrt() const;
	BigDecimal pow(const BigDecimal&) const;
	BigDecimal acos() const;
	BigDecimal asin() const;
	BigDecimal atan() const;
	BigDecimal atan2(const BigDecimal&) const;
	BigDecimal cos() const;
	BigDecimal sin() const;
	BigDecimal tan() const;
};
const BigDecimal b_pi = (std::string) "3.1415926535897932384626433832795028841971693993751058209749445923";
const BigDecimal b_ln2 = (std::string) "0.693147180559945309417232121458176568075500134360255254120680009";
const BigDecimal b_sqrt2 = (std::string) "1.4142135623730950488016887242096980785696718753769480731766797";
const BigDecimal b_one = 1;
BigDecimal& BigDecimal::operator= (const BigInteger& x) {
	e = 0, val = x, reset();
	return *this;
}
BigDecimal& BigDecimal::operator= (const BigDecimal& x) {
	e = x.e, val = x.val, reset();
	return *this;
}
BigDecimal& BigDecimal::operator= (const long double& x) {
	std::stringstream ss;
	ss << std::setprecision(24) << x;
	return *this = ss.str();
}
BigDecimal& BigDecimal::operator= (const std::string& s) {
	e = 0, val = 0, reset();
	if (s.empty()) return *this;
	int n = s.size(), pos = 0, flag = 1;
	while (pos < n && s[pos] == '-') flag = !flag, pos++;
	std::string a, b;
	while (pos < n && std::isdigit(s[pos])) a.push_back(s[pos]), pos++;
	if (s[pos] != '.') {
		val = a, e = 0;
		if (!flag) val = -val;
		return reset(), *this;
	}
	for (pos++; pos < n && std::isdigit(s[pos]); pos++) b.push_back(s[pos]);
	e = -static_cast<long long>(b.size()), val = a + b;
	if (!flag) val = -val;
	return reset(), *this;
}
std::string BigDecimal::to_string_base(std::function<std::string(std::string&)> solve) const {
	bool neg = false;
	std::string res = val.to_string();
	if (val < 0) res = res.substr(1), neg = true;
	auto add_neg = [&](const std::string& s) -> std::string {
		return neg ? ('-' + s) : s;
	};
	if (e == 0) return add_neg(solve(res));
	if (e > 0) {
		for (int i = 1; i <= e; i++) res += '0';
		return add_neg(solve(res));
	}
	int n = res.size(), point = n + e;
	if (point >= 0) {
		res.insert(res.begin() + point, '.');
		if (res[0] == '.') res = '0' + res;
		if ((int) res.size() >= 2 && res[0] == '-' && res[1] == '.') res.insert(res.begin() + 1, '0');
		return add_neg(solve(res));
	}
	std::string pre = "0.";
	for (int i = 1; i <= -point; i++) pre += '0';
	pre += res;
	return add_neg(solve(pre));
}
std::string BigDecimal::to_string(int precision) const {
	switch (__global_round_type) {
		case ROUND_CEIL:
			return ceil(precision);
		case ROUND_FLOOR:
			return floor(precision);
		case ROUND_04:
			return round_04(precision);
		case ROUND_05:
			return round_05(precision);
	}
	return "";
}
std::string BigDecimal::floor(int precision) const {
	int p = precision != INT_MAX ? -precision + 1 : -__global_precision;
	return to_string_base([&](std::string& s) -> std::string {
		int pos = s.find('.');
		if (pos == (int) std::string::npos) s += ".0", pos = s.find('.');
		int n = s.size(), sep = n - pos;
		if (sep == p) return s;
		if (sep > p) {
			for (int i = 1; i <= sep - p; i++) s.pop_back();
			if (s.back() == '.') s.pop_back();
			return s;
		}
		for (int i = 1; i <= p - sep; i++) s += '0';
		return s;
	});
}
std::string BigDecimal::ceil(int precision) const {
	int p = precision != INT_MAX ? -precision + 2 : -__global_precision + 1;
	return to_string_base([&](std::string& s) -> std::string {
		int pos = s.find('.');
		if (pos == (int) std::string::npos) s += ".0", pos = s.find('.');
		int n = s.size(), sep = n - pos;
		if (sep > p) {
			for (int i = 1; i <= sep - p; i++) s.pop_back();
		} else if (sep < p) {
			for (int i = 1; i <= p - sep; i++) s += '0';
		}
		int i = s.size() - 1;
		for (i--; s[i] == '.' || s[i] == '9'; i--) {
			if (s[i] == '.') continue;
			s[i] = '0';
			if (i == 0) {
				s.insert(s.begin(), '0');
				break;
			}
		}
		s[i]++, s.pop_back();
		if (s.back() == '.') s.pop_back();
		return s;
	});
}
std::string BigDecimal::round_045(int precision, bool round04) const {
	int p = precision != INT_MAX ? -precision + 2 : -__global_precision + 1;
	return to_string_base([&](std::string& s) -> std::string {
		int pos = s.find('.');
		if (pos == (int) std::string::npos) s += ".0", pos = s.find('.');
		int n = s.size(), sep = n - pos;
		if (sep > p) {
			for (int i = 1; i <= sep - p; i++) s.pop_back();
		} else if (sep < p) {
			for (int i = 1; i <= p - sep; i++) s += '0';
		}
		int i = s.size() - 1;
		bool flag = round04 ? (s[i] >= '5') : (s[i] > '5');
		if (flag) {
			for (i--; s[i] == '.' || s[i] == '9'; i--) {
				if (s[i] == '.') continue;
				s[i] = '0';
				if (i == 0) {
					s.insert(s.begin(), '0');
					break;
				}
			}
			s[i]++;
		}
		s.pop_back();
		if (s.back() == '.') s.pop_back();
		return s;
	});
}
std::string BigDecimal::round_04(int precision) const {
	return round_045(precision, true);
}
std::string BigDecimal::round_05(int precision) const {
	return round_045(precision, false);
}
BigInteger BigDecimal::to_bigint() const {
	BigInteger res = val;
	if (e < 0) res = move_r(res, -e);
	else if (e > 0) res = move_l(res, e);
	return res;
}
BigDecimal& BigDecimal::operator+= (const BigDecimal& x) {
	return reset(), val += x.val, *this;
}
BigDecimal& BigDecimal::operator-= (const BigDecimal& x) {
	return reset(), val -= x.val, *this;
}
BigDecimal& BigDecimal::operator*= (const BigDecimal& x) {
	return reset(), val *= x.val, e += x.e, reset(), *this;
}
BigDecimal& BigDecimal::operator/= (const BigDecimal& x) {
	if (x.is_zero()) throw NaNError();
	int n = x.len();
	reset(), val = move_l(val, n), e -= n;
	return val /= x.val, e -= x.e, reset(), *this;
}
BigDecimal& BigDecimal::operator%= (const BigDecimal& x) {
	if (x.is_zero()) throw NaNError();
	return *this -= x * (*this / x).to_bigint();
}
BigDecimal BigDecimal::operator+ (const BigDecimal& x) const {
	return BigDecimal(val + x.val, e);
}
BigDecimal BigDecimal::operator- (const BigDecimal& x) const {
	return BigDecimal(val - x.val, e);
}
BigDecimal BigDecimal::operator* (const BigDecimal& x) const {
	return BigDecimal(val * x.val, e + x.e);
}
BigDecimal BigDecimal::operator/ (const BigDecimal& x) const {
	if (x.is_zero()) throw NaNError();
	BigDecimal res = *this;
	res /= x;
	return res;
}
BigDecimal BigDecimal::operator% (const BigDecimal& x) const {
	if (x.is_zero()) throw NaNError();
	return *this - x * (*this / x).to_bigint();
}
bool BigDecimal::is_zero() const {
	return val == 0;
}
BigDecimal BigDecimal::abs() const {
	BigDecimal res = *this;
	res.val = res.val.abs();
	return res;
}
BigDecimal BigDecimal::operator- () const {
	BigDecimal res = *this;
	res.val = -res.val;
	return res;
}
int BigDecimal::compare(const BigDecimal& y) const {
	if (is_zero() && y.is_zero()) return 0;
	if (val < 0 && y.val > 0) return -1;
	if (val > 0 && y.val < 0) return 1;
	int sgn = 1;
	if (val < 0 && y.val < 0) sgn = -1;
	if (e != y.e) {
		if (e < y.e) return -sgn;
		return sgn;
	}
	if (val < y.val) return -sgn;
	if (val > y.val) return sgn;
	return 0;
}
bool BigDecimal::operator== (const BigDecimal& x) const {return compare(x) == 0;}
auto BigDecimal::operator<=> (const BigDecimal& x) const {return compare(x);}
BigDecimal BigDecimal::exp() const {
	__global_precision -= 2;
	BigDecimal a = 1, res = 0, h = *this;
	int k = 0;
	for (BigInteger fac = 1; ; k++, fac *= k) {
		BigDecimal b = a / fac;
		if (b.is_zero()) break;
		res += b, a *= h;
	}
	__global_precision += 2, res.reset();
	return res;
}
BigDecimal BigDecimal::ln() const {
	if (*this < 1) return -((b_one / *this).ln());
	__global_precision -= 2;
	BigDecimal coef = 0, x = *this, ln2 = b_ln2, sqrt2 = b_sqrt2;
	while (x > sqrt2) {
		x.val = x.val.div2();
		coef += ln2;
	}
	x = (x - 1) / (x + 1);
	BigDecimal cnt = x * 2, res = cnt + coef;
	x *= x;
	for (int i = 3; ; i += 2) {
		cnt *= x;
		BigDecimal b = cnt / i;
		if (b.is_zero()) break;
		res += b;
	}
	__global_precision += 2, res.reset();
	return res;
}
BigDecimal BigDecimal::pow(const BigDecimal& x) const {
	return (ln() * x).exp();
}
BigDecimal BigDecimal::sqrt() const {
	BigDecimal x = *this;
	int n = x.len();
	x.val = move_l(x.val, n), x.e -= n;
	if (x.e & 1) x.val = move_r(x.val, 1), x.e++;
	x.val = x.val.sqrt(), x.e /= 2, x.reset();
	return x;
}
BigDecimal BigDecimal::asin() const {
	BigDecimal x = *this;
	return (x / (b_one - x * x).sqrt()).atan();
}
BigDecimal BigDecimal::acos() const {
	return b_pi / 2 - asin();
}
BigDecimal BigDecimal::atan() const {
	if (val < 0) return -(abs().atan());
	BigDecimal x = *this;
	if (x > 1) return b_pi / 2 - (b_one / x).atan();
	if (x >= 0.1) return (x / ((x * x + 1).sqrt() + 1)).atan() * 2;
	__global_precision -= 2, x = *this;
	BigDecimal sgn = x, res = x;
	x = -x * x;
	for (int i = 3; ; i += 2) {
		sgn *= x;
		BigDecimal b = sgn / i;
		if (b.is_zero()) break;
		res += b;
	}
	__global_precision += 2, res.reset();
	return res;
}
BigDecimal BigDecimal::atan2(const BigDecimal& x) const {
	BigDecimal y = *this;
	if (x > 0) return (y / x).atan();
	if (y >= 0 && x < 0) return (y / x).atan() + b_pi;
	if (y < 0 && x < 0) return (y / x).atan() - b_pi;
	if (y > 0) return b_pi / 2;
	if (y < 0) return -b_pi / 2;
	throw NaNError();
}
BigDecimal BigDecimal::sin() const {
	if (*this < 0) return -(abs().sin());
	__global_precision -= 2;
	BigDecimal x = *this, y = b_pi * 2;
	x = (x % y + y) % y;
	BigDecimal val = x;
	BigDecimal res = x, cur = x;
	x *= x;
	for (int i = 3; ; i += 2) {
		cur = cur * x / (1LL * i * (i - 1));
		if (cur.is_zero()) break;
		cur = -cur, res += cur;
	}
	__global_precision += 2, res.reset();
	return res;
}
BigDecimal BigDecimal::cos() const {
	if (*this < 0) return abs().cos();
	__global_precision -= 2;
	BigDecimal x = *this, y = b_pi * 2;
	x = (x % y + y) % y;
	BigDecimal val = x;
	BigDecimal res = 1, cur = 1;
	x *= x;
	for (int i = 2; ; i += 2) {
		cur = cur * x / (1LL * i * (i - 1));
		if (cur.is_zero()) break;
		cur = -cur, res += cur;
	}
	__global_precision += 2, res.reset();
	return res;
}
BigDecimal BigDecimal::tan() const {
	return sin() / cos();
}
BigDecimal b_exp(const BigDecimal& x) {return x.exp();}
BigDecimal b_ln(const BigDecimal& x) {return x.ln();}
BigDecimal b_sqrt(const BigDecimal& x) {return x.sqrt();}
BigDecimal b_pow(const BigDecimal& x, const BigDecimal& y) {return x.pow(y);}
BigDecimal b_abs(const BigDecimal& x) {return x.abs();}
BigDecimal b_sin(const BigDecimal& x) {return x.sin();}
BigDecimal b_cos(const BigDecimal& x) {return x.cos();}
BigDecimal b_tan(const BigDecimal& x) {return x.tan();}
BigDecimal b_asin(const BigDecimal& x) {return x.asin();}
BigDecimal b_acos(const BigDecimal& x) {return x.acos();}
BigDecimal b_atan(const BigDecimal& x) {return x.atan();}
BigDecimal b_atan2(const BigDecimal& y, const BigDecimal& x) {return y.atan2(x);}
