#pragma once
#define FASTIO_H
#include <cmath>
#include <cstring>
#include <string>
#include <cstdint>
namespace FastIO {
	static constexpr size_t DEFAULT_BUFFER_SIZE = 1 << 16;
	static constexpr size_t MAX_NUMERIC_BUFFER_SIZE = 64;
	enum class IOError {
		None = 0,
		EndOfFile,
		InvalidInput,
		OutOfRange,
		IOError
	};
	class istream {
	private:
		char* buffer;
		char* current;
		char* end;
		const size_t bufferSize;
		IOError lastError = IOError::None;
		[[nodiscard]] constexpr inline __attribute__((always_inline)) bool fillBuffer() noexcept {
			size_t bytesRead = fread(buffer, 1, bufferSize, stdin);
			current = buffer;
			end = buffer + bytesRead;
			return bytesRead > 0;
		}
		[[nodiscard]] constexpr inline __attribute__((always_inline)) char getChar() noexcept {
			if (current == end) {
				if (!fillBuffer()) {
					lastError = feof(stdin) ? IOError::EndOfFile : IOError::IOError;
					return EOF;
				}
			}
			return *current++;
		}
		[[nodiscard]] constexpr inline __attribute__((always_inline)) char peekChar() noexcept {
			if (current == end && !fillBuffer()) {
				return EOF;
			}
			return *current;
		}
		template <typename T>
		constexpr inline __attribute__((always_inline)) void skipWhitespace() noexcept {
			char ch;
			while ((ch = peekChar()) != EOF && (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r' || ch == '\v' || ch == '\f')) {
				++current;
			}
		}
	public:
		explicit istream(size_t bufSize = DEFAULT_BUFFER_SIZE) : buffer(new char[bufSize]), current(buffer), end(buffer), bufferSize(bufSize) {}
		~istream() noexcept {
			delete[] buffer;
		}
		[[nodiscard]] IOError error() const noexcept {
			return lastError;
		}
		void clear() noexcept {
			lastError = IOError::None;
		}
		template <typename T, typename std::enable_if_t<std::is_unsigned_v<T>, bool> = true>
		inline __attribute__((always_inline)) istream & operator>>(T& x) noexcept {
			clear();
			skipWhitespace<T>();
			x = 0;
			char ch = getChar();
			if (ch == EOF) {
				lastError = IOError::EndOfFile;
				return *this;
			}
			while (ch == '0') {
				ch = getChar();
				if (ch == EOF) {
					lastError = IOError::EndOfFile;
					x = 0;
					return *this;
				}
			}
			if (ch < '0' || ch > '9') {
				lastError = IOError::InvalidInput;
				--current;
				return *this;
			}
			T maxLimit = std::numeric_limits<T>::max() / 10;
			T maxDigit = std::numeric_limits<T>::max() % 10;
			try {
				x = 0;
				do {
					if (x > maxLimit || (x == maxLimit && (ch - '0') > maxDigit)) {
						lastError = IOError::OutOfRange;
						while ((ch = getChar()) >= '0' && ch <= '9');
						return *this;
					}
					x = x * 10 + (ch - '0');
					ch = getChar();
				} while (ch >= '0' && ch <= '9');
				if (ch != EOF) {
					--current;
				}
			} catch (...) {
				lastError = IOError::IOError;
			}
			return *this;
		}
		template < typename T, typename std::enable_if_t < std::is_signed_v<T> && std::is_integral_v<T>, bool > = true >
		inline __attribute__((always_inline)) istream & operator>>(T& x) noexcept {
			clear();
			skipWhitespace<T>();
			x = 0;
			char ch = getChar();
			if (ch == EOF) {
				lastError = IOError::EndOfFile;
				return *this;
			}
			bool negative = false;
			if (ch == '-') {
				negative = true;
				ch = getChar();
			} else if (ch == '+') {
				ch = getChar();
			}
			while (ch == '0') {
				ch = getChar();
				if (ch == EOF) {
					lastError = IOError::EndOfFile;
					x = 0;
					return *this;
				}
			}
			if (ch < '0' || ch > '9') {
				lastError = IOError::InvalidInput;
				--current;
				return *this;
			}
			using UnsignedT = typename std::make_unsigned_t<T>;
			UnsignedT maxLimit = std::numeric_limits<UnsignedT>::max() / 10;
			UnsignedT maxDigit = std::numeric_limits<UnsignedT>::max() % 10;
			try {
				UnsignedT tmp = 0;
				do {
					if (tmp > maxLimit || (tmp == maxLimit && (ch - '0') > maxDigit)) {
						lastError = IOError::OutOfRange;
						while ((ch = getChar()) >= '0' && ch <= '9');
						return *this;
					}
					tmp = tmp * 10 + (ch - '0');
					ch = getChar();
				} while (ch >= '0' && ch <= '9');
				if (ch != EOF) {
					--current;
				}
				if (negative) {
					if (tmp > static_cast<UnsignedT>(std::numeric_limits<T>::max()) + 1) {
						lastError = IOError::OutOfRange;
						x = std::numeric_limits<T>::min();
					} else {
						x = -static_cast<T>(tmp);
					}
				} else {
					if (tmp > static_cast<UnsignedT>(std::numeric_limits<T>::max())) {
						lastError = IOError::OutOfRange;
						x = std::numeric_limits<T>::max();
					} else {
						x = static_cast<T>(tmp);
					}
				}
			} catch (...) {
				lastError = IOError::IOError;
			}
			return *this;
		}
		template <typename T, typename std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
		inline __attribute__((always_inline)) istream & operator>>(T& x) noexcept {
			clear();
			skipWhitespace<T>();
			x = 0;
			char ch = getChar();
			if (ch == EOF) {
				lastError = IOError::EndOfFile;
				x = std::numeric_limits<T>::quiet_NaN();
				return *this;
			}
			bool negative = false;
			if (ch == '-') {
				negative = true;
				ch = getChar();
			} else if (ch == '+') {
				ch = getChar();
			}
			T integerPart = 0;
			bool hasDigits = false;
			while (ch >= '0' && ch <= '9') {
				hasDigits = true;
				integerPart = integerPart * 10 + (ch - '0');
				ch = getChar();
			}
			T fraction = 0;
			T divisor = 1;
			bool hasFraction = false;
			if (ch == '.') {
				ch = getChar();
				while (ch >= '0' && ch <= '9') {
					hasFraction = true;
					fraction = fraction * 10 + (ch - '0');
					divisor *= 10;
					ch = getChar();
				}
			}
			int exponent = 0;
			if (ch == 'e' || ch == 'E') {
				ch = getChar();
				bool expNegative = false;
				if (ch == '-') {
					expNegative = true;
					ch = getChar();
				} else if (ch == '+') {
					ch = getChar();
				}
				while (ch >= '0' && ch <= '9') {
					exponent = exponent * 10 + (ch - '0');
					ch = getChar();
				}
				if (expNegative) {
					exponent = -exponent;
				}
			}
			if (!hasDigits && !hasFraction) {
				lastError = IOError::InvalidInput;
				--current;
				return *this;
			}
			try {
				x = (integerPart + (hasFraction ? fraction / divisor : 0)) * std::pow(10.0, exponent);
				if (negative) x = -x;
			} catch (...) {
				lastError = IOError::IOError;
				x = std::numeric_limits<T>::quiet_NaN();
			}
			return *this;
		}
		template <typename T, typename std::enable_if_t<std::is_same_v<T, std::string>, bool> = true>
		inline __attribute__((always_inline)) istream & operator>>(T& str) noexcept {
			clear();
			skipWhitespace<T>();
			str.clear();
			char ch;
			while ((ch = getChar()) != EOF && !std::isspace(ch)) {
				str.push_back(ch);
			}
			if (ch == EOF && str.empty()) {
				lastError = IOError::EndOfFile;
			}
			return *this;
		}
		inline __attribute__((always_inline)) istream& operator>>(char& ch) noexcept {
			clear();
			skipWhitespace<char>();
			ch = getChar();
			if (ch == EOF) {
				lastError = IOError::EndOfFile;
			}
			return *this;
		}
		inline __attribute__((always_inline)) istream& operator>>(char* str) noexcept {
			clear();
			skipWhitespace<char>();
			char ch;
			size_t i = 0;
			while ((ch = getChar()) != EOF && !(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r' || ch == '\v' || ch == '\f')) {
				if (i < std::numeric_limits<size_t>::max() - 1) {
					str[i++] = ch;
				}
			}
			if (ch == EOF && i == 0) {
				lastError = IOError::EndOfFile;
			}
			str[i] = '\0';
			return *this;
		}
		inline __attribute__((always_inline)) istream& operator>>(bool& b) noexcept {
			clear();
			skipWhitespace<int>();
			b = false;
			char ch = getChar();
			if (ch == EOF) {
				lastError = IOError::EndOfFile;
				return *this;
			}
			if (ch == '-' || ch == '+') ch = getChar();
			while (ch == '0') {
				ch = getChar();
				if (ch == EOF) {
					lastError = IOError::EndOfFile;
					b = 0;
					return *this;
				}
			}
			if (ch < '0' || ch > '9') {
				lastError = IOError::InvalidInput;
				--current;
				return *this;
			}
			do {
				if (ch != '0') {
					b = true;
					return *this;
				}
				ch = getChar();
			} while (ch >= '0' && ch <= '9');
			if (ch != EOF) {
				--current;
			}
			return *this;
		}
	};
	class ostream {
	private:
		char* buffer;
		char* current;
		char* end;
		const size_t bufferSize;
		unsigned int Base = 10;
		unsigned int Precision = 6;
		unsigned int Width = 0;
		char Fill = ' ';
		bool Align = false;
		bool ShowPos = false;
		bool ShowBase = false;
		bool Case = false;
		bool Scientific = false;
		bool Fixed = false;
		bool BoolAlpha = false;
		bool ShowPoint = true;
		constexpr inline __attribute__((always_inline)) void flushBuffer() noexcept {
			if (current != buffer) {
				fwrite(buffer, 1, current - buffer, stdout);
				current = buffer;
			}
		}
		constexpr inline __attribute__((always_inline)) void putChar(char ch) noexcept {
			if (current == end) {
				flushBuffer();
			}
			*current++ = ch;
		}
		template <typename T, typename std::enable_if_t<std::is_integral_v<T>, bool> = true>
		inline __attribute__((always_inline)) void outputInteger(T x) noexcept {
			char numBuf[MAX_NUMERIC_BUFFER_SIZE];
			char* p = numBuf + MAX_NUMERIC_BUFFER_SIZE - 1;
			*p = '\0';
			const char* digits = Case ? "0123456789ABCDEF" : "0123456789abcdef";
			bool isNegative = std::is_signed_v<T>&& x < 0;
			using UnsignedT = typename std::make_unsigned_t<T>;
			UnsignedT ux = isNegative ? -static_cast<UnsignedT>(x) : static_cast<UnsignedT>(x);
			if (ux == 0) {
				*--p = '0';
			} else {
				while (ux > 0) {
					*--p = digits[ux % Base];
					ux /= Base;
				}
			}
			if (ShowBase) {
				switch (Base) {
					case 16:
						*--p = Case ? 'X' : 'x';
						*--p = '0';
						break;
					case 8:
						*--p = Case ? 'O' : 'o';
						*--p = '0';
						break;
					case 2:
						*--p = Case ? 'B' : 'b';
						*--p = '0';
						break;
				}
			}
			if (isNegative) {
				*--p = '-';
			} else if (ShowPos) {
				*--p = '+';
			}
			size_t numLen = (numBuf + MAX_NUMERIC_BUFFER_SIZE - 1) - p;
			if (Width > numLen) {
				size_t padding = Width - numLen;
				if (!Align) {
					for (size_t i = 0; i < padding; ++i) putChar(Fill);
				}
				while (*p) putChar(*p++);
				if (Align) {
					for (size_t i = 0; i < padding; ++i) putChar(Fill);
				}
			} else {
				while (*p) putChar(*p++);
			}
		}
		template <typename T, typename std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
		inline __attribute__((always_inline)) void outputFloating(T x) noexcept {
			if (std::isnan(x)) {
				outputString(Case ? "NAN" : "nan");
				return;
			}
			if (std::isinf(x)) {
				if (x < 0) putChar('-');
				else if (ShowPos) putChar('+');
				outputString(Case ? "INF" : "inf");
				return;
			}
			if (x < 0) {
				putChar('-');
				x = -x;
			} else if (ShowPos) {
				putChar('+');
			}
			long double roundingFactor = 0.5000001l * std::pow(10.0l, -static_cast<int>(Precision));
			x += roundingFactor;
			if (Scientific) {
				int exponent = 0;
				if (x != 0) {
					exponent = static_cast<int>(std::log10(x));
					x /= std::pow(10.0, exponent);
				}
				putChar('0' + static_cast<int>(x));
				x -= static_cast<int>(x);
				if (Precision > 0) {
					if (ShowPoint) putChar('.');
					for (unsigned int i = 0; i < Precision; ++i) {
						x *= 10;
						putChar('0' + static_cast<int>(x));
						x -= static_cast<int>(x);
					}
				}
				putChar(Case ? 'E' : 'e');
				if (exponent < 0) {
					putChar('-');
					exponent = -exponent;
				} else {
					putChar('+');
				}
				if (exponent < 10) putChar('0');
				outputInteger(exponent);
			} else if (Fixed) {
				T integerPart = std::floor(x);
				T fractionalPart = x - integerPart;
				outputInteger(static_cast<long long>(integerPart));
				if (Precision > 0) {
					if (ShowPoint) putChar('.');
					for (unsigned int i = 0; i < Precision; ++i) {
						fractionalPart *= 10;
						int digit = static_cast<int>(fractionalPart);
						putChar('0' + digit);
						fractionalPart -= digit;
					}
				}
			} else {
				if (x == 0) {
					putChar('0');
					if (Precision > 0) {
						if (ShowPoint) putChar('.');
						for (unsigned int i = 0; i < Precision; ++i) putChar('0');
					}
					return;
				}
				int exponent = static_cast<int>(std::log10(x));
				if (exponent < -4 || (exponent > Precision && x != 0)) {
					x /= std::pow(10.0, exponent);
					putChar('0' + static_cast<int>(x));
					x -= static_cast<int>(x);
					if (Precision > 0) {
						if (ShowPoint) putChar('.');
						for (unsigned int i = 0; i < Precision; ++i) {
							x *= 10;
							putChar('0' + static_cast<int>(x));
							x -= static_cast<int>(x);
						}
					}
					putChar(Case ? 'E' : 'e');
					if (exponent < 0) {
						putChar('-');
						exponent = -exponent;
					} else {
						putChar('+');
					}
					if (exponent < 10) putChar('0');
					outputInteger(exponent);
				} else {
					if (exponent < 0) {
						putChar('0');
						if (ShowPoint) putChar('.');
						for (int i = -1; i > exponent; --i) putChar('0');
						T scaled = x * std::pow(10.0, -exponent - 1);
						for (unsigned int i = 0; i < static_cast<unsigned int>(-exponent - 1) + Precision; ++i) {
							scaled *= 10;
							putChar('0' + static_cast<int>(scaled));
							scaled -= static_cast<int>(scaled);
						}
					} else {
						T integerPart = std::floor(x);
						T fractionalPart = x - integerPart;
						outputInteger(static_cast<long long>(integerPart));
						if (Precision > 0) {
							if (ShowPoint) putChar('.');
							for (unsigned int i = 0; i < Precision; ++i) {
								fractionalPart *= 10;
								int digit = static_cast<int>(fractionalPart);
								putChar('0' + digit);
								fractionalPart -= digit;
							}
						}
					}
				}
			}
		}
		inline __attribute__((always_inline)) void outputString(const char* str) noexcept {
			size_t len = strlen(str);
			if (Width > len) {
				size_t padding = Width - len;
				if (!Align) {
					for (size_t i = 0; i < padding; ++i) putChar(Fill);
				}
				while (*str) putChar(*str++);
				if (Align) {
					for (size_t i = 0; i < padding; ++i) putChar(Fill);
				}
			} else {
				while (*str) putChar(*str++);
			}
		}
	public:
		explicit ostream(size_t bufSize = DEFAULT_BUFFER_SIZE) : buffer(new char[bufSize]), current(buffer), end(buffer + bufSize), bufferSize(bufSize) {}
		~ostream() noexcept {
			flushBuffer();
			delete[] buffer;
		}
		constexpr inline __attribute__((always_inline)) void flush() noexcept {
			flushBuffer();
			fflush(stdout);
		}
		inline __attribute__((always_inline)) ostream& setbase(unsigned int b) noexcept {
			Base = b;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& setprecision(unsigned int p) noexcept {
			Precision = p;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& setw(unsigned int w) noexcept {
			Width = w;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& setfill(char c) noexcept {
			Fill = c;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& left() noexcept {
			Align = true;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& right() noexcept {
			Align = false;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& showpos() noexcept {
			ShowPos = true;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& noshowpos() noexcept {
			ShowPos = false;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& showbase() noexcept {
			ShowBase = true;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& noshowbase() noexcept {
			ShowBase = false;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& uppercase() noexcept {
			Case = true;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& nouppercase() noexcept {
			Case = false;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& scientific() noexcept {
			Scientific = true;
			Fixed = false;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& fixed() noexcept {
			Fixed = true;
			Scientific = false;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& defaultfloat() noexcept {
			Scientific = false;
			Fixed = false;
			return *this;
		}
		inline __attribute__((always_inline)) ostream& boolalpha() noexcept {
		    BoolAlpha = true;
		    return *this;
		}
		inline __attribute__((always_inline)) ostream& noboolalpha() noexcept {
		    BoolAlpha = false;
		    return *this;
		}
		inline __attribute__((always_inline)) ostream& showpoint() noexcept {
		    ShowPoint = true;
		    return *this;
		}
		inline __attribute__((always_inline)) ostream& noshowpoint() noexcept {
		    ShowPoint = false;
		    return *this;
		}
		template <typename T, typename std::enable_if_t<std::is_integral_v<T>, bool> = true>
		inline __attribute__((always_inline)) ostream & operator<<(T x) noexcept {
			outputInteger(x);
			return *this;
		}
		template <typename T, typename std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
		inline __attribute__((always_inline)) ostream & operator<<(T x) noexcept {
			outputFloating(x);
			return *this;
		}
		inline __attribute__((always_inline)) ostream& operator<<(const std::string& str) noexcept {
			outputString(str.c_str());
			return *this;
		}
		inline __attribute__((always_inline)) ostream& operator<<(char ch) noexcept {
			putChar(ch);
			return *this;
		}
		inline __attribute__((always_inline)) ostream& operator<<(const char* str) noexcept {
			outputString(str);
			return *this;
		}
		inline __attribute__((always_inline)) ostream& operator<<(bool b) noexcept {
		    if (BoolAlpha) {
		        outputString(b ? (Case ? "TRUE" : "true") : (Case ? "FALSE" : "false"));
		    } else {
		        putChar(b ? '1' : '0');
		    }
		    return *this;
		}
		inline __attribute__((always_inline)) ostream& operator<<(const void* p) noexcept {
		    if (p == nullptr) {
		        outputString("nullptr");
		    } else {
		        unsigned int oldBase = Base;
		        bool oldShowBase = ShowBase;
		        Base = 16;
		        ShowBase = true;
		        outputInteger(reinterpret_cast<uintptr_t>(p));
		        Base = oldBase;
		        ShowBase = oldShowBase;
		    }
		    return *this;
		}
		inline __attribute__((always_inline)) ostream& operator<<(ostream & (*manip)(ostream&)) noexcept {
			return manip(*this);
		}
	};
	inline __attribute__((always_inline)) ostream& flush(ostream& os) noexcept {
		os.flush();
		return os;
	}
	inline __attribute__((always_inline)) ostream& left(ostream& os) noexcept {
		os.left();
		return os;
	}
	inline __attribute__((always_inline)) ostream& right(ostream& os) noexcept {
		os.right();
		return os;
	}
	inline __attribute__((always_inline)) ostream& showpos(ostream& os) noexcept {
		os.showpos();
		return os;
	}
	inline __attribute__((always_inline)) ostream& noshowpos(ostream& os) noexcept {
		os.noshowpos();
		return os;
	}
	inline __attribute__((always_inline)) ostream& showbase(ostream& os) noexcept {
		os.showbase();
		return os;
	}
	inline __attribute__((always_inline)) ostream& noshowbase(ostream& os) noexcept {
		os.noshowbase();
		return os;
	}
	inline __attribute__((always_inline)) ostream& uppercase(ostream& os) noexcept {
		os.uppercase();
		return os;
	}
	inline __attribute__((always_inline)) ostream& nouppercase(ostream& os) noexcept {
		os.nouppercase();
		return os;
	}
	inline __attribute__((always_inline)) ostream& scientific(ostream& os) noexcept {
		os.scientific();
		return os;
	}
	inline __attribute__((always_inline)) ostream& fixed(ostream& os) noexcept {
		os.fixed();
		return os;
	}
	inline __attribute__((always_inline)) ostream& defaultfloat(ostream& os) noexcept {
		os.defaultfloat();
		return os;
	}
	inline __attribute__((always_inline)) ostream& boolalpha(ostream& os) noexcept {
	    os.boolalpha();
	    return os;
	}
	inline __attribute__((always_inline)) ostream& noboolalpha(ostream& os) noexcept {
	    os.noboolalpha();
	    return os;
	}
	inline __attribute__((always_inline)) ostream& showpoint(ostream& os) noexcept {
	    os.showpoint();
	    return os;
	}
	inline __attribute__((always_inline)) ostream& noshowpoint(ostream& os) noexcept {
	    os.noshowpoint();
	    return os;
	}
	class setbase {
		unsigned int base;
	public:
		explicit setbase(unsigned int p) noexcept : base(p) {}
		friend ostream& operator<<(ostream& os, const setbase& manip) noexcept {
			return os.setbase(manip.base);
		}
	};
	class setprecision {
		unsigned int precision;
	public:
		explicit setprecision(unsigned int p) noexcept : precision(p) {}
		friend ostream& operator<<(ostream& os, const setprecision& manip) noexcept {
			return os.setprecision(manip.precision);
		}
	};
	class setw {
		unsigned int width;
	public:
		explicit setw(unsigned int w) noexcept : width(w) {}
		friend ostream& operator<<(ostream& os, const setw& manip) noexcept {
			return os.setw(manip.width);
		}
	};
	class setfill {
		char fill;
	public:
		explicit setfill(char c) noexcept : fill(c) {}
		friend ostream& operator<<(ostream& os, const setfill& manip) noexcept {
			return os.setfill(manip.fill);
		}
	};
	class hex {
	public:
		friend ostream& operator<<(ostream& os, const hex&) noexcept {
			return os.setbase(16).showbase();
		}
	};
	class oct {
	public:
		friend ostream& operator<<(ostream& os, const oct&) noexcept {
			return os.setbase(8).showbase();
		}
	};
	class bin {
	public:
		friend ostream& operator<<(ostream& os, const bin&) noexcept {
			return os.setbase(2).showbase();
		}
	};
	class dec {
	public:
		friend ostream& operator<<(ostream& os, const dec&) noexcept {
			return os.setbase(10).noshowbase();
		}
	};
	inline __attribute__((always_inline)) ostream& endl(ostream& os) noexcept {
		os << '\n';
		os.flush();
		return os;
	}
	inline __attribute__((always_inline)) ostream& ends(ostream& os) noexcept {
		os << '\0';
		return os;
	}
	inline istream cin;
	inline ostream cout;
}
