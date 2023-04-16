#pragma GCC optimize("O3")
#include <bits/stdc++.h>

#define ll long long
#define ld double
#define pb push_back
#define ff first
#define ss second
#define files(in, out) freopen(in, "r", stdin); freopen(out, "w", stdout)
#define rep(i, n) for (int i = 0; i < (int)(n); ++i)
#define all(a) a.begin(), a.end()
#define pii pair<int, int>

using namespace std;

const double eps = 1e-5;

using comp = complex<double>;

enum Sign {
    NEG = -1,
    ZERO = 0,
    POS = 1
};

class BigInteger {
private:
    vector<int> nums;
    int base;
    Sign sign = ZERO;

    Sign get_sign(int n) const {
        if (n > 0) return POS;
        if (n < 0) return NEG;
        return ZERO;
    }

    int get_mod(int n) {
        if (n >= 0) return n % base;
        return (-((-n) % base) + base) % base;
    }

    int get_div(int n)  {
        if (n >= 0) return n / base;
        return -(((-n) + base - 1) / base);
    }

    void change_sign() {
        if (sign == NEG) sign = POS;
        else if (sign == POS) sign = NEG;
    }

    void remove_lead_zeroes() {
        while (nums.size() > 1 && nums.back() == 0) {
            nums.pop_back();
        }
    }

    static void remove_lead_zeroes(string& s) {
        while (s.size() > 1 && s.back() == '0') {
            s.pop_back();
        }
    }

    void plus_nums(const vector<int>& numbers, short dop_sign) {
        for (size_t i = 0; i < max(numbers.size(), nums.size()); ++i) {
            if (i >= nums.size()) nums.push_back(0);
            nums[i] *= sign;
            if (i < numbers.size()) nums[i] += dop_sign * numbers[i];
        }
        remove_lead_zeroes();
        sign = get_sign(nums.back());
        for (size_t i = 0; i < nums.size(); ++i) {
            nums[i] *= sign;
        }
        for (size_t i = 0; i < nums.size(); ++i) {
            int add = get_div(nums[i]);
            if (add != 0) {
                if (i + 1 >= nums.size()) nums.push_back(0);
                nums[i + 1] += add;
            }
            nums[i] = get_mod(nums[i]);
        }
        remove_lead_zeroes();
    }

    static void fft (vector<comp>& a, bool invert) {
        int n = static_cast<int>(a.size());
        int k = 0;
        for (int i = 1; i < n; ++i) {
            int bit = n >> 1;
            for (;k >= bit; bit >>= 1){
                k -= bit;
            }
            k += bit;
            if (i < k) {
                swap(a[i], a[k]);
            }
        }
        double pi = acos(-1);
        for (int len = 2; len <= n; len *= 2) {
            double ang = 2 * pi / len * (invert ? -1 : 1);
            comp wi(cos(ang), sin(ang));
            for (int i = 0; i < n; i += len) {
                comp w(1);
                for (int j = 0; j < len / 2; ++j) {
                    comp u = a[i + j];
                    comp v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wi;
                }
            }
        }
        if (invert){
            for (int i = 0; i < n; ++i) {
                a[i] /= n;
            }
        }
    }

    void multiply(vector<int>& a, const vector<int>& b) {
        vector<comp> fa(a.begin(), a.end()), fb(b.begin(), b.end());
        int n = 1;
        while (n < static_cast<int>(max(a.size(), b.size()))) n *= 2;
        n *= 2;
        fa.resize(n);
        fb.resize(n);
        fft(fa, false);
        fft(fb, false);
        for (int i = 0; i < n; ++i) {
            fa[i] *= fb[i];
        }
        fft(fa, true);
        a.clear();
        a.resize(n);
        int carry = 0;
        for (size_t i = 0; i < a.size(); ++i) {
            long long el = static_cast<long long>(fa[i].real() + 0.5);
            el += carry;
            carry = el / base;
            if (carry > 0 && i + 1 >= a.size())
                a.push_back(0);
            el %= base;
            a[i] = el;
        }
    }

    void shift(int n) {
        nums.push_back(0);
        for (int i = static_cast<int>(nums.size()) - 2; i >= 0; --i) {
            nums[i + 1] = nums[i];
        }
        nums[0] = n;
        remove_lead_zeroes();
        if (nums.size() == 1 && nums.back() == 0) {
            sign = ZERO;
        } else {
            sign = POS;
        }
    }

    string digitToString(int n) const {
        string s = "";
        if (n == 0) {
            string res = "";
            int temp = base;
            while (temp > 1) {
                res += '0';
                temp /= 2;
            }
            return res;
        }
        int old = n;
        while (n > 0) {
            s += (n % 2) + '0';
            n /= 2;
        }
        n = old;
        while (2 * n < base) {
            s += '0';
            n *= 2;
        }
        return s;
    }

public:

    BigInteger() : nums({0}), sign(ZERO) {}

    BigInteger(int base) : base(base) {
        sign = ZERO;
        nums = {0};
    }

    BigInteger(int n, int base) : base(base) {
        if (n == 0) {
            sign = ZERO;
            nums = {0};
            return;
        }
        if (n < 0) {
            sign = NEG;
            n = -n;
        } else {
            sign = POS;
        }
        while (n > 0) {
            nums.push_back(n % base);
            n /= base;
        }
    }

    friend bool operator<(const BigInteger&, const BigInteger&);

    explicit operator bool() const {
        return sign != 0;
    }

    int get_sign() const {
        return sign;
    }

    bool is_even() const {
        return nums[0] % 2 == 0;
    }

    BigInteger& operator+=(const BigInteger& number) {
        BigInteger copy = number;
        plus_nums(copy.nums, copy.sign);
        return *this;
    }

    BigInteger& operator-=(const BigInteger& number) {
        BigInteger copy = number;
        plus_nums(copy.nums, -copy.sign);
        return *this;
    }

    BigInteger& operator*=(const BigInteger& number) {
        if (number.nums.size() == 1) {
            int num = number.nums[0];
            num *= static_cast<int>(number.sign);
            *this *= num;
            return *this;
        }
        BigInteger copy = number;
        multiply(nums, copy.nums);
        sign = static_cast<Sign>(static_cast<int>(copy.sign) * static_cast<int>(sign));
        remove_lead_zeroes();
        return *this;
    }

    BigInteger& operator*=(int num) {
        if (num < 0) {
            num *= -1;
            change_sign();
        }
        long long add = 0;
        for (size_t i = 0; i < nums.size(); ++i) {
            long long el = nums[i] * num;
            if (add > 0) el += add;
            add = el / base;
            if (add > 0 && i + 1 == nums.size()) nums.push_back(0);
            el %= base;
            nums[i] = el;
        }
        remove_lead_zeroes();
        if (nums.size() == 1 && nums[0] == 0) sign = ZERO;
        return *this;
    }

    BigInteger& operator/=(const BigInteger& number) {
        if (number.nums.size() == 1) {
            int num = number.nums[0];
            num *= static_cast<int>(number.sign);
            *this /= num;
            return *this;
        }
        if (number.sign == 0) {
            sign = ZERO;
            nums = {0};
            return *this;
        }
        BigInteger remain;
        for (int i = static_cast<int>(nums.size()) - 1; i >= 0; --i) {
            remain.shift(nums[i]);
            int last = 0;
            int l = 0;
            int r = base;
            while (r - l > 1) {
                int mid = (l + r) / 2;
                if (number > 0) {
                    if (number * mid <= remain) l = mid;
                    else r = mid;
                } else {
                    if ((-number) * mid <= remain) l = mid;
                    else r = mid;
                }
            }
            last = l;
            if (number > 0) remain -= number * last;
            else remain += number * last;
            nums[i] = last;
        }
        sign = static_cast<Sign>(static_cast<int>(number.sign) * static_cast<int>(sign));
        remove_lead_zeroes();
        if (nums.size() == 1 && nums.back() == 0) sign = ZERO;
        return *this;
    }

    BigInteger& operator/=(int num) {
        if (num < 0) {
            change_sign();
            num *= -1;
        }
        for (int i = static_cast<int>(nums.size()) - 1; i >= 0; --i) {
            if (i != 0) nums[i - 1] += base * (nums[i] % num);
            nums[i] /= num;
        }
        remove_lead_zeroes();
        if (nums.size() == 1 && nums[0] == 0) sign = ZERO;
        return *this;
    }

    BigInteger& operator%=(const BigInteger& number) {
        *this = (*this - (*this / number) * number);
        return *this;
    }

    BigInteger& operator++() {
        *this += 1;
        return *this;
    }

    BigInteger& operator--() {
        *this -= 1;
        return *this;
    }

    string toString() const {
        string res = "";
        for (size_t i = 0; i < nums.size(); ++i) {
            res += digitToString(nums[i]);
        }
        remove_lead_zeroes(res);
        if (sign == NEG) res += '-';
        reverse(res.begin(), res.end());
        return res;
    }

    BigInteger transform_recursive(int to_base, int n, int k, const vector<BigInteger> pows) {
        if (n == 1) {
            return BigInteger(nums[0], to_base);
        }
        BigInteger left(base), right(base);
        left.nums.resize(n / 2);
        right.nums.resize(n / 2);
        rep(i, n) {
            if (i < n / 2) left.nums[i] = nums[i];
            else right.nums[i - n / 2] = nums[i];
        }
        left = left.transform_recursive(to_base, n / 2, k - 1, pows);
        right = right.transform_recursive(to_base, n / 2, k - 1, pows);
        return left + right * pows[k];
    }

    BigInteger transform(int to_base) {
        // BigInteger pow(base);
        int n = 1;
        int k = 0;
        while (n < nums.size()) {
            n *= 2;
            k++;
        }
        nums.resize(n);
        vector<BigInteger> pows(23);
        rep(i, k) {
            if (i == 0) pows[i] = BigInteger(base, to_base);
            else pows[i] = pows[i - 1] * pows[i - 1];
        }
        return transform_recursive(to_base, nums.size(), k - 1, pows);
    }

    friend bool operator<(const BigInteger&, const BigInteger&);

    friend bool operator>(const BigInteger&, const BigInteger&);

    friend bool operator==(const BigInteger&, const BigInteger&);

    friend bool operator!=(const BigInteger&, const BigInteger&);

    friend bool operator<=(const BigInteger&, const BigInteger&);

    friend bool operator>=(const BigInteger&, const BigInteger&);

    friend BigInteger operator-(const BigInteger&);

    friend BigInteger operator+(const BigInteger&, const BigInteger&);

    friend BigInteger operator-(const BigInteger&, const BigInteger&);

    friend BigInteger operator*(const BigInteger&, const BigInteger&);

    friend BigInteger operator/(const BigInteger&, const BigInteger&);

    friend BigInteger operator%(const BigInteger&, const BigInteger&);

    friend istream& operator>>(istream&, BigInteger&);

    friend ostream& operator<<(ostream&, const BigInteger&);
};

bool operator<(const BigInteger& first, const BigInteger& second) {
    if (first.sign * static_cast<int>(first.nums.size()) < second.sign * static_cast<int>(second.nums.size())) return true;
    if (first.sign * static_cast<int>(first.nums.size()) > second.sign * static_cast<int>(second.nums.size())) return false;
    for (int i = static_cast<int>(first.nums.size()) - 1; i >= 0; --i) {
        if (first.nums[i] == second.nums[i]) continue;
        return first.sign * first.nums[i] < second.sign * second.nums[i];
    }
    return false;
}

bool operator>(const BigInteger& first, const BigInteger& second) {
    return second < first;
}

bool operator==(const BigInteger& first, const BigInteger& second) {
    return !(first < second) && !(second < first);
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
    return !(first == second);
}

bool operator<=(const BigInteger& first, const BigInteger& second) {
    return first == second || first < second;
}

bool operator>=(const BigInteger& first, const BigInteger& second) {
    return first == second || first > second;
}

BigInteger operator-(const BigInteger& first) {
    BigInteger copy = first;
    copy *= -1;
    return copy;
}

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
    BigInteger res = first;
    res += second;
    return res;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
    BigInteger res = first;
    res -= second;
    return res;
}

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
    BigInteger res = first;
    res *= second;
    return res;
}

BigInteger operator*(const BigInteger& first, int num) {
    BigInteger res = first;
    res *= num;
    return res;
}

BigInteger operator*(int num, const BigInteger& first) {
    return first * num;
}

BigInteger operator/(const BigInteger& first, int num) {
    BigInteger res = first;
    res /= num;
    return res;
}

BigInteger operator/(const BigInteger& first, const BigInteger& second) {
    BigInteger res = first;
    res /= second;
    return res;
}

BigInteger operator%(const BigInteger& first, const BigInteger& second) {
    BigInteger res = first;
    res %= second;
    return res;
}

istream& operator>>(istream& in, BigInteger& number) {
    string s;
    in >> s;
    reverse(s.begin(), s.end());
    if (s.back() == '-') {
        number.sign = NEG;
        s.pop_back();
    } else {
        number.sign = POS;
    }
    number.nums.clear();
    int cur = 0;
    int pow = 1;
    for (size_t i = 0; i < s.size(); ++i) {
        if (pow >= number.base) {
            number.nums.push_back(cur);
            cur = 0;
            pow = 1;
        }
        cur += (s[i] - '0') * pow;
        pow *= 10;
    }
    number.nums.push_back(cur);
    number.remove_lead_zeroes();
    if (number.nums.size() == 1 && number.nums.back() == 0) {
        number.sign = ZERO;
    }
    return in;
}

ostream& operator<<(ostream& out, const BigInteger& number) {
    out << number.toString();
    return out;
}
