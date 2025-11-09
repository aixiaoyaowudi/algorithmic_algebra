#ifndef __ALGORITHMIC_ALGEBRA_CLASSES_H__
#define __ALGORITHMIC_ALGEBRA_CLASSES_H__

#include <gmpxx.h>
#include <iostream>
#include <type.h>
#include <queue>

class Z
{
public:
    mpz_class value;

    Z() : value(0){}
    
    template<typename T>
    Z(const T &_value) : value(_value){}

    Z(const Z &_value) : value(_value.value){}

    template<typename T>
    void set(const T &_value)
    {
        value = _value;
    }
    
    friend std::ostream& operator<<(std::ostream &out, const Z &v)
    {
        out << v.value;
        return out;
    }
    
    static const Z identity, zero_const;

    Z operator -() const
    {
        return -value;
    }

    friend Z operator + (const Z &a,const Z &b)
    {
        return a.value + b.value;
    }

    friend Z operator * (const Z &a,const Z &b)
    {
        return a.value * b.value;
    }
    
    static Z get_identity(){return identity;}
    static Z get_zero(){return zero_const;}
    bool is_zero() const {return mpz_sgn(value.get_mpz_t()) == 0;}
};

class Q
{
public:
    mpq_class value;

    Q() : value(0){}
    
    template<typename T>
    Q(const T &_value) : value(_value){}

    Q(const Q &_value) : value(_value.value){}

    template<typename T>
    void set(const T &_value)
    {
        value = _value;
    }
    
    friend std::ostream& operator<<(std::ostream &out, const Q &v)
    {
        out << v.value;
        return out;
    }
    
    static const Q identity, zero_const;

    Q operator -() const
    {
        return -value;
    }

    friend Q operator + (const Q &a,const Q &b)
    {
        return a.value + b.value;
    }

    friend Q operator * (const Q &a,const Q &b)
    {
        return a.value * b.value;
    }

    static Q get_identity(){return identity;}
    static Q get_zero(){return zero_const;}
    bool is_zero() const {return mpq_sgn(value.get_mpq_t()) == 0;}

    Q inverse() const
    {
        if(is_zero()) return zero_const;
        else return identity.value / value;
    }
};

const Z Z::identity(1), Z::zero_const(0);
const Q Q::identity(1), Q::zero_const(0);

template<> struct is_ring<Z> : std::true_type {};
template<> struct is_ring<Q> : std::true_type {};
template<> struct has_identity<Q> : std::true_type {};
template<> struct has_identity<Z> : std::true_type {};

template<> struct is_field<Q> : std::true_type {};

template<> struct is_strongly_computable<Z, void>:std::true_type {};

namespace tools
{
    mpz_class ex_gcd(const mpz_class &a, const mpz_class &b, mpz_class &x, mpz_class &y)
    {
        if(mpz_sgn(a.get_mpz_t()) == 0)
        {
            x = 0;
            y = 1;
            return b;
        }
        else if(mpz_sgn(b.get_mpz_t()) == 0)
        {
            x = 1;
            y = 0;
            return a;
        }
        else if(a < b)
        {
            mpz_class d = b / a;
            mpz_class g = ex_gcd(a, b - d * a, x, y);
            x -= d * y;
            return g;
        }
        else
        {
            mpz_class d = a / b;
            mpz_class g = ex_gcd(a - d * b, b, x, y);
            y -= d * x;
            return g;
        }
    }
}

template<>
class ideal_helper<Z> : public std::true_type
{
public:
    const std::vector<Z> s;
    
    std::vector<Z> res;
    std::int64_t count;
    mpz_class g;
    ideal_helper(const std::vector<Z> &s0):s(s0),g(0),count(0)
    {
        std::vector<std::int64_t> father(s.size());
        std::vector<mpz_class> value(s.size());
        std::vector<mpz_class> coef(s.size());
        std::vector<int> sgns;
        auto cmp_func = [&](std::size_t a, std::size_t b) -> bool {return value[a] < value[b];};    
        std::priority_queue<std::size_t, std::vector<std::size_t>, decltype(cmp_func)> pq(cmp_func);
        for(std::size_t _(0); _ < s.size(); ++ _)
        {
            const auto &v = s[_];
            auto sgn = mpz_sgn(v.value.get_mpz_t());
            if(sgn > 0) value[_] = v.value, ++ count;
            else if(sgn < 0) value[_] = - v.value, ++ count;
            else value[_] = 0, father[_] = 0;
            sgns.emplace_back(sgn);
            if(sgn != 0)
            {
                father[_] = _;
                pq.push(_);
            }
        }
        while(count > 1)
        {
            std::size_t x = pq.top(); pq.pop();
            std::size_t y = pq.top(); pq.pop();
            mpz_class s, t;
            mpz_class g = tools::ex_gcd(value[x], value[y], s, t);
            std::size_t cur = value.size();
            value.emplace_back(g);
            pq.push(cur);
            father.emplace_back(cur);
            coef.emplace_back();
            coef[x] = s;
            coef[y] = t;
            father[x] = cur;
            father[y] = cur;
            -- count;
        }
        if(count > 0)
        {
            g = value[pq.top()];
        }
        res.resize(s.size());
        for(std::size_t _(0); _ < s.size(); ++ _) if(sgns[_])
        {
            if(sgns[_] > 0) res[_] = Z::get_identity();
            else res[_] = - Z::get_identity();
            for(std::int64_t ptr = _; ptr != father[ptr]; ptr = father[ptr])
            {
                res[_].value *= coef[ptr];
            }
        }
    }
    std::pair<std::vector<Z>, bool> detach(const Z &p) const
    {
        if(p.is_zero())
        {
            return {std::vector<Z>(s.size()), true};
        }
        if(count == 0) return {{}, false};
        mpz_class val = p.value;
        mpz_class d = val / g, r = val - g * d;
        if(mpz_sgn(r.get_mpz_t()) != 0) return {{}, false};
        std::vector<Z> p_res(res);
        Z Zd(d);
        for(std::size_t _(0); _ < s.size(); ++ _)
            p_res[_] = p_res[_] * Zd;
        return {p_res, true};
    }
};

#endif