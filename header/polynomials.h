#ifndef __ALGORITHMIC_ALGEBRA_POLYNOMIALS_H__
#define __ALGORITHMIC_ALGEBRA_POLYNOMIALS_H__

#include <vector>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <type.h>
#include <set>
#include <numeric>
#include <type_traits>
#include <functional>

// template<typename field, std::enable_if_t<is_field<field>::value, bool> = true>
// std::pair<std::vector<field>, bool> detach(const std::vector<field> &s, const field &p)
// {
//     if(p.is_zero())
//     {
//         return {std::vector<field>(s.size()), true};
//     }
//     bool all_zero(true); std::size_t idx = 0;
//     for(std::size_t i = 0; i < s.size(); ++ i) if(!s[i].is_zero())
//     {
//         idx = i; all_zero = false;
//         break;
//     }
//     if(all_zero)
//     {
//         return {{}, false};
//     }
//     std::vector<field> res(s.size());
//     res[idx] = s[idx].inverse() * p;
//     return {res, true};
// }
template<typename field>
class ideal_helper<field, std::enable_if_t<is_field<field>::value>> : std::true_type
{
public:
    const std::vector<field> &s;
    
    ideal_helper(const std::vector<field> &s0):s(s0)
    {
    }

    std::pair<std::vector<field>, bool> detach(const field &p) const
    {
        if(p.is_zero())
        {
            return {std::vector<field>(s.size()), true};
        }
        bool all_zero(true); std::size_t idx = 0;
        for(std::size_t i = 0; i < s.size(); ++ i) if(!s[i].is_zero())
        {
            idx = i; all_zero = false;
            break;
        }
        if(all_zero)
        {
            return {{}, false};
        }
        std::vector<field> res(s.size());
        res[idx] = s[idx].inverse() * p;
        return {res, true};
    }
};

template<typename field>
class syzygy_helper<field, std::enable_if_t<is_field<field>::value>>
{
private:
    std::vector<std::size_t> indices;
    std::size_t nonzero_count;
public:
    std::vector<std::vector<field>> basis;
    std::vector<field> s;
    std::int64_t base_num;
    syzygy_helper(const std::vector<field> &_s) : s(_s)
    {
        std::vector<std::size_t> idx;
        for(std::size_t i = 0; i < s.size(); ++ i) if(!s[i].is_zero())
        {
            idx.emplace_back(i);
        }
        nonzero_count = 0;
        for(std::size_t i = 0; i + 1 < idx.size(); ++ i)
        {
            ++ nonzero_count;
            std::size_t x(idx[i]), y(idx[i + 1]);
            std::vector<field> row(s.size());
            row[x] = s[x].inverse();
            row[y] = - s[y].inverse();
            indices.emplace_back(x);
            basis.emplace_back(row);
        }
        for(std::size_t i = 0; i < s.size(); ++ i) if(s[i].is_zero())
        {
            indices.emplace_back(i);
        }
        base_num = basis.size();
    }
    std::pair<std::vector<field>, bool> solve_syzygy(std::vector<field> p)
    {
        p.resize(s.size()); field sum;
        for(std::size_t i = 0; i < s.size(); ++ i) sum += s[i] * p[i];
        if(!sum.is_zero()) return {{}, false};
        field total; std::vector<field> res(basis.size());
        for(std::size_t i = 0; i < basis.size(); ++i)
        {
            if(i < nonzero_count)
            {
                std::size_t x(indices[i]);
                total = total + s[x] * p[x];
                res[i] = total;
            }
            else
            {
                std::size_t x(indices[i]);
                res[i] = p[x];
            }
        }
        return res;
    }
};

class Hterm
{
public:
    std::vector<std::int64_t> powers;
    Hterm() = default;
    template<typename T> Hterm(const T& v) : powers(v){};

    friend Hterm operator * (const Hterm &a, const Hterm &b)
    {
        Hterm r;
        r.powers.resize(a.powers.size());
        std::copy(a.powers.begin(), a.powers.end(), r.powers.begin());
        std::size_t l(b.powers.size());
        for(std::size_t i(0);i<l;++i) r.powers[i] += b.powers[i];
        return r;
    }

    friend Hterm operator / (const Hterm &a, const Hterm &b)
    {
        Hterm r(a);
        std::size_t l = std::max(a.powers.size(), b.powers.size());
        r.powers.resize(l);
        for(std::size_t i(0); i < b.powers.size(); ++ i) r.powers[i] = std::max<std::int64_t>(r.powers[i] - b.powers[i], 0);
        return r;
    }

    bool is_zero() const
    {
        for(const auto &v: powers) if(v != 0) return false;
        return true;
    }
    decltype(powers) zero() const {return decltype(powers)(powers.size());}

    friend bool operator <= (Hterm a, Hterm b)
    {
        std::size_t l = std::max(a.powers.size(), b.powers.size());
        a.powers.resize(l);b.powers.resize(l);
        bool capable(true);
        for(std::size_t i(0); i < l; ++ i) capable &= (a.powers[i] <= b.powers[i]);
        return capable;
    }

    friend bool operator == (Hterm a, Hterm b)
    {
        std::size_t l = std::max(a.powers.size(), b.powers.size());
        a.powers.resize(l);b.powers.resize(l);
        bool capable(true);
        for(std::size_t i(0); i < l; ++ i) capable &= (a.powers[i] == b.powers[i]);
        return capable;
    }

    friend bool operator < (Hterm a, Hterm b)
    {
        return (a <= b) && (!(a == b));
    }

    friend Hterm hterm_max(const Hterm &a, const Hterm &b)
    {
        return a * (b / a);
    }
};

namespace term_ordering
{
    bool equal(const Hterm &a, const Hterm &b)
    {
        std::size_t len = std::min(a.powers.size(), b.powers.size());
        bool result = true;
        for(std::size_t i(0); i < len; ++i) result &= (a.powers[i] == b.powers[i]);
        for(std::size_t i(len); i < a.powers.size(); ++i) result &= (a.powers[i] == 0);
        for(std::size_t i(len); i < b.powers.size(); ++i) result &= (b.powers[i] == 0);
        return result;
    }

    bool lexicographic_ordering_func(const Hterm &a, const Hterm &b)
    {
        std::size_t len = std::min(a.powers.size(), b.powers.size());
        for(std::size_t i(0); i < len; ++i) if(a.powers[i] != b.powers[i])
            return a.powers[i] < b.powers[i];
        std::int64_t total_a = std::accumulate(a.powers.begin() + len, a.powers.end(), 0);
        std::int64_t total_b = std::accumulate(b.powers.begin() + len, b.powers.end(), 0);
        return total_a < total_b;
    }

    bool lexicographic_ordering_inferior_to_sum_func(const Hterm &a, const Hterm &b)
    {
        std::int64_t sum_a = 0, sum_b = 0;
        sum_a = std::accumulate(a.powers.begin(), a.powers.end(), 0);
        sum_b = std::accumulate(b.powers.begin(), b.powers.end(), 0);
        if(sum_a != sum_b) return sum_a < sum_b;
        else return lexicographic_ordering_func(a,b);
    }
    struct lexicographic_ordering
    {
        bool operator()(const Hterm &a, const Hterm &b) const {return lexicographic_ordering_func(a, b);}
        // bool le(const Hterm &a, const Hterm &b) const {return equal(a, b) || lexicographic_ordering_func(a,b);}
    };
    struct lexicographic_ordering_inferior_to_sum
    {
        bool operator()(const Hterm &a, const Hterm &b) const {return lexicographic_ordering_inferior_to_sum_func(a, b);}
        // bool le(const Hterm &a, const Hterm &b) const {return equal(a, b) || lexicographic_ordering_inferior_to_sum_func(a,b);}
    };
}


template<typename base_ring, typename func, std::enable_if_t<is_unitary_ring_v<base_ring>, bool> = true>
class multivariate_polynomial_ring
{
private:
    struct compare_func
    {
        bool operator()(const std::pair<Hterm, base_ring> &a, const std::pair<Hterm, base_ring> &b) const {return func{}(a.first, b.first);}
    };
public:
    std::set<std::pair<Hterm, base_ring>, compare_func> data;
    void insert(const std::pair<Hterm, base_ring> &a)
    {
        auto it = data.find(a);
        if(it == data.end())
        {
            if(data.size() == 1 && data.begin() -> second.is_zero())
            {
                data.clear();
                data.insert(a);
            }
            else data.insert(a);
        }
        else
        {
            it -> second = it -> second + a -> second;
            if(it -> second.is_zero())
            {
                if(data.size() == 1)
                {
                    data.clear();
                    data.insert({a -> first.zero(), base_ring::get_zero()});
                }
                else data.erase(it);
            }
        }
    }
    multivariate_polynomial_ring() : data({std::make_pair(Hterm({0}), base_ring::get_zero())})
    {
    }

    template<typename T>
    multivariate_polynomial_ring(const T &_data)
    {
        data = _data;
    }

    static multivariate_polynomial_ring<base_ring, func> get_zero()
    {
        return multivariate_polynomial_ring<base_ring, func>();
    }

    static multivariate_polynomial_ring<base_ring, func> get_identity()
    {
        return multivariate_polynomial_ring<base_ring, func>({Hterm({0}), base_ring::get_identity()});
    }

    bool is_zero() const
    {
        if(data.size() == 1 && data.begin() -> second.is_zero())
        {
            return true;
        }
        else return false;
    }
};

template<typename base_ring, typename func, std::enable_if_t<is_unitary_ring_v<base_ring>, bool> = true>
multivariate_polynomial_ring<base_ring,func> operator +(const multivariate_polynomial_ring<base_ring,func> &a, const multivariate_polynomial_ring<base_ring,func> &b)
{
    multivariate_polynomial_ring<base_ring, func> res;
    res.data = (decltype(res.data))(a.data);
    for(const auto &v : b.data)
        res.insert(v);
    return res;
}

template<typename base_ring, typename func, std::enable_if_t<is_unitary_ring_v<base_ring>, bool> = true>
multivariate_polynomial_ring<base_ring,func> operator -(const multivariate_polynomial_ring<base_ring,func> &a, const multivariate_polynomial_ring<base_ring,func> &b)
{
    multivariate_polynomial_ring<base_ring, func> res;
    res.data = (decltype(res.data))(a.data);
    for(auto v : b.data)
    {
        v.second = - v.second;
        res.insert(v);
    }
    return res;
}

template<typename base_ring, typename func, std::enable_if_t<is_unitary_ring_v<base_ring>, bool> = true>
multivariate_polynomial_ring<base_ring,func> operator *(const multivariate_polynomial_ring<base_ring,func> &a, const multivariate_polynomial_ring<base_ring,func> &b)
{
    multivariate_polynomial_ring<base_ring, func> res;
    for(const auto &v_a : a.data)
    {
        for(const auto &v_b : b.data)
            res.insert({v_a.first * v_b.first, v_a.second * v_b.second});
    }
    return res;
}

template<typename base_ring, typename func> struct is_ring<multivariate_polynomial_ring<base_ring, func>> : std::true_type {};
template<typename base_ring, typename func> struct has_identity<multivariate_polynomial_ring<base_ring, func>> : std::true_type {};

template<typename base_ring, typename func> struct is_strongly_computable<multivariate_polynomial_ring<base_ring, func>,
    std::enable_if_t<is_strongly_computable<base_ring>::value, void> > : std::true_type {};


template<typename base_ring, typename func, std::enable_if_t<is_strongly_computable<base_ring>::value, bool> = true>
std::vector<std::pair<multivariate_polynomial_ring<base_ring, func>, std::vector<multivariate_polynomial_ring<base_ring, func>>>>
grobner_basis(const std::vector<multivariate_polynomial_ring<base_ring,func>> &g0)
{
    using poly = multivariate_polynomial_ring<base_ring, func>;
    std::vector<std::pair<poly, std::vector<poly>>> g;
    std::int64_t l(g0.size());
    if(!l) return g;
    for(std::int64_t i(0); i < l; ++i) if(!g0[i].is_zero())
    {
        std::vector<poly> coef(l);
        coef[i].insert(std::make_pair(Hterm({0}), base_ring::get_identity()));
        g.emplace_back(g0[i], coef[i]);
    }
    auto g_reduction = [&](poly &a) -> std::vector<poly>
    {
        std::vector<poly> coefs(l);
        while(!a.is_zero())
        {
            std::vector<base_ring> c;
            std::vector<std::int64_t> idx;
            for(std::int64_t i(0); i < g.size(); ++i) if(g[i].first.data.rbegin()->first <= a.data.rbegin()->first)
            {
                c.emplace_back(g[i].first.data.rbegin() -> second);
                idx.emplace_back(i);
            }
            ideal_helper<base_ring> c_ideal(c);
            auto [res, cap] = c_ideal.detach(a.data.rbegin()->second);
            if(!cap) break;
            for(std::int64_t i(0); i < idx.size(); ++i)
            {
                poly multiplicator({std::make_pair(a.data.rbegin()->first / g[idx[i]].first.data.rbegin()->first, res[i])});
                for(std::int64_t j(0); j < l; ++j)
                {
                    coefs[j] = coefs[j] - g[idx[i]].second[j] * multiplicator;
                }
                a = a - g[idx[i]].first * multiplicator;
            }
        }
        return coefs;
    };
    std::function<bool(const std::vector<std::int64_t>&, std::int64_t pos)> dfs = [&](std::vector<std::int64_t> idx, std::int64_t pos) -> bool
    {
        if(pos >= g.size()) return false;
        if(dfs(idx, pos + 1)) return true;
        if(!g[pos].first.is_zero())
        {
            idx.emplace_back(pos);
            if(idx.size() > 1)
            {
                std::vector<base_ring> c;
                Hterm lcm;
                for(std::int64_t i : idx)
                {
                    c.emplace_back(g[i].first.rbegin()->second);
                    lcm = hterm_max(lcm, g[i].first.rbegin()->first);
                }
                syzygy_helper<base_ring> c_syzygy(c);
                std::int64_t li(idx.size());
                for(std::int64_t i(0); i < c_syzygy.base_num; ++ i)
                {
                    poly a;
                    for(std::int64_t j(0); j < li; ++ j)
                    {
                        poly multiplicator({std::make_pair(lcm / g[idx[j]].first.rbegin()->first, c_syzygy.basis[i][j])});
                        a = a + g[idx[j]].second * multiplicator;
                    }
                    auto coefs = g_reduction(a);
                    if(!a.is_zero())
                    {
                        for(std::int64_t j(0); j < li; ++j)
                        {
                            poly multiplicator({std::make_pair(lcm / g[idx[j]].first.rbegin()->first, c_syzygy.basis[i][j])});
                            for(std::int64_t k(0); k < l; ++k)
                            {
                                coefs[k] = coefs[k] + multiplicator * g[idx[j]].second[k];
                            }
                        }
                        g.emplace_back(a, coefs);
                        return true;
                    }
                }
            }
            if(dfs(idx, pos + 1)) return true;
        }
        return false;
    };
    while(true)
    {
        if(!dfs(std::vector<std::int64_t>(), 0))
        {
            break;
        }
    }
    return g;
}

#endif