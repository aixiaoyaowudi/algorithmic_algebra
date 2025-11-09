#ifndef __ALGORITHMIC_ALGEBRA_TYPE_H__
#define __ALGORITHMIC_ALGEBRA_TYPE_H__

#include <utility>

template<typename T> struct is_field:std::false_type {};

template<typename T> struct is_ring:std::false_type {};

template<typename T> struct has_identity:std::false_type {};

template <typename T>
struct is_unitary_ring : std::bool_constant<std::conjunction_v<is_ring<T>, has_identity<T>> > {};

template <typename T>
constexpr bool is_unitary_ring_v = is_unitary_ring<T>::value;

template<typename T, typename = void>
struct is_strongly_computable:std::false_type {};

template<typename T>
struct is_strongly_computable<T, std::enable_if_t<is_field<T>::value>>:std::true_type {};

template <typename F>
auto make_functor(F&& f)
{
    return [f = std::forward<F>(f)]<typename... Args>(Args&&... args) -> decltype(auto)
    {
        return f(std::forward<Args>(args)...);
    };
}

template<typename T, typename = void>
class syzygy_helper:public std::false_type
{
};

template<typename T, typename = void>
class ideal_helper:public std::false_type
{
};

template<typename T, std::enable_if_t<ideal_helper<T>::value, bool> = true>
class ideal_solver:public ideal_helper<T>
{
    using ideal_helper<T>::ideal_helper;
};

template<typename T, std::enable_if_t<syzygy_helper<T>::value, bool> = true>
class syzygy_solver:public syzygy_helper<T>
{
    using syzygy_helper<T>::syzygy_helper;
};

#endif