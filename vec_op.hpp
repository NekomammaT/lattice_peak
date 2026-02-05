#ifndef INCLUDED_vec_op_hpp_
#define INCLUDED_vec_op_hpp_

#include <vector>
#include <cassert>

template <typename T>
struct is_std_vector : std::false_type {};

template <typename T, typename A>
struct is_std_vector<std::vector<T, A>> : std::true_type {};

template <typename T>
inline constexpr bool is_std_vector_v = is_std_vector<T>::value;


template <typename T1, typename T2>
void add_inplace_recursive(T1& a, const T2& b)
{
    if constexpr (is_std_vector_v<T1> && is_std_vector_v<T2>) {
        static_assert(std::is_same_v<
            typename T1::value_type,
            typename T2::value_type>,
            "vector value_type mismatch");

        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) {
            add_inplace_recursive(a[i], b[i]);  // 再帰！
        }
    } else {
        a += b;  // スカラー
    }
}

template <typename T1, typename T2>
std::vector<T1>& operator+=(std::vector<T1>& v1,
                            const std::vector<T2>& v2)
{
    assert(v1.size() == v2.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        add_inplace_recursive(v1[i], v2[i]);
    }
    return v1;
}

template <typename T1, typename T2>
std::vector<T1> operator+(std::vector<T1>& v1,
                            const std::vector<T2>& v2)
{
    auto result = v1;
    result += v2;
    return result;
}


template <typename T, typename Scalar>
void scale_inplace_recursive(T& a, const Scalar& alpha)
{
    if constexpr (is_std_vector_v<T>) {
        for (auto& ai : a) {
            scale_inplace_recursive(ai, alpha);  // 再帰
        }
    } else {
        a *= alpha;  // スカラー
    }
}

template <typename T, typename Scalar>
std::vector<T>& operator*=(std::vector<T>& v, const Scalar& alpha)
{
    for (auto& vi : v)
        scale_inplace_recursive(vi, alpha);
    return v;
}

template <typename T, typename Scalar>
std::vector<T> operator*(const std::vector<T>& v, const Scalar& alpha)
{
    auto out = v;
    out *= alpha;
    return out;
}



#endif
