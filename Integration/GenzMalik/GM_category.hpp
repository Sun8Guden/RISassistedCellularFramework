#pragma once

#include <limits>

namespace GM
{
    template <class T>
    struct genz_malik_category
    {
           static const unsigned value =
                (std::numeric_limits<T>::is_specialized == 0) ? 999 :
                (std::numeric_limits<T>::radix == 2) ?
                (
                    (std::numeric_limits<T>::digits <= std::numeric_limits<float>::digits) && std::is_convertible<float, T>::value ? 0 :
                    (std::numeric_limits<T>::digits <= std::numeric_limits<double>::digits) && std::is_convertible<double, T>::value ? 1 :
                    (std::numeric_limits<T>::digits <= std::numeric_limits<long double>::digits) && std::is_convertible<long double, T>::value ? 2 :
            #ifdef BOOST_HAS_FLOAT128
                    (std::numeric_limits<T>::digits <= 113) && std::is_constructible<__float128, T>::value ? 3 :
            #endif
                    (std::numeric_limits<T>::digits10 <= 110) ? 4 : 999
                ) : (std::numeric_limits<T>::digits10 <= 110) ? 4 : 999;
    };
     
} // namespace GM
