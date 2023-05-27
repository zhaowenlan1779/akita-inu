add_library(src STATIC
    ASHE.cpp
    ASHE.h
    big_int.h
    common.h
    fast_evaluation.cpp
    fast_evaluation.h
    interpolation.cpp
    interpolation.h
    kissfft.hpp
    modular_arithmetic.hpp
    ring_polynomial.hpp
)

target_link_libraries(src PUBLIC boost)
target_link_libraries(src PRIVATE thread-pool)
