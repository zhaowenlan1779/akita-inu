add_library(src STATIC
    ASHE.cpp
    ASHE.h
    big_int.h
    modular_arithmetic.hpp
    ring_polynomial.hpp
)

target_link_libraries(src PRIVATE kissfft)
