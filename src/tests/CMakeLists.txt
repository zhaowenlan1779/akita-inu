add_executable(tests
    ASHE.cpp
    big_int.cpp
    fast_evaluation.cpp
    interpolation.cpp
    modular_arithmetic.cpp
    ring_polynomial.cpp
)

target_link_libraries(tests PRIVATE src Catch2::Catch2WithMain)

add_test(NAME tests COMMAND tests)
