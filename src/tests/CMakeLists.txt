add_executable(tests
    big_int.cpp
    modular_arithmetic.cpp
    ring_polynomial.cpp
)

target_link_libraries(tests PRIVATE src Catch2::Catch2WithMain)

add_test(NAME tests COMMAND tests)
