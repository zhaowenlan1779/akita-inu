# Catch2
option(CATCH_INSTALL_DOCS "" OFF)
option(CATCH_INSTALL_EXTRAS "" OFF)
add_subdirectory(catch2)

option(KISSFFT_TEST "" OFF)
option(KISSFFT_TOOLS "" OFF)
set(KISSFFT_DATATYPE double)
option(KISSFFT_STATIC "" ON)
option(BUILD_SHARED_LIBS "" OFF)
add_subdirectory(kissfft)
