# Catch2
set(CATCH_INSTALL_DOCS OFF)
set(CATCH_INSTALL_EXTRAS OFF)
add_subdirectory(catch2)

set(KISSFFT_TEST OFF)
set(KISSFFT_TOOLS OFF)
set(KISSFFT_DATATYPE double)
set(KISSFFT_STATIC ON)
set(BUILD_SHARED_LIBS OFF)
add_subdirectory(kissfft)
