add_subdirectory(boost)

# Catch2
option(CATCH_INSTALL_DOCS "" OFF)
option(CATCH_INSTALL_EXTRAS "" OFF)
add_subdirectory(catch2)

add_subdirectory(thread-pool)
