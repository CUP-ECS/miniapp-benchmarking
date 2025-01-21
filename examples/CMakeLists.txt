include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

blt_add_executable( NAME main
                SOURCES main.cpp 
                OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY})

# Tell CMAKE we need to link with these libraries. Our "FindMPI_Advance" script
# does not define an explicit target, but it does set the variable below.
# We also need to link with MPI, which is done through the target they provide
# (It's not C++ MPI, but rather the libraries for a C++ program using MPI)

# Finally, we also need to tell CMAKE where the library header is. 

target_link_libraries( main PRIVATE myproject )


install(TARGETS main DESTINATION ${CMAKE_INSTALL_BINDIR})
