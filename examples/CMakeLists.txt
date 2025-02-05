blt_add_executable( NAME main
        SOURCES main.cpp
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY})


#TCLAP_start
find_path(TCLAP_INCLUDE_DIR tclap/CmdLine.h
        PATHS
        ${SPACK_PREFIX}/include
)

# Check if we found TCLAP
if (TCLAP_INCLUDE_DIR)
    message(STATUS "Found TCLAP: ${TCLAP_INCLUDE_DIR}")
    set(TCLAP_FOUND TRUE)
else ()
    message(FATAL_ERROR "TCLAP not found!")
endif ()
set(TCLAP_INCLUDE_DIRS ${TCLAP_INCLUDE_DIR})


include_directories(${TCLAP_INCLUDE_DIR})
#TCLAP_end



# nlohmann_json_start
find_path(NLOHMANN_JSON_INCLUDE_DIR nlohmann/json.hpp
        PATHS
        ${SPACK_PREFIX}/include
)

# Check if we found nlohmann/json
if (NLOHMANN_JSON_INCLUDE_DIR)
    message(STATUS "Found nlohmann/json: ${NLOHMANN_JSON_INCLUDE_DIR}")
    set(NLOHMANN_JSON_FOUND TRUE)
else ()
    message(FATAL_ERROR "nlohmann/json not found!")
endif ()
set(NLOHMANN_JSON_INCLUDE_DIRS ${NLOHMANN_JSON_INCLUDE_DIR})

include_directories(${NLOHMANN_JSON_INCLUDE_DIR})
# nlohmann_json_end





target_link_libraries( main PRIVATE myproject )


install(TARGETS main DESTINATION ${CMAKE_INSTALL_BINDIR})
