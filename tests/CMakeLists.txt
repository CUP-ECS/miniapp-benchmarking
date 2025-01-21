if (ENABLE_TESTS)

    # tstMesh needs debugging for >1 processes
    # blt_add_executable(NAME tstMesh
    #                    SOURCES tstMesh.cpp 
    #                    HEADERS tstMesh.hpp TestingBase.hpp 
    #                    DEPENDS_ON beatnik gtest
    #                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY})
    # blt_add_test(NAME MeshTests
    #              COMMAND tstMesh)
    
                 
    # blt_add_executable( NAME     tstBoundaryCondition
    #              SOURCES tstBoundaryCondition.cpp
    #              HEADERS tstBoundaryCondition.hpp TestingBase.hpp
    #              DEPENDS_ON beatnik gtest
    #              OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY})
    # blt_add_test(NAME BoundaryConditionTests
    #              COMMAND tstBoundaryCondition)

endif()
