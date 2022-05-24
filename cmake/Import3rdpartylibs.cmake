if(NOT BUILD_WITH_SIMBODY)
    FINDLIBHOME(SIMBODY_HOME)
    FIND_PACKAGE(Simbody_custom)
    IF(NOT Simbody_custom_FOUND)
        MESSAGE("Finding Simbody from system")
        FIND_PACKAGE(Simbody)
        IF(NOT Simbody_FOUND)
            MESSAGE(FATAL_ERROR "Simbody library not found")
        ENDIF()
        set(Simbody_DEBUG_LIBRARIES ${Simbody_LIBRARIES})
        set(Simbody_RELEASE_LIBRARIES ${Simbody_LIBRARIES})
    ENDIF()
    INCLUDE_DIRECTORIES("${Simbody_INCLUDE_DIR}")
    LINK_DIRECTORIES("${Simbody_LIB_DIR}")
    MESSAGE("Simbody found:")
    MESSAGE("  Simbody_INCLUDE_DIR:${Simbody_INCLUDE_DIR}")
    MESSAGE("  Simbody_LIB_DIR:${Simbody_LIB_DIR}")
    MESSAGE("  Simbody_LIBRARIES:${Simbody_LIBRARIES}")
    MESSAGE("  Simbody_RELEASE_LIBRARIES:${Simbody_RELEASE_LIBRARIES}")
    MESSAGE("  Simbody_DEBUG_LIBRARIES:${Simbody_DEBUG_LIBRARIES}")
endif()

if(NOT EMSCRIPTEN)
    FINDLIBHOME(TBB_HOME)
    FIND_PACKAGE(TBB_custom)
    IF(NOT TBB_custom_FOUND)
        MESSAGE("Finding TBB from system")
        FIND_PACKAGE(TBB)
        IF(NOT TBB_FOUND)
            MESSAGE(FATAL_ERROR "TBB library not found")
        ENDIF()
        set(TBB_LIBRARYS TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy)
    ENDIF()
    INCLUDE_DIRECTORIES("${TBB_INCLUDE_DIR}")
    LINK_DIRECTORIES("${TBB_LIB_DIR}")
    MESSAGE("TBB found:")
    MESSAGE("  TBB_INCLUDE_DIR:${TBB_INCLUDE_DIR}")
    MESSAGE("  TBB_LIB_DIR:${TBB_LIB_DIR}")
    MESSAGE("  TBB_LIBRARYS:${TBB_LIBRARYS}")

    if(DEFINED BOOST_AVAILABLE)
        FINDLIBHOME(BOOST_HOME)
        set(BOOST_ROOT "${BOOST_HOME}")
        IF(MSVC)
        FIND_PACKAGE(Boost REQUIRED 
            COMPONENTS program_options)
        ELSE(MSVC)
            FIND_PACKAGE(Boost REQUIRED 
            COMPONENTS program_options filesystem system )
        ENDIF(MSVC)
        IF(Boost_FOUND)
            INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}")
            LINK_DIRECTORIES("${Boost_LIBRARY_DIRS}")
            MESSAGE("${Boost_INCLUDE_DIRS}")
            MESSAGE("${Boost_LIBRARY_DIRS}")
            MESSAGE("${Boost_LIBRARIES}")
        ELSE(Boost_FOUND)
            MESSAGE(FATAL_ERROR "Boost library not found")
        ENDIF(Boost_FOUND)
    endif()
endif()

IF(MSVC)
    FINDLIBHOME(GTEST_HOME)
    set(GTest_ROOT "${GTEST_HOME}") 
    find_package(GTest CONFIG REQUIRED)
    IF(GTest_FOUND)
       INCLUDE_DIRECTORIES(${GTEST_HOME}/include)
       LINK_DIRECTORIES("${GTEST_HOME}/lib")
    ELSE(GTest_FOUND)
        MESSAGE(FATAL_ERROR "GTest library not found")
    ENDIF(GTest_FOUND)
ENDIF(MSVC)
