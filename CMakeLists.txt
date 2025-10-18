CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

PROJECT("SmartCGMS")

message("========================")
message("Building dependencies...")
message("========================")

#########################################################################################
# ExcelFormat

execute_process(
        COMMAND cmake -S . -B ./build/ -G ${CMAKE_GENERATOR}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_DEBUG_POSTFIX=d # Without this debug build is unable to find the lib file
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/deps/ExcelFormat
)

message(STATUS "> ExcelFormat has been configured.")


execute_process(
        COMMAND cmake --build ./build/ --parallel
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/deps/ExcelFormat
)

message(STATUS "> ExcelFormat has been built.")

#########################################################################################
# XLNT

execute_process(
        COMMAND cmake -S . -B ./build/ -G ${CMAKE_GENERATOR}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_DEBUG_POSTFIX=d
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/deps/xlnt
)

message(STATUS "> xlnt has been configured.")

execute_process(
        COMMAND cmake --build ./build/ --parallel
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/deps/xlnt
)

message(STATUS "> xlnt has been built.")

#########################################################################################
# EIGEN 3

execute_process(
        COMMAND cmake -S . -B ./build/ -G ${CMAKE_GENERATOR}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/deps/Eigen3
)

message(STATUS "> Eigen3 has been configured.")

# Copy dependency DLLs to build dir (dll name is different depending on if it's a debug build or not)
if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/deps/xlnt/build/source/xlntd.dll ${CMAKE_CURRENT_BINARY_DIR}/compiled/xlntd.dll COPYONLY)
else ()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/deps/xlnt/build/source/xlnt.dll ${CMAKE_CURRENT_BINARY_DIR}/compiled/xlnt.dll COPYONLY)
endif ()

message("======================")
message("Configuring Project...")
message("======================")

ADD_SUBDIRECTORY("common")

ADD_SUBDIRECTORY("core")

IF (NOT DEFINED BUILD_CORE_ONLY OR NOT BUILD_CORE_ONLY)
    ADD_SUBDIRECTORY("console")
    ADD_SUBDIRECTORY("desktop")
ENDIF ()

ADD_SUBDIRECTORY("wrappers")

ADD_SUBDIRECTORY("PathFinder_test")

