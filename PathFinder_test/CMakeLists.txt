CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

SET(PROJ "PathFinder_test")

PROJECT(${PROJ})
SET(CMAKE_CXX_STANDARD 17)

IF (NOT DEFINED SMARTCGMS_COMMON_DIR)
	SET(SMARTCGMS_COMMON_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../common/" CACHE PATH "SmartCGMS 'common' directory location")
ENDIF()

INCLUDE_DIRECTORIES("${SMARTCGMS_COMMON_DIR}")

SET(PAGMO_DIR "" CACHE PATH "Pagmo2 directory")

FILE(GLOB src_files src/*.cpp src/*.h)
IF(WIN32)
	FILE(GLOB_RECURSE win_src_files src/win/*.cpp src/win/*.c src/win/*.def src/win/*.rc)
	SET(src_files ${src_files} ${win_src_files})
ENDIF()

SET(COMMON_FILES )

LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/rtl/scgmsLib.cpp")
LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/rtl/referencedImpl.cpp")
LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/rtl/FilterLib.cpp")
LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/rtl/DeviceLib.cpp")
LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/rtl/SolverLib.cpp")
LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/rtl/UILib.cpp")
LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/rtl/FilesystemLib.cpp")
LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/rtl/hresult.c")
LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/lang/dstrings.cpp")
LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/iface/DeviceIface.cpp")
LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/utils/string_utils.cpp")
LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/rtl/Dynamic_Library.cpp")
LIST(APPEND COMMON_FILES "${SMARTCGMS_COMMON_DIR}/scgms/utils/winapi_mapping.c")

ADD_EXECUTABLE(pathfinder_test ${src_files} ${COMMON_FILES})

set_target_properties(pathfinder_test
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/compiled/"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/compiled/"
)
