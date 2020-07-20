if(NOT CPPCHECK_TOOL_PATH)
    return()
endif()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
configure_file(${PROJECT_SOURCE_DIR}/scripts/test/cppcheck.in.sh ${PROJECT_BINARY_DIR}/cppcheck.sh)

if(DEFINED ENV{NUM_THREADS})
    set(CPPCHECK_THREADS -j $ENV{NUM_THREADS})
endif()

add_custom_target(cppcheck
    COMMAND ${CPPCHECK_TOOL_PATH}
        --project=${PROJECT_BINARY_DIR}/compile_commands.json
        --language=c++
        --std=c++17
        --enable=all
        --inconclusive
        ${CPPCHECK_THREADS}
        -i ${PROJECT_BINARY_DIR}/CMakeFiles
        -i ${PROJECT_SOURCE_DIR}/ThirdParty
        -i ${PROJECT_SOURCE_DIR}/Applications/DataExplorer
        -i ${PROJECT_SOURCE_DIR}/Tests
        --xml --xml-version=2
        --output-file=${PROJECT_BINARY_DIR}/cppcheck.log
        ${PROJECT_SOURCE_DIR}
)
