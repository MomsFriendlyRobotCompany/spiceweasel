# GTest -----------------
set(DOWNLOAD_EXTRACT_TIMESTAMP true)

include(FetchContent)

FetchContent_Declare(
    gtest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# list(APPEND libs gtest)
FetchContent_MakeAvailable( gtest )

enable_testing()

list(APPEND gtests
    main-gtest.cpp
)

# run all tests
add_executable(test-all ${gtests})
target_link_libraries(test-all GTest::gtest_main ${PROJECT_NAME})
include(GoogleTest)
gtest_discover_tests(test-all)