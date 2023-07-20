LIST(APPEND apps
    simple
    jacobian
    extended
)

foreach(APP ${apps})
    add_executable(${APP} ${APP}.cpp)
    # target_include_directories(${APP} PUBLIC ${CMAKE_SOURCE_DIR}/src)
    target_link_libraries(${APP} PUBLIC ${PROJECT_NAME})
endforeach()