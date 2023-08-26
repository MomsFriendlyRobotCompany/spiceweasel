LIST(APPEND apps
    simple
    jacobian
    extended
    nav
)

foreach(APP ${apps})
    add_executable(${APP} ${APP}.cpp)
    target_link_libraries(${APP} PRIVATE ${PROJECT_NAME})
endforeach()