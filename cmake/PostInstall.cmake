# PostInstall.cmake

include(cmake/TestInstall.cmake)

add_custom_target(post_install)

add_custom_target(
    install_all
    COMMAND ${CMAKE_COMMAND} --build . --target install
    COMMAND ${CMAKE_COMMAND} --build . --target post_install
    COMMENT "Executing `make install` and `make post install`."
    VERBATIM)