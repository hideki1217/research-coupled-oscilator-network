
add_executable(main main.cpp)
target_link_libraries(main boost mine thread_pool)

add_executable(order_checker order_checker.cpp)
target_link_libraries(order_checker boost mine)