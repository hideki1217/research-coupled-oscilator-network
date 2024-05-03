
add_executable(main main.cpp)
target_link_libraries(main boost mine thread_pool)
target_compile_definitions(main PRIVATE MINE_UNIT_N=${MINE_UNIT_N})

add_executable(order_checker order_checker.cpp)
target_link_libraries(order_checker boost mine)
target_compile_definitions(order_checker PRIVATE MINE_UNIT_N=${MINE_UNIT_N})

add_executable(test_order_evaluator test_order_evaluator.cpp)
target_link_libraries(test_order_evaluator boost mine)
target_compile_definitions(test_order_evaluator PRIVATE MINE_UNIT_N=2)
add_test(
    NAME test_order_evaluator
    COMMAND test_order_evaluator
)