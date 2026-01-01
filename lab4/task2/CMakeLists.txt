add_executable(lab4task2 main.cpp functions.cpp main.h ${MATRIX_SRC} ../../lab1/task2/functions.cpp)

add_custom_target(plots_lab4task2
        COMMAND ${VENV_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/plots.py
        DEPENDS lab4task2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)