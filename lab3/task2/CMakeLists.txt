add_executable(lab3task2 main.cpp functions.cpp main.h ../../lab1/task2/functions.cpp)

add_custom_target(plots_lab3task2
        COMMAND ${VENV_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/plots.py
        DEPENDS lab3task2
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)