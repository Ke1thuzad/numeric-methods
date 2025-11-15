add_executable(lab3task1 main.cpp main.h
        functions.cpp)

add_custom_target(plots_lab3task1
        COMMAND ${VENV_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/plots.py
        DEPENDS lab3task1
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)