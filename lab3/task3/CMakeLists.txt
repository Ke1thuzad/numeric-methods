add_executable(lab3task3 main.h main.cpp functions.cpp)

add_custom_target(plots_lab3task3
        COMMAND ${VENV_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/plots.py
        DEPENDS lab3task3
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

