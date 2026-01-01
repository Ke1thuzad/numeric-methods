add_executable(lab4task1 main.h main.cpp functions.cpp)

add_custom_target(plots_lab4task1
        COMMAND ${VENV_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/plots.py
        DEPENDS lab4task1
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)