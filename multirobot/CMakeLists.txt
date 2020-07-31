# create an executable
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEMPLATE_BINARY_DIR}/multirobot)
ADD_EXECUTABLE (simviz_demo simviz_demo.cpp ${TEMPLATE_COMMON_SOURCE})

# and link the library against the executable
TARGET_LINK_LIBRARIES (simviz_demo ${TEMPLATE_COMMON_LIBRARIES})

# export resources such as model files.
# NOTE: this requires an install build
SET(APP_RESOURCE_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources)
FILE(MAKE_DIRECTORY ${APP_RESOURCE_DIR})
FILE(COPY world.urdf DESTINATION ${APP_RESOURCE_DIR})
FILE(COPY panda_arm.urdf DESTINATION ${APP_RESOURCE_DIR})


