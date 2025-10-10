cmake_minimum_required(VERSION 3.8)
project(csf_pcl)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(PCL 1.3 REQUIRED COMPONENTS filters)

add_library(${PROJECT_NAME} SHARED
  src/cloth_simulation_filter.cpp
)
target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)
target_link_libraries(${PROJECT_NAME}
  ${PCL_LIBRARIES}
)

install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}_export
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)
install(
  DIRECTORY include/
  DESTINATION include/${PROJECT_NAME}
)
install(
  EXPORT ${PROJECT_NAME}_export
  FILE ${PROJECT_NAME}-config.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION lib/cmake/${PROJECT_NAME}
)
