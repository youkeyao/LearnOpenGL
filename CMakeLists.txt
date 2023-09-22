cmake_minimum_required( VERSION 3.16 )
set( CMAKE_CXX_STANDARD 17 )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/bin )

project( OpenGL )

set( GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE )
set( GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE )
set( GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE )

set( BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE )
set( ASSIMP_BUILD_TESTS OFF CACHE BOOL "" FORCE )
set( ASSIMP_INSTALL OFF CACHE BOOL "" FORCE )
set( ASSIMP_BUILD_ASSIMP_VIEW OFF CACHE BOOL "" FORCE )

add_subdirectory( vendor/glfw )
add_subdirectory( vendor/glad )
add_subdirectory( vendor/assimp )
add_subdirectory( vendor/glm )
add_subdirectory( vendor/stb_image )

include_directories(
    PUBLIC
        src/
        vendor/glfw/include
        vendor/glad/include
        vendor/assimp/include
        vendor/stb_image
)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.cpp)
add_executable( main ${SOURCES} )
target_link_libraries( main
    PRIVATE
        glfw
        glad
        assimp
        glm::glm
        stb_image
)