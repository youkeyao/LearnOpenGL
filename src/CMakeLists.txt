file( GLOB SOURCES main.cpp )

include_directories( . )

add_executable( main ${SOURCES} )
target_link_libraries( main
    PRIVATE
        glfw
        glad
        assimp
        glm::glm
        stb_image
)