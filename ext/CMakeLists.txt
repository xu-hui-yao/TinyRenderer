target_include_directories(tiny-renderer PRIVATE stb)
target_include_directories(tiny-renderer PRIVATE imgui)
target_include_directories(tiny-renderer PRIVATE pugixml)
target_include_directories(tiny-renderer PRIVATE tinyexr)

add_library(pugixml STATIC pugixml/src/pugixml.cpp)
target_link_libraries(tiny-renderer PRIVATE pugixml)
