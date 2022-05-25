file(GLOB ALGLIB_HEADERS "src/*.h")
file(GLOB ALGLIB_SOURCES "src/*.cpp")

file(COPY ${ALGLIB_HEADERS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/alglib)
file(COPY ${ALGLIB_SOURCES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/src)

if (WIN32)
    add_library(alglib STATIC ${ALGLIB_SOURCES} ${ALGLIB_HEADERS})
else()
    add_library(alglib SHARED ${ALGLIB_SOURCES} ${ALGLIB_HEADERS})
endif()

target_include_directories(alglib PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include/)
target_include_directories(alglib INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include/alglib)
