#note: this utility rebuild infrequently; change MY_DEPLOYMENT_EXE rather manually (when needs)
cmake_minimum_required(VERSION 3.13)
project( VaPrcFluTst )
set(CMAKE_GENERATOR_PLATFORM x64)
set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable( I16Prc )

target_sources( I16Prc PRIVATE 
   I16Prc.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/../I16_Cmn/I16_Intf.cpp
)

target_include_directories( I16Prc PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/.
  ${CMAKE_CURRENT_SOURCE_DIR}/../I16_Cmn
)

if (WIN32) # -- Windows
  target_compile_definitions( I16Prc PRIVATE
    _WIN32
    _CONSOLE
    $<$<CONFIG:Debug>:_DEBUG> )
  target_compile_options( I16Prc PRIVATE
    $<$<CONFIG:Debug>: /Od /Zi>
    $<$<CONFIG:Release>: /O2 /Oy-> )
  set( MY_DEPLOYMENT_EXE "e:/Img16Src/Tst/Debug64/Img16Prc.exe" )
else()  # -- Linux
  target_compile_definitions( I16Prc PRIVATE
    $<$<CONFIG:Debug>:_DEBUG> )
  target_compile_options( I16Prc PRIVATE
    -Wall
    $<$<CONFIG:Debug>:-O0 -g>
    $<$<CONFIG:Release>:-O3> )
  set( MY_DEPLOYMENT_EXE "${CMAKE_CURRENT_SOURCE_DIR}/../Tst/Release64/Img16Prc" )
endif()

add_custom_command(
  TARGET I16Prc POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:I16Prc> ${MY_DEPLOYMENT_EXE}
  COMMENT "Copying executable to default deployment folder..." )   
 