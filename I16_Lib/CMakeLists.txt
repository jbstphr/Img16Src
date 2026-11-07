cmake_minimum_required(VERSION 3.10)
#set(CMAKE_CUDA_ARCHITECTURES 61 75 80 86 87 89 90 100 110 120)
set(CMAKE_CUDA_ARCHITECTURES 75)
project( VAPrcFluDll LANGUAGES C CXX CUDA )
set(CMAKE_GENERATOR_PLATFORM x64)
set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

# don't do this under ~modern CMake #find_package( CUDAToolkit REQUIRED )

if (ENABLE_CUDA_DEBUG) # enable cuda-gdb (may significantly affect performance)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G")
else() # add line information to all builds for debug tools (exclusive to -G option)
  ##set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo")
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G")
endif()

add_library( I16Lib SHARED )

target_compile_features( I16Lib PRIVATE cxx_std_17 cuda_std_17 )
# if older nvcc for cuda architecture "61" 
add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:-Wno-deprecated-gpu-targets>)
add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:-D__NO_MATH_INLINES>)

#set_target_properties( I16Lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON )
target_sources( I16Lib PRIVATE 
  kernel.cu
  BS_Cmn.cpp
  BS_GIni.cpp
  BS_Ut0.cpp
  BS_Ut1.cpp
  BS_Ut2.cpp
  BS_Ut2a.cpp
  BS_Ut3.cpp
  BS_Ut4.cpp
  BS_Ut5.cpp
  BS_Ut5a.cpp
  BS_Ut5c.cpp
  BS_Ut5d.cpp
  BS_Ut6.cpp
  BS_Ut7.cpp
  BS_Ut7a.cpp
  BS_Ut7b.cpp
  BS_Ut7c.cpp
  I16Lib.cpp )

target_include_directories( I16Lib PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/.
  ${CMAKE_CURRENT_SOURCE_DIR}/../I16_Cmn )

if (WIN32)
  message( STATUS "CONFIG WIN32" )
  set( MY_DEPLOYMENT_EXE "D:/Img16Src/Tst/Img16Prc.dll" )
  include_directories( "/UTNVIDIA/CUDA132/include" )
  target_compile_definitions( I16Lib PRIVATE
    _WIN32
    _CONSOLE
    $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>
    $<$<CONFIG:Debug>:_DEBUG>)
  target_compile_options( I16Lib PRIVATE 
      $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Debug>>: -Xcompiler="/Od /Zi">
      $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Release>>: -Xcompiler="/O2 /Oy-">
      $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>: /Od /Zi>
      $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>: /O2 /Oy->
  )        
else()  #-Wextra 
  message( STATUS "CONFIG LINUX" )
  
  set( MY_DEPLOYMENT_EXE "${CMAKE_CURRENT_SOURCE_DIR}/../Tst/Img16Prc.so" )
  include_directories( "/usr/local/cuda/include" )
  target_compile_definitions( I16Lib PRIVATE
    $<$<CONFIG:Debug>:_DEBUG> )
  target_compile_options( I16Lib PRIVATE
    -Wall
    $<$<CONFIG:Debug>:-O0 -g>
    $<$<CONFIG:Release>:-O3> )
endif()

add_custom_command(
  TARGET I16Lib POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:I16Lib> ${MY_DEPLOYMENT_EXE}
  COMMENT "Copying executable to default deplyment folder..." )   
