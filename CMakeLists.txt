cmake_minimum_required (VERSION 3.4)
project (pal-mem CXX)
   
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost COMPONENTS filesystem regex) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()

add_executable(pal-mem pal-mem.cpp) 
target_link_libraries(pal-mem ${Boost_LIBRARIES})
