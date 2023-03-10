add_library(GenomicPIR globals.hpp client.hpp client.cpp server.hpp server.cpp comparator.cpp comparator.hpp tools.cpp tools.hpp)
target_link_libraries(GenomicPIR helib)

add_executable(main main.cpp)

target_link_libraries(main GenomicPIR)

install(TARGETS GenomicPIR
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
)
