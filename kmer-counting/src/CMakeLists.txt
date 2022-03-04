add_library(PolyassemblyLib SHARED 
				#polyassembly_countkmers.cpp
				polyassembly_findkmers.cpp
				argumentparser.cpp
				jellyfishkmercounter.cpp
				uniquekmers.cpp
								)
				
				
#add_executable(polyassembly_countkmers polyassembly_countkmers.cpp)
add_executable(polyassembly_findkmers polyassembly_findkmers.cpp)

#target_link_libraries(polyassembly_countkmers PolyassemblyLib ${JELLYFISH_LDFLAGS_OTHER})
#target_link_libraries(polyassembly_countkmers PolyassemblyLib ${JELLYFISH_LIBRARIES})

target_link_libraries(polyassembly_findkmers PolyassemblyLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(polyassembly_findkmers PolyassemblyLib ${JELLYFISH_LIBRARIES})
