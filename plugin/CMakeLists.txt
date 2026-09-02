# CDP CLAP plugin collection: one binary, one s_plugins[] entry per effect.
# CLAP headers are vendored in ./clap (free-audio/clap v1.2.7, MIT).
add_library(cdp_fx MODULE src/cdp_gain.c src/cdp_blur.c src/cdp_distort.c)
target_include_directories(cdp_fx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/clap/include)

# Factory smoke-test: enumerates the built .clap and runs process() on each plugin.
add_executable(scan_check src/scan_check.c)
target_include_directories(scan_check PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/clap/include)
set_target_properties(cdp_fx PROPERTIES
    C_STANDARD 11
    PREFIX ""
    SUFFIX ".clap"
    OUTPUT_NAME "CDP FX Collection"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/NewRelease
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/NewRelease)
