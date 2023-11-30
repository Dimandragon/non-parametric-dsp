add_rules("mode.debug", "mode.release")

add_requires("brew::matplotplusplus")

set_languages("c++23")
--set_toolset("cxx", "clang++")
--set_toolset("ld", "clang++")

target("pocketfft")

    set_kind("headeronly")
    --add_headerfiles("pocketfft/pocketfft_hdronly.h", {public = true})
    add_includedirs("pocketfft", {public = true})


package("matplotplusplus")
    add_deps("cmake")
    set_sourcedir(path.join(os.scriptdir(), "matplotplusplus"))
    on_install(function (package)
        local configs = {"-DMATPLOTPP_BUILD_EXAMPLES=OFF                  \
                             -DMATPLOTPP_BUILD_SHARED_LIBS=OFF                \
                             -DMATPLOTPP_BUILD_TESTS=OFF                     \
                             -DCMAKE_BUILD_TYPE=Release            \
                             -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON"}
        --table.insert(configs, "-DCMAKE_BUILD_TYPE=Release")
        --table.insert(configs, "-DMATPLOTPP_BUILD_EXAMPLES=OFF")
        --table.insert(configs, "-DMATPLOTPP_BUILD_SHARED_LIBS=ON")
        --table.insert(configs, "-DMATPLOTPP_BUILD_TESTS=OFF")
        --table.insert(configs, "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON")
        import("package.tools.cmake").install(package, configs)
    end)
    on_test(function (package)
        assert(package:has_cfuncs("add", {includes = "matplot/matplot.h"}))
    end)
package_end()

--add_requires("matplotplusplus")


target("non-parametric_dsp")
    --set_toolset("cxx", "clang")
    --set_toolset("ld", "clang++")
    set_kind("static")
    add_files("src/npdsp_concepts.ixx", "src/signals.ixx", "src/derivators.ixx"
           , "src/integrators.ixx", "src/filters.ixx", "src/inst_freq_computers.ixx",
           "src/mode_grabbers.ixx", "src/utility_math.ixx", "src/approximators.ixx")
    add_packages("brew::matplotplusplus")
    add_deps("pocketfft")


target("fft-example")
    --set_toolset("cxx", "clang")
    --set_toolset("ld", "clang++")
    set_kind("binary")
    add_files("examples/utility/fft/main.cpp")
    add_deps("pocketfft")
    add_packages("brew::matplotplusplus")
--
-- If you want to known more usage about xmake, please see https://xmake.io
--
-- ## FAQ
--
-- You can enter the project directory firstly before building project.
--
--   $ cd projectdir
--
-- 1. How to build project?
--
--   $ xmake
--
-- 2. How to configure project?
--
--   $ xmake f -p [macosx|linux|iphoneos ..] -a [x86_64|i386|arm64 ..] -m [debug|release]
--
-- 3. Where is the build output directory?
--
--   The default output directory is `./build` and you can configure the output directory.
--
--   $ xmake f -o outputdir
--   $ xmake
--
-- 4. How to run and debug target after building project?
--
--   $ xmake run [targetname]
--   $ xmake run -d [targetname]
--
-- 5. How to install target to the system directory or other output directory?
--
--   $ xmake install
--   $ xmake install -o installdir
--
-- 6. Add some frequently-used compilation flags in xmake.lua
--
-- @code
--    -- add debug and release modes
--    add_rules("mode.debug", "mode.release")
--
--    -- add macro definition
--    add_defines("NDEBUG", "_GNU_SOURCE=1")
--
--    -- set warning all as error
--    set_warnings("all", "error")
--
--    -- set language: c99, c++11
--    set_languages("c99", "c++11")
--
--    -- set optimization: none, faster, fastest, smallest
--    set_optimize("fastest")
--
--    -- add include search directories
--    add_includedirs("/usr/include", "/usr/local/include")
--
--    -- add link libraries and search directories
--    add_links("tbox")
--    add_linkdirs("/usr/local/lib", "/usr/lib")
--
--    -- add system link libraries
--    add_syslinks("z", "pthread")
--
--    -- add compilation and link flags
--    add_cxflags("-stdnolib", "-fno-strict-aliasing")
--    add_ldflags("-L/usr/local/lib", "-lpthread", {force = true})
--
-- @endcode
--

