add_rules("mode.debug", "mode.release")

set_languages("c++23")

--xmake f --cxx=clang++ --cc=clang

task("build matplot++")
    on_run(function()
        fl = os.exec("sh matplot_build.sh")
    end)

target("matplot++_external")
    set_kind("headeronly")
    on_load(function(target)
        import("core.base.task")
        task.run("build matplot++")
    end)
    add_headerfiles("$(projectdir)/matplotplusplus/source/matplot/matplot.h", {public = true})
    add_headerfiles("$(projectdir)/matplotplusplus/source/3rd_party/nodesoup/include/nodesoup.h", {public = true})

    add_links("$(projectdir)/matplotplusplus/build/local/source/matplot/libmatplot.a", {public = true})
    add_links("$(projectdir)/matplotplusplus/build/local/source/3rd_party/libnodesoup.a", {public = true})


target("pocketfft")
    set_kind("headeronly")
    --add_headerfiles("pocketfft/pocketfft_hdronly.h", {public = true})
    add_includedirs("pocketfft", {public = true})


target("non-parametric_dsp")
    set_kind("static")
    add_files("src/npdsp_concepts.ixx", "src/signals.ixx", "src/derivators.ixx"
           , "src/integrators.ixx", "src/filters.ixx", "src/inst_freq_computers.ixx",
           "src/mode_grabbers.ixx", "src/utility_math.ixx", "src/approximators.ixx")
    add_deps("pocketfft")
    add_deps("matplot++_external")


target("fft-example")
    set_kind("binary")
    add_files("examples/utility/fft/main.cpp")
    add_deps("non-parametric_dsp")
    add_deps("pocketfft")
    add_deps("matplot++_external")


target("fft-experiments")
    set_kind("binary")
    add_files("examples/utility/fft/experiments.cpp")
    add_deps("non-parametric_dsp")
    add_deps("pocketfft")
    add_deps("matplot++_external")

target("theta_loss_check")
    set_kind("binary")
    add_files("examples/utility/fft/cos_theta_loss.cpp")
    add_deps("matplot++_external")

target("ampl_loss_check")
    set_kind("binary")
    add_files("examples/utility/fft/cos_ampl_loss.cpp")
    add_deps("matplot++_external")
--[[
target("matplot++")

    set_kind("static")
    add_files("matplotplusplus/source/matplot/axes_objects/**.cpp")
    add_files("matplotplusplus/source/matplot/backend/**.cpp")
    add_files("matplotplusplus/source/matplot/core/**.cpp")
    add_files("matplotplusplus/source/matplot/detail/**.cpp")
    add_files("matplotplusplus/source/matplot/freestanding/**.cpp")
    add_files("matplotplusplus/source/matplot/util/**.cpp")
    add_files("matplotplusplus/source/3rd_party/nodesoup/src/**.cpp")
    add_includedirs("matplotplusplus/source/matplot/axes_objects")
    add_includedirs("matplotplusplus/source/matplot/backend")
    add_includedirs("matplotplusplus/source/matplot/core")
    add_includedirs("matplotplusplus/source/matplot/detail")
    add_includedirs("matplotplusplus/source/matplot/freestanding")
    add_includedirs("matplotplusplus/source/matplot/util")
    add_includedirs("matplotplusplus/source/3rd_party/nodesoup/src")
    add_includedirs("matplotplusplus/source/3rd_party/nodesoup/include")
    add_includedirs("matplotplusplus/source/3rd_party/simg")
    add_includedirs("matplotplusplus/source/matplot", {public = true})
]]--

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