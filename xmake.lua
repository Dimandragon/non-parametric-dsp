add_rules("mode.debug", "mode.release")

set_languages("c++23")


--xmake f --cxx=clang++ --cc=clang
--xmake project -k compile_commands

target("icecream")
    set_kind("headeronly")
    --add_includedirs("$(projectdir)/icecream-cpp", {public = true})
    add_includedirs("$(projectdir)/my_icecream", {public = true})

task("build matplot++")
    on_run(function()
        fl = os.exec("sh matplot_build.sh")
    end)

target("matplot++_external")
    set_kind("headeronly")
    --on_load(function(target)
        --import("core.base.task")
        --task.run("build matplot++")
    --end)
    add_headerfiles("$(projectdir)/matplotplusplus/source/matplot/matplot.h", {public = true})
    add_headerfiles("$(projectdir)/matplotplusplus/source/3rd_party/nodesoup/include/nodesoup.h", {public = true})
    --add_headerfiles("$(projectdir)/matplotplusplus/build/local/_deps/glad-build/include/glad/glad.h", {public = true})
    --add_headerfiles("$(projectdir)/matplotplusplus/build/local/_deps/glad-build/include/KHR/khrplatform.h", {public = true})
    

    add_links("$(projectdir)/matplotplusplus/build/local/source/matplot/libmatplot.a", {public = true})
    --add_links("$(projectdir)/matplotplusplus/build/local/source/matplot/libmatplot_opengl.a", {public = true})
    add_links("$(projectdir)/matplotplusplus/build/local/source/3rd_party/libnodesoup.a", {public = true})
    --add_links("$(projectdir)/matplotplusplus/build/local/_deps/glad-build/libglad.a", {public = true})


target("pocketfft")
    set_kind("headeronly")
    --add_headerfiles("pocketfft/pocketfft_hdronly.h", {public = true})
    add_includedirs("pocketfft", {public = true})


target("non-parametric_dsp")
    set_kind("static")
    add_files("src/npdsp_concepts.ixx", "src/signals.ixx", "src/derivators.ixx"
           , "src/integrators.ixx", "src/filters.ixx", "src/inst_freq_computers.ixx",
           "src/mode_grabbers.ixx", "src/utility_math.ixx", "src/approximators.ixx", "src/config.ixx")
    add_deps("pocketfft")
    add_deps("matplot++_external")
    add_deps("icecream")


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

target("generic_signal_example")
    set_kind("binary")
    add_files("examples/signals/generic_signal.cpp")
    add_deps("matplot++_external")
    add_deps("non-parametric_dsp")
    add_deps("icecream")

target("riman_integral_example")
    set_kind("binary")
    add_files("examples/integrators/riman.cpp")
    add_deps("matplot++_external")
    add_deps("non-parametric_dsp")
    add_deps("icecream")

target("finnite_difference_example")
    set_kind("binary")
    add_files("examples/derivators/finnite_difference.cpp")
    add_deps("matplot++_external")
    add_deps("non-parametric_dsp")
    add_deps("icecream")

target("fs_based_approx_example")
    set_kind("binary")
    add_files("examples/approximators/fourier_series_based.cpp")
    add_deps("matplot++_external")
    add_deps("non-parametric_dsp")
    add_deps("icecream")

--target("opengl_plotting_test")
    --set_kind("binary")
    --add_files("examples/plotting/opengl_test.cpp")
    --add_deps("matplot++_external")