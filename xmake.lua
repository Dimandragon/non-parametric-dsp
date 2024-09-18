add_rules("mode.debug", "mode.release")

set_languages("c++23")


--xmake f --cxx=clang++ --cc=clang -m debug --debugger=lldb-16
--xmake project -k compile_commands

--add_requires("libomp")
--add_requires("libtorch")

target("boost_interpolate")
    set_kind("headeronly")
    add_includedirs("interpolate_from_boost", {public = true})

target("pocketfft")
    set_kind("headeronly")
    --add_headerfiles("pocketfft/pocketfft_hdronly.h", {public = true})
    add_includedirs("pocketfft", {public = true})


target("non-parametric_dsp")
    set_kind("headeronly")
    --add_headerfiles("src/npdsp_concepts.hpp", "src/signals.hpp", "src/derivators.hpp",
    --        "src/integrators.hpp", "src/filters.hpp", "src/inst_freq_computers.hpp",
    --        "src/utility_math.hpp", "src/approximators.hpp", "src/config.hpp", 
    --        "src/phase_computers.hpp", "src/inst_ampl_computers.hpp"
    --        ,"src/modes_extractors.hpp", {public = true}
    --        )
    add_includedirs("np_dsp", {public = true})

    add_deps("pocketfft")
    add_deps("boost_interpolate")



target("main_extractor_example")
    set_kind("binary")
    add_files("example/main.cpp")
    add_deps("non-parametric_dsp")
    --add_deps("my-icecream")
