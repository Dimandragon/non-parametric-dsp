cd matplotplusplus

export CC=/usr/bin/clang
export CCX=/usr/bin/clang++
cmake -B build/local \
            -DMATPLOTPP_BUILD_EXAMPLES=OFF \
            -DMATPLOTPP_BUILD_SHARED_LIBS=ON \
            -DMATPLOTPP_BUILD_TESTS=OFF \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \

cmake --build build/local

cd ../