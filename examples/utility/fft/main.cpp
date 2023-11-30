#include <complex>
#include <cmath>
#include <vector>
#include <iostream>
#include "pocketfft_hdronly.h"
#include <matplot/matplot.h>
#include <set>

using namespace std;
using namespace pocketfft;

// floating point RNG which is good enough for simple demos
// Do not use for anything important!
inline double simple_drand()
{
    constexpr double norm = 1./RAND_MAX;
    return rand()*norm;
}

template<typename T> void crand(vector<complex<T>> &v)
{
    for (auto & i:v)
        i = complex<T>(simple_drand()-0.5, simple_drand()-0.5);
}

template<typename T> void crand(vector<T> &v)
{
    for (auto & i:v)
        i = simple_drand()-0.5;
}

template<typename T1, typename T2> long double l2err
        (const vector<T1> &v1, const vector<T2> &v2)
{
    long double sum1=0, sum2=0;
    for (size_t i=0; i<v1.size(); ++i)
    {
        long double dr = v1[i].real()-v2[i].real(),
                di = v1[i].imag()-v2[i].imag();
        long double t1 = sqrt(dr*dr+di*di), t2 = abs(v1[i]);
        sum1 += t1*t1;
        sum2 += t2*t2;
    }
    return sqrt(sum1/sum2);
}


int main()
{
    /*
    size_t len = 8192;

    shape_t shape{len};
    stride_t stridef(shape.size());
    stride_t strided(shape.size());
    stride_t stridel(shape.size());

    size_t tmpf = sizeof(complex<float>);
    size_t tmpd = sizeof(complex<double>);
    size_t tmpl=sizeof(complex<long double>);

    for (int i=shape.size()-1; i>=0; --i)
    {
        stridef[i]=tmpf;
        tmpf*=shape[i];
        strided[i]=tmpd;
        tmpd*=shape[i];
        stridel[i]=tmpl;
        tmpl*=shape[i];
    }
    size_t ndata=1;
    for (size_t i=0; i<shape.size(); ++i){
        ndata*=shape[i];
    }


    vector<complex<float>> dataf(ndata);
    vector<complex<double>> datad(ndata);
    vector<complex<long double>> datal(ndata);
    crand(dataf);
    for (size_t i=0; i<ndata; ++i)
    {
        datad[i] = dataf[i];
        datal[i] = dataf[i];
    }
    shape_t axes;
    for (size_t i=0; i<shape.size(); ++i)
        axes.push_back(i);
    auto resl = datal;
    auto resd = datad;
    auto resf = dataf;


    c2c(shape, stridel, stridel, axes, FORWARD,
        datal.data(), resl.data(), 1.L);
    c2c(shape, strided, strided, axes, FORWARD,
        datad.data(), resd.data(), 1.);
    c2c(shape, stridef, stridef, axes, FORWARD,
        dataf.data(), resf.data(), 1.f);

    //matplot::show();
//    c2c(shape, stridel, stridel, axes, POCKETFFT_BACKWARD,
//        resl.data(), resl.data(), 1.L/ndata);
    cout << l2err(resl, resf) << endl;
     */
    size_t len = 8192;

    shape_t shape{len};
    stride_t stridef(shape.size());
    stride_t strided(shape.size());
    stride_t stridel(shape.size());

    size_t tmpf = sizeof(float);
    size_t tmpd = sizeof(double);
    size_t tmpl = sizeof(long double);

    for (int i=shape.size()-1; i>=0; --i)
    {
        stridef[i]=tmpf;
        tmpf*=shape[i];
        strided[i]=tmpd;
        tmpd*=shape[i];
        stridel[i]=tmpl;
        tmpl*=shape[i];
    }

    size_t ndata=1;
    for (size_t i=0; i<shape.size(); ++i){
        ndata*=shape[i];
    }

    vector<float> dataf(ndata);
    vector<double> datad(ndata);
    vector<long double> datal(ndata);
    crand(dataf);
    for (size_t i=0; i<ndata; ++i)
    {
        datad[i] = dataf[i];
        datal[i] = dataf[i];
    }
    shape_t axes;
    for (size_t i=0; i<shape.size(); ++i){
        axes.push_back(i);
    }

    auto resl = datal;
    auto resd = datad;
    auto resf = dataf;

    //matplot::plot(datad, "k");

    //r2r_fftpack(shape, stridel, stridel, axes, true, FORWARD,
     //   datal.data(), resl.data(), 1.L, 1);

    r2r_fftpack(shape, strided, strided, axes, true, FORWARD,
        datad.data(), resd.data(), 1., 1);
    //r2r_fftpack(shape, stridef, stridef, axes, true, FORWARD,
        //dataf.data(), resf.data(), 1.f, 1);

    //matplot::plot(resd, "k");

    using namespace matplot;

    std::set<std::vector<double>> Y = {
            {16, 5, 9, 4}, {2, 11, 7, 14}, {3, 10, 6, 15}, {13, 8, 12, 1}};
    plot(Y);

    show();
    return 0;

}
