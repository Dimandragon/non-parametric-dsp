#include <utility_math.hpp>
#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <complex>
#include <cstdlib>


float dotProduct(std::vector<float> a, std::vector<float> b) {
  float ans = 0;

  for (std::size_t i = 0; i < b.size(); i++) {
    ans += a[i] * b[i];
  }
  return ans;
}

std::vector<float> conv(std::vector<float> a, std::vector<float> b) {
  std::vector<float> cv = std::vector<float>();

  std::vector<float> n = std::vector<float>();
  for (std::size_t i = 0; i <= a.size() - b.size() - 1; i++) {
    n.clear();
    for (std::size_t j = 0; j < b.size(); j++) {
      n.push_back(a[i + j]);
    }
    cv.push_back(dotProduct(n, b));
  }
  return cv;
}



int main(){
    int a_len = 5000;
    int b_len = 400;

    std::vector<float> a = {};
    std::vector<float> b = {};

    for (int i = 0; i < a_len; i++){
        a.push_back(std::sin(static_cast<float>(i) / 63.69) + 1. + std::sin(static_cast<float>(i) / 500.));
    }
    for (int i = 0; i < b_len; i++){
        b.push_back(1);
    }

    auto average = 0.0;
    for (int i = 0; i < b_len; i++){
      average += b[i]; 
    }
    for (int i = 0; i < b_len; i++){
      b[i] /= average; 
    }
    
    for (int i = 0; i < b_len; i++){
      //b[i] -= 1.0/b_len; 
    }
    


    auto c_old = conv(a, b);

    std::vector<float> c = {};
    for(int i = 0; i < b.size()/2; i++){
        c.push_back(0.0);
    }

    for (int i = 0; i < c_old.size(); i++){
        c.push_back(c_old[i]);
    }

    for(int i = 0; i < b.size()/2; i++){
        c.push_back(0.0);
    }

    matplot::plot(b);
    matplot::show();

    matplot::plot(a);
    matplot::hold(1);
    matplot::plot(c);
    matplot::hold(0);
    matplot::show();



    NP_DSP::ONE_D::UTILITY_MATH::fastConvolution
      <std::vector<float>, std::vector<float>, std::vector<float>, float>(a, b, c);

    matplot::plot(a);
    matplot::hold(1);
    matplot::plot(c);
    matplot::hold(0);
    matplot::show();
    /*

    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT compute_buffer;

    for (int i = 0; i < 50; i++) {
        signal1.base->vec->push_back(std::rand());
        signal3.base->vec->push_back(0.0);
    }
    for (int i = 0; i < 10; i++) {
        signal2.base->vec->push_back(0.1);
    }

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::UTILITY_MATH::fastConvolution(signal1, signal2, signal1);

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);

    return 0;

    */
}