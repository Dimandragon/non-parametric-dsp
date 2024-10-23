// Copyright Nick Thompson, 2020
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_CUBIC_HERMITE_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_CUBIC_HERMITE_DETAIL_HPP
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>
#include <cmath>

namespace boost {
namespace math {
namespace interpolators {
namespace detail {

template<class RandomAccessContainer>
class cubic_hermite_detail {
public:
    using Real = typename RandomAccessContainer::value_type;
    using Size = typename RandomAccessContainer::size_type;

    cubic_hermite_detail(RandomAccessContainer && x, RandomAccessContainer && y, RandomAccessContainer dydx)
     : x_{std::move(x)}, y_{std::move(y)}, dydx_{std::move(dydx)}
    {
        using std::abs;
        using std::isnan;
        if (x_.size() != y_.size())
        {
            throw std::domain_error("There must be the same number of ordinates as abscissas.");
        }
        if (x_.size() != dydx_.size())
        {
            throw std::domain_error("There must be the same number of ordinates as derivative values.");
        }
        if (x_.size() < 2)
        {
            throw std::domain_error("Must be at least two data points.");
        }
        Real x0 = x_[0];
        for (size_t i = 1; i < x_.size(); ++i)
        {
            Real x1 = x_[i];
            if (x1 <= x0)
            {
                std::ostringstream oss;
                oss.precision(std::numeric_limits<Real>::digits10+3);
                oss << "Abscissas must be listed in strictly increasing order x0 < x1 < ... < x_{n-1}, ";
                oss << "but at x[" << i - 1 << "] = " << x0 << ", and x[" << i << "] = " << x1 << ".\n";
                throw std::domain_error(oss.str());
            }
            x0 = x1;
        }
    }

    void push_back(Real x, Real y, Real dydx)
    {
        using std::abs;
        using std::isnan;
        if (x <= x_.back())
        {
             throw std::domain_error("Calling push_back must preserve the monotonicity of the x's");
        }
        x_.push_back(x);
        y_.push_back(y);
        dydx_.push_back(dydx);
    }

    Real operator()(Real x) const
    {
        if  (x < x_[0] || x > x_.back())
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x_[0] << ", " << x_.back() << "]";
            throw std::domain_error(oss.str());
        }
        // We need t := (x-x_k)/(x_{k+1}-x_k) \in [0,1) for this to work.
        // Sadly this neccessitates this loathesome check, otherwise we get t = 1 at x = xf.
        if (x == x_.back())
        {
            return y_.back();
        }

        auto it = std::upper_bound(x_.begin(), x_.end(), x); 
        //it - итератор на первый элемент x_, который больше x
        auto i = std::distance(x_.begin(), it) -1;
        //i - номер элемента под итератором it
        Real x0 = *(it-1); //элемент x_ до it 
        Real x1 = *it; //элемент x_ под it
        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real s0 = dydx_[i];
        Real s1 = dydx_[i+1];
        Real dx = (x1-x0);
        Real t = (x-x0)/dx;

        // See the section 'Representations' in the page
        // https://en.wikipedia.org/wiki/Cubic_Hermite_spline
        Real y = (1-t)*(1-t)*(y0*(1+2*t) + s0*(x-x0))
              + t*t*(y1*(3-2*t) + dx*s1*(t-1));
        /*Real a_ = (s0/pow(dx, 2) + s1/pow(dx, 2) + 2*y0/pow(dx, 3) 
            - 2*y1/pow(dx, 3)) ; // * pow(x, 3)
        Real b_ = (-2*s0/dx - s1/dx - 3*s0*x0/pow(dx, 2) - 3*s1*x0/pow(dx, 2) 
            - 3*y0/pow(dx, 2) + 3*y1/pow(dx, 2) - 6*x0*y0/pow(dx, 3) 
            + 6*x0*y1/pow(dx, 3)); // * pow(x, 2)
        Real c_ = (s0 + 4*s0*x0/dx + 2*s1*x0/dx + 3*s0*pow(x0, 2)/pow(dx, 2) 
            + 3*s1*pow(x0, 2)/pow(dx, 2) + 6*x0*y0/pow(dx, 2) 
            - 6*x0*y1/pow(dx, 2) + 6*pow(x0, 2)*y0/pow(dx, 3) 
            - 6*pow(x0, 2)*y1/pow(dx, 3)); // * x
        Real d_ =  y0 - 2*s0*pow(x0, 2)/dx - s1*pow(x0, 2)/dx - s0*pow(x0, 3)/pow(dx, 2) 
            - s1*pow(x0, 3)/pow(dx, 2) - 3*pow(x0, 2)*y0/pow(dx, 2) + 3*pow(x0, 2)*y1/pow(dx, 2) 
            - 2*pow(x0, 3)*y0/pow(dx, 3) + 2*pow(x0, 3)*y1/pow(dx, 3) - s0*x0; 
        double y_ = a_ * pow(x, 3) + b_ * pow(x, 2) + c_ * x + d_;
        IC(y, y_);*/
        return y;
    }

    Real derive(Real x, double a) const
    {
        if  (x < x_[0] || x > x_.back())
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x_[0] << ", " << x_.back() << "]";
            throw std::domain_error(oss.str());
        }
        if (x == x_.back())
        {
            return y_.back();
        } //что делать с последней точкой? нужна лт для нее дробная производная?

        //double tgamma(double a); //надо нормально подключить функцию гамма
        auto it = std::upper_bound(x_.begin(), x_.end(), x);
        auto i = std::distance(x_.begin(), it) -1;
        Real x0 = *(it-1);
        Real x1 = *it;
        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real s0 = dydx_[i];
        Real s1 = dydx_[i+1];
        Real dx = (x1-x0);
        Real t = (x-x0)/dx;


        Real a_ = (s0/pow(dx, 2) + s1/pow(dx, 2) + 2*y0/pow(dx, 3) 
            - 2*y1/pow(dx, 3)) ; // * pow(x, 3)
        Real b_ = (-2*s0/dx - s1/dx - 3*s0*x0/pow(dx, 2) - 3*s1*x0/pow(dx, 2) 
            - 3*y0/pow(dx, 2) + 3*y1/pow(dx, 2) - 6*x0*y0/pow(dx, 3) 
            + 6*x0*y1/pow(dx, 3)); // * pow(x, 2)
        Real c_ = (s0 + 4*s0*x0/dx + 2*s1*x0/dx + 3*s0*pow(x0, 2)/pow(dx, 2) 
            + 3*s1*pow(x0, 2)/pow(dx, 2) + 6*x0*y0/pow(dx, 2) 
            - 6*x0*y1/pow(dx, 2) + 6*pow(x0, 2)*y0/pow(dx, 3) 
            - 6*pow(x0, 2)*y1/pow(dx, 3)); // * x
        Real d_ =  y0 - 2*s0*pow(x0, 2)/dx - s1*pow(x0, 2)/dx - s0*pow(x0, 3)/pow(dx, 2) 
            - s1*pow(x0, 3)/pow(dx, 2) - 3*pow(x0, 2)*y0/pow(dx, 2) + 3*pow(x0, 2)*y1/pow(dx, 2) 
            - 2*pow(x0, 3)*y0/pow(dx, 3) + 2*pow(x0, 3)*y1/pow(dx, 3) - s0*x0; 

        //IC(a_, pow(x, 3.0-a), tgamma(4.0-a), b_, pow(x, 2.0-a), tgamma(3.0-a), c_, pow(x, 1.0-a), tgamma(2.0-a));
        /*return a_ * pow(x, 3-a)/tgamma(4-a) + b_ * pow(x, 2-a) / tgamma(3-a)
            + c_ * pow(x, 1-a) / tgamma(2-a);*/
        //IC(x*x, pow(x, 2));
        /*Real derive = a_ * pow(x, 3.0-a) * tgamma(4) / tgamma(4.0 - a) * 
            + b_ * pow(x, 2.0-a) * tgamma(3) / tgamma(3.0 - a)
            + c_ * pow(x, 1.0 - a) * tgamma(2) / tgamma(2.0 - a);
            + d_ * pow(x, 0.0 - a) * tgamma(1) / tgamma(1.0 - a);*/
        
        /*Real derive = 6*a_*pow(x, 3.0)/(pow(x, a)*tgamma(4.0 - a)) + 
            2*b_*pow(x, 2.0)/(pow(x, a)*tgamma(3.0 - a)) + 
            c_*x/(pow(x, a)*tgamma(2.0 - a)) + 
            d_/(pow(x, a)*tgamma(1 - a));*/
        //return a_ * x*x * 3. + b_ * 2. * x + c_;
        Real tc_ = dx*s0; //*t + 
        Real ta_ = (dx*s0 + dx*s1 + 2*y0 - 2*y1); // *t**3 
        Real tb_ = (-2*dx*s0 - dx*s1 - 3*y0 + 3*y1); // * t**2 
        Real td_ = y0;
        Real derive_t = 
            ta_ * pow(t, 3.0-a) * tgamma(4) / tgamma(4.0 - a) * 
            + tb_ * pow(t, 2.0-a) * tgamma(3) / tgamma(3.0 - a)
            + tc_ * pow(t, 1.0 - a) * tgamma(2) / tgamma(2.0 - a);
            //+ td_ * pow(t, 0.0 - a) * tgamma(1) / tgamma(1.0 - a);
        Real tx_c = dx; //*x
        Real tx_d = -x0 * dx;

        Real derive_t_x = tx_c * pow(x, 1.0 - a) * tgamma(2) / tgamma(2.0 - a);
            //+ tx_d * pow(x, 0.0 - a) * tgamma(1) / tgamma(1.0 - a);
        //Real derive_x = pow(x, 1.0 - a) * tgamma(2) / tgamma(2.0 - a);
        //return derive * derive_x;
        //Real dtdx = 
        return derive_t;// * derive_t_x;
    }

    Real prime(Real x) const
    {
        if  (x < x_[0] || x > x_.back())
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x_[0] << ", " << x_.back() << "]";
            throw std::domain_error(oss.str());
        }
        if (x == x_.back())
        {
            return dydx_.back();
        }
        auto it = std::upper_bound(x_.begin(), x_.end(), x);
        auto i = std::distance(x_.begin(), it) -1;
        Real x0 = *(it-1);
        Real x1 = *it;
        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real s0 = dydx_[i];
        Real s1 = dydx_[i+1];
        Real dx = (x1-x0);

        Real d1 = (y1 - y0 - s0*dx)/(dx*dx);
        Real d2 = (s1 - s0)/(2*dx);
        Real c2 = 3*d1 - 2*d2;
        Real c3 = 2*(d2 - d1)/dx;
        return s0 + 2*c2*(x-x0) + 3*c3*(x-x0)*(x-x0); 
    }


    friend std::ostream& operator<<(std::ostream & os, const cubic_hermite_detail & m)
    {
        os << "(x,y,y') = {";
        for (size_t i = 0; i < m.x_.size() - 1; ++i)
        {
            os << "(" << m.x_[i] << ", " << m.y_[i] << ", " << m.dydx_[i] << "),  ";
        }
        auto n = m.x_.size()-1;
        os << "(" << m.x_[n] << ", " << m.y_[n] << ", " << m.dydx_[n] << ")}";
        return os;
    }

    Size size() const
    {
        return x_.size();
    }

    int64_t bytes() const
    {
        return 3*x_.size()*sizeof(Real) + 3*sizeof(x_);
    }

    std::pair<Real, Real> domain() const
    {
        return {x_.front(), x_.back()};
    }

    RandomAccessContainer x_;
    RandomAccessContainer y_;
    RandomAccessContainer dydx_;
};

template<class RandomAccessContainer>
class cardinal_cubic_hermite_detail {
public:
    using Real = typename RandomAccessContainer::value_type;
    using Size = typename RandomAccessContainer::size_type;

    cardinal_cubic_hermite_detail(RandomAccessContainer && y, RandomAccessContainer dydx, Real x0, Real dx)
    : y_{std::move(y)}, dy_{std::move(dydx)}, x0_{x0}, inv_dx_{1/dx}
    {
        using std::abs;
        using std::isnan;
        if (y_.size() != dy_.size())
        {
            throw std::domain_error("There must be the same number of derivatives as ordinates.");
        }
        if (y_.size() < 2)
        {
            throw std::domain_error("Must be at least two data points.");
        }
        if (dx <= 0)
        {
            throw std::domain_error("dx > 0 is required.");
        }

        for (auto & dy : dy_)
        {
            dy *= dx;
        }
    }

    // Why not implement push_back? It's awkward: If the buffer is circular, x0_ += dx_.
    // If the buffer is not circular, x0_ is unchanged.
    // We need a concept for circular_buffer!

    inline Real operator()(Real x) const
    {
        const Real xf = x0_ + (y_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf)
        {
            return y_.back();
        }
        return this->unchecked_evaluation(x);
    }

    inline Real unchecked_evaluation(Real x) const
    {
        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(y_.size())>(ii);
        Real t = s - ii;
        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real dy0 = dy_[i];
        Real dy1 = dy_[i+1];

        Real r = 1-t;
        return r*r*(y0*(1+2*t) + dy0*t)
              + t*t*(y1*(3-2*t) - dy1*r);
    }

    inline Real prime(Real x) const
    {
        const Real xf = x0_ + (y_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf)
        {
            return dy_.back()*inv_dx_;
        }
        return this->unchecked_prime(x);
    }

    inline Real unchecked_prime(Real x) const
    {
        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(y_.size())>(ii);
        Real t = s - ii;
        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real dy0 = dy_[i];
        Real dy1 = dy_[i+1];

        Real dy = 6*t*(1-t)*(y1 - y0)  + (3*t*t - 4*t+1)*dy0 + t*(3*t-2)*dy1;
        return dy*inv_dx_;
    }


    Size size() const
    {
        return y_.size();
    }

    int64_t bytes() const
    {
        return 2*y_.size()*sizeof(Real) + 2*sizeof(y_) + 2*sizeof(Real);
    }

    std::pair<Real, Real> domain() const
    {
        Real xf = x0_ + (y_.size()-1)/inv_dx_;
        return {x0_, xf};
    }

private:

    RandomAccessContainer y_;
    RandomAccessContainer dy_;
    Real x0_;
    Real inv_dx_;
};


template<class RandomAccessContainer>
class cardinal_cubic_hermite_detail_aos {
public:
    using Point = typename RandomAccessContainer::value_type;
    using Real = typename Point::value_type;
    using Size = typename RandomAccessContainer::size_type;

    cardinal_cubic_hermite_detail_aos(RandomAccessContainer && dat, Real x0, Real dx)
    : dat_{std::move(dat)}, x0_{x0}, inv_dx_{1/dx}
    {
        if (dat_.size() < 2)
        {
            throw std::domain_error("Must be at least two data points.");
        }
        if (dat_[0].size() != 2)
        {
            throw std::domain_error("Each datum must contain (y, y'), and nothing else.");
        }
        if (dx <= 0)
        {
            throw std::domain_error("dx > 0 is required.");
        }

        for (auto & d : dat_)
        {
            d[1] *= dx;
        }
    }

    inline Real operator()(Real x) const
    {
        const Real xf = x0_ + (dat_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf)
        {
            return dat_.back()[0];
        }
        return this->unchecked_evaluation(x);
    }

    inline Real unchecked_evaluation(Real x) const
    {
        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(dat_.size())>(ii);

        Real t = s - ii;
        // If we had infinite precision, this would never happen.
        // But we don't have infinite precision.
        if (t == 0)
        {
            return dat_[i][0];
        }
        Real y0 = dat_[i][0];
        Real y1 = dat_[i+1][0];
        Real dy0 = dat_[i][1];
        Real dy1 = dat_[i+1][1];

        Real r = 1-t;
        return r*r*(y0*(1+2*t) + dy0*t)
              + t*t*(y1*(3-2*t) - dy1*r);
    }

    inline Real prime(Real x) const
    {
        const Real xf = x0_ + (dat_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf)
        {
            return dat_.back()[1]*inv_dx_;
        }
        return this->unchecked_prime(x);
    }

    inline Real unchecked_prime(Real x) const
    {
        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(dat_.size())>(ii);
        Real t = s - ii;
        if (t == 0)
        {
            return dat_[i][1]*inv_dx_;
        }
        Real y0 = dat_[i][0];
        Real dy0 = dat_[i][1];
        Real y1 = dat_[i+1][0];
        Real dy1 = dat_[i+1][1];

        Real dy = 6*t*(1-t)*(y1 - y0)  + (3*t*t - 4*t+1)*dy0 + t*(3*t-2)*dy1;
        return dy*inv_dx_;
    }


    Size size() const
    {
        return dat_.size();
    }

    int64_t bytes() const
    {
        return dat_.size()*dat_[0].size()*sizeof(Real) + sizeof(dat_) + 2*sizeof(Real);
    }

    std::pair<Real, Real> domain() const
    {
        Real xf = x0_ + (dat_.size()-1)/inv_dx_;
        return {x0_, xf};
    }


private:
    RandomAccessContainer dat_;
    Real x0_;
    Real inv_dx_;
};

}
}
}
}
#endif
