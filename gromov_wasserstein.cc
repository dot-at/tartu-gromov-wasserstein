// gromov_wasserstein.cc

#include "gromov_wasserstein.hh"

Gromov_Wasserstein::MMSpace_from_Stream::MMSpace_from_Stream(std::istream & in, bool scale_diam)
{
    char c;
    int  i;
    double q;
    double max_d;

    in>>std::ws>>c;
    if (c!='n')                           throw GW_exception(std::string("MMSpace_from_Stream: Expected n, got ")+c);
    in>>i;
    if (in.fail())                        throw GW_exception(std::string("MMSpace_from_Stream: Could not read number of points"));
    if (i<=2 || i > 999999)               throw GW_exception(std::string("MMSpace_from_Stream: Expected number of points, got ")+std::to_string(i));
    _n = i;

    in>>std::ws>>c;
    switch(c) {

    case 'u':
        uniform = true;
        break;

    case 'g':
        uniform = false;
        mu.resize(_n);
        for (Point x; x<end(); ++x) {
            in >> q;
            if (in.fail())                throw GW_exception(std::string("MMSpace_from_Stream: Could not read p(")+std::to_string(x.x)+")");
            if (q<0 || q >= 1)            throw GW_exception(std::string("MMSpace_from_Stream: Expected p(")+std::to_string(x.x)+"), got "+std::to_string(q));
            mu[x.x] = q;
        }
        break;


    default:                              throw GW_exception(std::string("MMSpace_from_Stream: Expected 'u' (for uniform distribution) or 'g' (for explicit distribution)"));
    }

    D.resize(_n*_n);
    max_d = -1;
    Point x;
    d(x.x,x.x) = 0.;
    while (++x<end()) {
        d(x.x,x.x) = 0.;
        in>>std::ws>>c;
        if (c!='d')                       throw GW_exception(std::string("MMSpace_from_Stream: Expected d, got ")+c);
        for (Point y; y<x; ++y) {
            in >> q;
            if (in.fail())                throw GW_exception(std::string("MMSpace_from_Stream: Could not read d(")+std::to_string(x.x)+","+std::to_string(y.x)+")");
            if (q<0 || q > 1.e20)         throw GW_exception(std::string("MMSpace_from_Stream: Expected d(")+std::to_string(x.x)+","+std::to_string(y.x)+"), got "+std::to_string(q));
            d(x.x,y.x) = q;
            d(y.x,x.x) = q;
            max_d = std::max(max_d,q);
        } //^ for y
    } //^ while x

    if (max_d < 1.e-8)                    throw GW_exception(std::string("MMSpace_from_Stream: Diameter is ")+std::to_string(max_d)+" --- that can't be right!");

    if (scale_diam) {
        for (Point x; x<end(); ++x) {
            for (Point y; y<x; ++y) {
                d(x.x,y.x) /= max_d;
                d(y.x,x.x) /= max_d;
            }
        }
    } //^ if scale_diam

} //^ MMSpace_from_Stream -- constructor
