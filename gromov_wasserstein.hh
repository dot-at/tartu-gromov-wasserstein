// gromov_wasserstein.hh
#ifndef __gromov_wasserstein_hh__
#define __gromov_wasserstein_hh__

// #include "dumb_digraph.hh"


#include <valarray>
#include <stdexcept>
#include <iostream>

namespace Gromov_Wasserstein {

// ##########################################################################################################################################################################
//  M e t r i c   M e a s u r e   S p a c e s
// ##########################################################################################################################################################################

    struct GW_exception: public std::runtime_error {
        GW_exception(std::string what): std::runtime_error(what) {}
    };

    // ******************************************************************************************************************************************************
    // MMSpace
    // ******************************************************************************************************************************************************
    class Simple_MMSpace
    {
    protected:
        unsigned                _n;
        std::valarray<double>   D;
        std::valarray<double>   mu;


    public:
        struct Point {
            unsigned x;
            inline Point(unsigned _x=0): x(_x) {}
            inline Point operator ++ () {  ++x; return *this; }
            inline bool operator == (Point q) const { return x == q.x; }
            inline bool operator <  (Point q) const { return x <  q.x; }
            inline bool operator <= (Point q) const { return x <= q.x; }
            inline bool operator >  (Point q) const { return x >  q.x; }
            inline bool operator >= (Point q) const { return x >= q.x; }
        };

        Simple_MMSpace(unsigned __n): _n(__n), D(-1.,_n*_n), mu(-1.,_n)   {}

        unsigned      n   ()                      const  { return _n; }
        Point         end ()                      const  { return Point(_n); }
        // size of underlying set

        double & d (unsigned x, unsigned y)         { if (x>=_n || y>=_n) throw std::domain_error("Gromov_Wasserstein::Simple_MMSpace::d(): out of range"); else return D[y*_n+x]; }
        double   d (unsigned x, unsigned y)  const  { if (x>=_n || y>=_n) throw std::domain_error("Gromov_Wasserstein::Simple_MMSpace::d(): out of range"); else return D[y*_n+x]; }
        double   d (Point    x, Point    y)  const  { return d(x.x,y.x); }
        // distance. x,y must be in {0,...,n-1}

        double & pr(unsigned x)                          { if (x>=_n) throw std::domain_error("Gromov_Wasserstein::Simple_MMSpace::p(): out of range"); else return mu[x]; }
        double   pr(unsigned x)                   const  { if (x>=_n) throw std::domain_error("Gromov_Wasserstein::Simple_MMSpace::p(): out of range"); else return mu[x]; }
        double   pr(Point    x)                   const  { return pr(x.x); }
        // measure of point x in {0,...,n-1}
        static bool     is_uniform()                     { return false; }
    };

    // ******************************************************************************************************************************************************
    // MSpace
    // ******************************************************************************************************************************************************
    class Simple_MSpace
    {
    protected:
        unsigned              _n;
        std::valarray<double> D;

    public:
        struct Point {
            unsigned x;
            inline Point(unsigned _x=0): x(_x) {}
            inline Point operator ++ () {  ++x; return *this; }
            inline bool operator == (Point q) const { return x == q.x; }
            inline bool operator <  (Point q) const { return x <  q.x; }
            inline bool operator <= (Point q) const { return x <= q.x; }
            inline bool operator >  (Point q) const { return x >  q.x; }
            inline bool operator >= (Point q) const { return x >= q.x; }
        };

        Simple_MSpace(unsigned __n): _n{__n}, D(-1.,_n*_n)   {}

        unsigned      n   ()                      const  { return _n; }
        Point         end ()                      const  { return Point(_n); }
        // size of underlying set

        double & d (unsigned x, unsigned y)              { if (x>=_n || y>=_n) throw std::domain_error("Gromov_Wasserstein::Simple_MMSpace::d(): out of range"); else return D[y*_n+x]; }
        double   d (unsigned x, unsigned y)       const  { if (x>=_n || y>=_n) throw std::domain_error("Gromov_Wasserstein::Simple_MMSpace::d(): out of range"); else return D[y*_n+x]; }
        double   d(Point     x, Point    y)       const  { return d(x.x,y.x); }
        // distance. x,y must be in {0,...,n-1}

        double   pr(unsigned x)                   const  { if (x>=_n) throw std::domain_error("Gromov_Wasserstein::Simple_MMSpace::p(): out of range"); else return 1./_n; }
        double   pr(Point    x)                   const  { return pr(x.x); }
        // measure of point x in {0,...,n-1}
        static bool     is_uniform()                     { return true; }
    };

    // ******************************************************************************************************************************************************
    // MSpace from Stream
    // ******************************************************************************************************************************************************
    class MMSpace_from_Stream
    {
    protected:
        unsigned                _n;
        std::valarray<double>   D;
        std::valarray<double>   mu;
        bool                    uniform;

        double & d (unsigned x, unsigned y)              { if (x>=_n || y>=_n) throw std::domain_error("Gromov_Wasserstein::Simple_MMSpace::d(): out of range"); else return D[y*_n+x]; }

    public:
        struct Point {
            unsigned x;
            inline Point(unsigned _x=0): x(_x) {}
            inline Point operator ++ () {  ++x; return *this; }
            inline bool operator == (Point q) const { return x == q.x; }
            inline bool operator <  (Point q) const { return x <  q.x; }
            inline bool operator <= (Point q) const { return x <= q.x; }
            inline bool operator >  (Point q) const { return x >  q.x; }
            inline bool operator >= (Point q) const { return x >= q.x; }
        };

        MMSpace_from_Stream(std::istream &, bool scale_diameter=false);

        unsigned n  ()                       const  { return _n; }
        Point    end()                       const  { return Point(_n); }
        // size of underlying set

        double   d(unsigned x, unsigned y)   const  { if (x>=_n || y>=_n) throw std::domain_error("Gromov_Wasserstein::Simple_MMSpace::d(): out of range"); else return D[y*_n+x]; }
        double   d(Point    x, Point    y)   const  { return d(x.x,y.x); }

        double   pr(unsigned x)              const  { if (x>=_n) throw std::domain_error("Gromov_Wasserstein::Simple_MMSpace::p(): out of range"); else return ( uniform? 1./_n : mu[x] ); }
        double   pr(Point    x)              const  { return pr(x.x); }

        bool     is_uniform()                const  { return uniform; }
    };


// ##########################################################################################################################################################################
//  C o m p u t a t i o n   o f    G - W   D i s t a n c e
// ##########################################################################################################################################################################

    // ******************************************************************************************************************************************************
    // GW_distance
    // ******************************************************************************************************************************************************
    template<class M1, class M2>
    class GW_distance
    {
    protected:
        typedef   typename M1::Point   Point1;
        typedef   typename M2::Point   Point2;

        struct Sparse_Entry {
            Point1 u1;
            Point2 u2;
            double val;
        };

        const M1 &   X1;
        const M2 &   X2;
        const double p;

        mutable double        distance;  // set to -1 to invalidate after changing u

    private:
        std::valarray<double>      u_dense;         // X2.n() by X1.n() matrix in row major (contiguous rows).
        std::vector<Sparse_Entry>  u_sparse;

        mutable double        total_marginal_violation;
        mutable double        total_negativity;
        mutable unsigned      n_nonzeros;

        enum class Sparse_Or_Dense {undef, sparse, dense};
        Sparse_Or_Dense            sparseornot;

    public:
        virtual int compute_it(std::string *p_return_info) = 0;  // return value positive if success; negative if failure.
                                                                 // The string contains info to be analyzed/printed by caller.

        virtual ~GW_distance() {}
        // may throw a GW_exception

        GW_distance(const M1 & _X1, const M2 & _X2, double _p):
            X1{_X1}, X2{_X2}, p{_p}, distance{-1.}, sparseornot(Sparse_Or_Dense::undef)  {if (p<=0) throw std::domain_error("Gromov_Wasserstein::GW_distance<M1,M2>[constructor]: p cannot be negative.");}

        inline
        long double pth_power_difference(unsigned x1, unsigned x2, unsigned y1, unsigned y2) const;
        long double pth_power_difference(typename M1::Point x1, typename M2::Point x2, typename M1::Point y1, typename M2::Point y2) const { return pth_power_difference(x1.x,x2.x,y1.x,y2.x); }

        inline
        unsigned    index (unsigned  x1, unsigned  x2)  const;
        unsigned    index (typename M1::Point x1, typename M2::Point x2)  const    { return index(x1.x,x2.x); }


        // dense solution access
        void        solutions_are_dense() { sparseornot = Sparse_Or_Dense::dense; u_dense.resize(X1.n()*X2.n()); }
        double      U     (unsigned  x1, unsigned  x2)                    const    { if(sparseornot==Sparse_Or_Dense::dense) return u_dense[index(x1,x2)]; else throw std::runtime_error("Gromov_Wasserstein::U(): not dense."); }
        double      U     (typename M1::Point x1, typename M2::Point x2)  const    { if(sparseornot==Sparse_Or_Dense::dense) return u_dense[index(x1,x2)]; else throw std::runtime_error("Gromov_Wasserstein::U(): not dense."); }
        double &    U     (typename M1::Point x1, typename M2::Point x2)           { if(sparseornot==Sparse_Or_Dense::dense) return u_dense[index(x1,x2)]; else throw std::runtime_error("Gromov_Wasserstein::U(): not dense."); }
        // const std::valarray<double> get_u()        const { return u; }

        // sparse solution access
        void        solutions_are_sparse()              { sparseornot = Sparse_Or_Dense::sparse; u_sparse.reserve(X1.n()*X2.n()); u_sparse.clear(); }
        void        U_clear()                           { if(sparseornot==Sparse_Or_Dense::sparse) u_sparse.clear();        else throw std::runtime_error("Gromov_Wasserstein::U(): not sparse."); }
        void        U_append(Sparse_Entry e)            { if(sparseornot==Sparse_Or_Dense::sparse) u_sparse.push_back(e);   else throw std::runtime_error("Gromov_Wasserstein::U(): not sparse."); }
        auto        U_allofit()                  const  { if(sparseornot==Sparse_Or_Dense::sparse) return u_sparse;         else throw std::runtime_error("Gromov_Wasserstein::U(): not sparse."); }
        auto        U_begin()                    const  { if(sparseornot==Sparse_Or_Dense::sparse) return u_sparse.begin(); else throw std::runtime_error("Gromov_Wasserstein::U(): not sparse."); }
        auto        U_end()                      const  { if(sparseornot==Sparse_Or_Dense::sparse) return u_sparse.end();   else throw std::runtime_error("Gromov_Wasserstein::U(): not sparse."); }

        double  the_distance       (double *p_tmv=0, double *p_ng=0, unsigned *p_nnz=0, double eps_nnz=1.e-99)   const;
        void    write_pdf          (std::ostream &)                                                              const;
    }; //^ GW_distance

} //^ namespace


// pth_power_difference()
template<class M1, class M2>
long double Gromov_Wasserstein::GW_distance<M1,M2>::pth_power_difference(unsigned x1, unsigned x2, unsigned y1, unsigned y2) const
{
    const long double d1 = X1.d(x1,y1);
    const long double d2 = X2.d(x2,y2);
    if (p==1)  return std::fabs( d2-d1 );
    if (p==2)  {
        long double a=std::fabs( d2-d1 );
        return a*a;
    }
    return std::pow( std::fabs( d2-d1 ) ,  (long double)p );
} //^ pth_power_difference()

// index()
template<class M1, class M2>
unsigned Gromov_Wasserstein::GW_distance<M1,M2>::index(unsigned x1, unsigned x2) const
{
#   ifdef Gromov_Wasserstein_DEBUG
    if (x1>=X1.n()) throw std::domain_error("Gromov_Wasserstein::GW_distance::index(): x1 out of range");
    if (x2>=X2.n()) throw std::domain_error("Gromov_Wasserstein::GW_distance::index(): x2 out of range");
#   endif
    return x1*X2.n() + x2;
} //^ index()


template<class M1, class M2>
double Gromov_Wasserstein::GW_distance<M1,M2>::the_distance(double *p_tmv, double *p_ng, unsigned *p_nnz, double eps) const
{
    if (distance < -.5) {
        long double sum;
        switch (sparseornot) {
        case Sparse_Or_Dense::dense:
        {
            sum = 0.;
            for(Point1 x1; x1<X1.end(); ++x1)
                for(Point2 x2; x2<X2.end(); ++x2)
                    for(Point1 y1; y1<X1.end(); ++y1)
                        for(Point2 y2; y2<X2.end(); ++y2) // note that the objective function diagonal is 0
                            sum += U(x1,x2)*U(y1,y2)*pth_power_difference(x1,x2,y1,y2);
            distance = std::pow( sum , 1/(long double)p );

            sum = 0.;
            for(Point1 x1; x1<X1.end(); ++x1) {
                long double marg = X1.pr(x1);
                for(Point2 x2=0; x2<X2.end(); ++x2)  marg -= U(x1,x2);
                sum += std::fabs(marg);
            }
            for(Point2 x2; x2<X2.end(); ++x2) {
                long double marg = X2.pr(x2);
                for(Point1 x1; x1<X1.end(); ++x1)  marg -= U(x1,x2);
                sum += std::fabs(marg);
            }
            total_marginal_violation = sum;

            sum = 0.;
            n_nonzeros = 0;
            for(Point1 x1; x1<X1.end(); ++x1) {
                for(Point2 x2; x2<X2.end(); ++x2) {
                    if ( U(x1,x2) > eps )    ++n_nonzeros;
                    else                     sum -= std::min(0.,U(x1,x2));
                }
            }
            total_negativity = sum;
        }
        break;
        case Sparse_Or_Dense::sparse:
        {
            sum = 0.;
            for (auto e : U_allofit())
                for (auto f: U_allofit())
                    sum += e.val*f.val*pth_power_difference(e.u1,e.u2,f.u1,f.u2);

            distance = std::pow( sum , 1/(long double)p );

            {
                std::vector<long double> p1 (X1.n());
                for(Point1 x1; x1<X1.end(); ++x1) p1[x1.x] = X1.pr(x1);

                std::vector<long double> p2 (X2.n());
                for(Point2 x2; x2<X2.end(); ++x2) p2[x2.x] = X2.pr(x2);

                for (auto e : U_allofit()) {
                    p1[e.u1.x] -= e.val;
                    p2[e.u2.x] -= e.val;
                }

                sum = 0;
                for(Point1 x1; x1<X1.end(); ++x1) sum += std::fabs( p1[x1.x] );
                for(Point2 x2; x2<X2.end(); ++x2) sum += std::fabs( p2[x2.x] );

                total_marginal_violation = sum;
            }

            sum = 0.;
            n_nonzeros = 0;
            for (auto e : U_allofit()) {
                if ( e.val > eps )    ++n_nonzeros;
                else                  sum -= std::min(0.,e.val);
            }
            total_negativity = sum;
        }
        break;

        default:
            throw std::runtime_error("Gromov_Wasserstein::the_distance(): not sparse not dense.");
        }
    }

    if (p_tmv) *p_tmv = total_marginal_violation;
    if (p_ng)  *p_ng  = total_negativity;
    if (p_nnz) *p_nnz = n_nonzeros;
    return distance;
} //^ the_distance()

template<class M1, class M2>
void Gromov_Wasserstein::GW_distance<M1,M2>::write_pdf(std::ostream & out) const
{
    switch (sparseornot) {
    case Sparse_Or_Dense::dense:
    {
        out<<"D "<<X1.n()<<' '<<X2.n()<<'\n';
        for(Point1 x1; x1<X1.end(); ++x1) {
            out<<"p ";
            for(Point2 x2; x2<X2.end(); ++x2) out << U(x1,x2) <<' ';
            out<<'\n';
        }
    }
    break;
    case Sparse_Or_Dense::sparse:
    {
        out<<"S "<<X1.n()<<' '<<X2.n()<<' '<<U_allofit().size()<<'\n';
        for (auto e : U_allofit())  out <<"q "<<e.u1.x<<' '<<e.u2.x<<' '<<e.val<<'\n';
    }
    break;

    default:
        throw std::runtime_error("Gromov_Wasserstein::write_pdf(): not sparse not dense.");

    }
} //^ write_solution_pdf()

#endif
// EOF gromov_wasserstein.hh
