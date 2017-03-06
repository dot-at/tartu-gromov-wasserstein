// memoli.hh
#ifndef __memoli_hh__
#define __memoli_hh__

#include "gromov_wasserstein.hh"

#include <valarray>
#include <stdexcept>

#include <iostream>

class GRBEnv;
class GRBModel;
class GRBVar;

namespace Gromov_Wasserstein {


    template<class M1, class M2>
    class Memoli: public GW_distance<M1,M2>
    {
        using GW_distance<M1,M2>::X1;
        using GW_distance<M1,M2>::X2;
        using GW_distance<M1,M2>::p;
        using GW_distance<M1,M2>::U;
        using GW_distance<M1,M2>::distance;
        using GW_distance<M1,M2>::pth_power_difference;
        using GW_distance<M1,M2>::index;

        typedef typename GW_distance<M1,M2>::Sparse_Entry Sparse_Entry;

    protected:
        typedef   typename M1::Point   Point1;
        typedef   typename M2::Point   Point2;


        std::ostream & log;
        const char *   gurobi_log_file;
        const int      gurobi_threads;
        const double   eps;
        const int      n_outer_iter;
        const int      n_inner_iter;
        const bool     warm_start_LPs;


        GRBEnv   * p_genv  =0;
        GRBModel * p_model =0;

        GRBVar * var       =0;

        std::vector< Sparse_Entry > current_LP_sol;
        double make_current_LP_sol();   // retrieves the current LP solution x from Gurobi and stores it into current_LP_sol;
                                        // then computes and returns x^T G x.

        void add_marginal_equations();  // add all marginal equations

        void make_initial_obj();        // writes initial objective function into ``obj[]''
        // void make_next_obj();        obsolete // reads GRB var values and writes next objective function into ``obj[]''
        void make_next_obj__sparse();   // same, for sparse LP solutions

        bool check_GRB_error(std::string *p_error);

        void store_current_LP_solution(); // stores the solution current_LP_sol into GW_distance::u

    public:
        virtual int compute_it(std::string *p_return_info);  // return value positive if success; negative if failure.

        Memoli(const M1 & _X1, const M2 & _X2, double p,
               int n_outer_iterations = 10, int n_inner_iterations = 100,
               bool warm_start_inner_LPs = false,
               std::ostream& _log=std::cout, const char * gurobi_log_file=0,
               int gurobi_threads=0,
               double eps=1.e-100);
        ~Memoli();
    };

} // namespace

#endif
// EOF memoli.hh
