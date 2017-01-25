// memoli.t.cc
#ifndef __memoli_t_cc__
#define __memoli_t_cc__

#include "memoli.hh"

#include "gurobi_c++.h"

#include <ctime>
#include <cstdlib>
#include <iomanip>


template<class M1, class M2>
Gromov_Wasserstein::Memoli<M1,M2>::Memoli(const M1 & _X1, const M2 & _X2, double p,
                                          int n_outer_iterations, int n_inner_iterations,
                                          bool warm_start,
                                          std::ostream& _log, const char * _gurobi_log_file,
                                          double _eps)
    :
    GW_distance<M1,M2>(_X1,_X2,p),
    log{_log},
    gurobi_log_file{_gurobi_log_file},
    eps{_eps},
    n_outer_iter{n_outer_iterations},
    n_inner_iter{n_inner_iterations},
    warm_start_LPs{warm_start}
{
    GW_distance<M1,M2>::solutions_are_sparse();
    current_LP_sol.reserve(X1.n()*X2.n());
    if (n_outer_iterations<=0 || n_inner_iterations<=0 ) throw std::invalid_argument("Gromov_Wasserstein::Memoli(constructor): invalid # iterations.");
}




template<class M1, class M2>
void Gromov_Wasserstein::Memoli<M1,M2>::make_initial_obj()
{
    // Initial objective is: absolute value of   || d_1(x_1,*) ||_p  - || d_2(x_2,*) ||_p
    // ... which BTW explains why it's so good: it matches the vertices of the same average distance to the other vertices.
    // Create the objective function for the initial solution u_0
    double s1[X1.n()];
    for (Point1 x1; x1<X1.end(); ++x1) {
        double sum = 0.;
        for (Point1 y1; y1<X1.end(); ++y1) {
            sum += X1.pr(y1) * std::pow( X1.d(x1,y1) , p) ;
        }
        s1[x1.x] = std::pow( sum, 1/p );
    }

    double s2[X2.n()];
    for (Point2 x2; x2<X2.end(); ++x2) {
        double sum = 0.;
        for (Point2 y2; y2<X2.end(); ++y2) {
            sum += X2.pr(y2) * std::pow( X2.d(x2,y2) , p) ;
        }
        s2[x2.x] = std::pow( sum, 1/p );
    }

    for (Point2 x2; x2<X2.end(); ++x2) {
        for (Point1 x1; x1<X1.end(); ++x1) {
            var[index(x1,x2)].set(GRB_DoubleAttr_Obj,  std::fabs( s1[x1.x] - s2[x2.x] ) );
        }
    }
} //^ make_initial_obj()


template<class M1, class M2>
double Gromov_Wasserstein::Memoli<M1,M2>::make_current_LP_sol()
{
    // make the sparse vector
    current_LP_sol.clear();
    for (Point1 y1; y1<X1.end(); ++y1) {
        for (Point2 y2; y2<X2.end(); ++y2) {
            const unsigned i = index(y1,y2);
            const double val = var[i].get(GRB_DoubleAttr_X);
            if (val > 0 ){ //eps
                Sparse_Entry e;
                e.u1  = y1;
                e.u2  = y2;
                e.val = val;
                current_LP_sol.push_back(e);
            } //^ if nonzero
        } //^ for y1
    } //^ for y2

    // compute the value
    long double sum = 0.;
    for (auto e : current_LP_sol) {
        for (auto f : current_LP_sol) {
            const long double c  = pth_power_difference(e.u1,e.u2,f.u1,f.u2);
            sum += c * e.val * f.val;
        } //^ for f
    } //^ for e

    // {
    //     long double sum_check = 0.;
    //     for(Point1 x1; x1<X1.end(); ++x1)
    //         for(Point2 x2; x2<X2.end(); ++x2)
    //             for(Point1 u1; u1<X1.end(); ++u1)
    //                 for(Point2 u2; u2<X2.end(); ++u2)
    //                     sum_check += var[index(x1,x2)].get(GRB_DoubleAttr_X)*var[index(u1,u2)].get(GRB_DoubleAttr_X)*pth_power_difference(u1,u2,x1,x2);
    //
    //     if (std::fabs(sum_check - sum) > 1.e-3) log<<"Gromov_Wasserstein::Memoli: VALUE MISMATCH: sum="<<sum<<" != "<<sum_check<<"=sum_check"<<std::endl;
    // }
    //
    // log<<"Gromov_Wasserstein::Memoli: LP solution w/ "<<current_LP_sol.size()<<" nonzeros; density="<<current_LP_sol.size()/(double)X1.n()/X2.n()*100.<<"%\n";
    // log<<"Gromov_Wasserstein::Memoli: Value = "<<sum<<'\n';

    // return the value
    return sum;
} //^ make_current_LP_sol()




template<class M1, class M2>
void Gromov_Wasserstein::Memoli<M1,M2>::make_next_obj__sparse()
{
    // Matrix-vector multiplication
    for (Point1 x1; x1<X1.end(); ++x1) {
        for (Point2 x2; x2<X2.end(); ++x2) {
            long double sum = 0.;
            for (auto e : current_LP_sol) {
                const long double c  = pth_power_difference(x1,x2,e.u1,e.u2);
                sum += c * e.val;
            } //^ current_LP_sol
            var[index(x1,x2)].set(GRB_DoubleAttr_Obj, sum);
        } //^ x1
    } //^ x2
} //^ make_next_obj__sparse()

template<class M1, class M2>
void Gromov_Wasserstein::Memoli<M1,M2>::add_marginal_equations()
{
    // X1 marginals
    for (Point1 x1; x1<X1.end(); ++x1) {
        // marginal: \sum_{x_2 \in X_2} \mu(x_1,x_2) = \mu(x_1)
        GRBLinExpr lhs;
        for (Point2 x2; x2<X2.end(); ++x2)   lhs += var[ index(x1,x2) ];
        p_model->addConstr( lhs == X1.pr(x1) );
    }

    // X2 marginals
    for (Point2 x2; x2<X2.end(); ++x2) {
        // marginal: \sum_{x_1 \in X_1} \mu(x_1,x_2) = \mu(x_2)
        GRBLinExpr lhs;
        for (Point1 x1; x1<X1.end(); ++x1)   lhs += var[ index(x1,x2) ];
        p_model->addConstr( lhs == X2.pr(x2) );
    }
} //^ add_marginal_equations()

template<class M1, class M2>
bool Gromov_Wasserstein::Memoli<M1,M2>::check_GRB_error(std::string *p_error)
{
    const int status = p_model->get(GRB_IntAttr_Status);
    switch (status) {
    case GRB_OPTIMAL: return false; // p_model->get(GRB_DoubleAttr_ObjVal);

    case GRB_LOADED:          *p_error="LOADED: Model is loaded, but no solution information is available."; return -1;
        // case GRB_OPTIMAL:         *p_error="OPTIMAL: Model was solved to optimality (subject to tolerances), and an optimal solution is available."; return -1;
    case GRB_INFEASIBLE:      *p_error="INFEASIBLE: Model was proven to be infeasible."; return -1;
    case GRB_INF_OR_UNBD:     *p_error="INF_OR_UNBD: Model was proven to be either infeasible or unbounded. To obtain a more definitive conclusion, set the DualReductions parameter to 0 and reoptimize."; return -1;
    case GRB_UNBOUNDED:       *p_error="UNBOUNDED: Model was proven to be unbounded."; return -1;
    case GRB_CUTOFF:          *p_error="CUTOFF: Optimal objective for model was proven to be worse than the value specified in the Cutoff parameter. No solution information is available."; return -1;
    case GRB_ITERATION_LIMIT: *p_error="ITERATION_LIMIT: Optimization terminated because the total number of simplex iterations performed exceeded the value specified in the IterationLimit parameter, or because the total number of barrier iterations exceeded the value specified in the BarIterLimit parameter."; return -1;
    case GRB_NODE_LIMIT:      *p_error="NODE_LIMIT: Optimization terminated because the total number of branch-and-cut nodes explored exceeded the value specified in the NodeLimit parameter."; return -1;
    case GRB_TIME_LIMIT:      *p_error="TIME_LIMIT: Optimization terminated because the time expended exceeded the value specified in the TimeLimit parameter."; return -1;
    case GRB_SOLUTION_LIMIT:  *p_error="SOLUTION_LIMIT: Optimization terminated because the number of solutions found reached the value specified in the SolutionLimit parameter."; return -1;
    case GRB_INTERRUPTED:     *p_error="INTERRUPTED: Optimization was terminated by the user."; return -1;
    case GRB_NUMERIC:         *p_error="NUMERIC: Optimization was terminated due to unrecoverable numerical difficulties."; return -1;
    case GRB_SUBOPTIMAL:      *p_error="SUBOPTIMAL: Unable to satisfy optimality tolerances; a sub-optimal solution is available."; return -1;
    case GRB_INPROGRESS:      *p_error="INPROGRESS: An asynchronous optimization call was made, but the associated optimization run is not yet complete."; return -1;

    default: *p_error="unknown status code"; return true;
    }
} //^ check_GRB_status()

template<class M1, class M2>
void Gromov_Wasserstein::Memoli<M1,M2>::store_current_LP_solution()
{
    GW_distance<M1,M2>::U_clear();
    for (auto e : current_LP_sol)  GW_distance<M1,M2>::U_append(e);
    distance = -1.; // invalidate distance
} //^ store_current_LP_solution()

template<class M1, class M2>
int Gromov_Wasserstein::Memoli<M1,M2>::compute_it(std::string *p_return_info)
{
    if (p_return_info) *p_return_info = "";

    // P r e p a r a t i o n s
    // =======================
    const unsigned sz = X1.n()*X2.n();

    p_genv  = new GRBEnv;
    p_model = new GRBModel{*p_genv};

    p_model->set(GRB_IntParam_LogToConsole, 0);
    p_model->set(GRB_StringParam_LogFile, (gurobi_log_file? gurobi_log_file : "") );

    // Prepare for initial LP
    var = p_model->addVars(nullptr, // lb=0
                           nullptr, // ub=infty
                           nullptr, // set objective later
                           nullptr, // type=continuous
                           nullptr, // no names
                           sz);
    p_model->update();
    add_marginal_equations();
    p_model->update();

    double      current_solution_value;
    double      best_solution_value = 1.e300;
    std::string grb_error_msg;

    // ****************************************************************************************************
    // M a i n   p r o c e d u r e   ---  n e w
    // =========================================
    for (int o=0; o<n_outer_iter; ++o) {
        for (int i=0; i<n_inner_iter; ++i) {
            log<<"Gromov_Wasserstein::Memoli: LP # "<<std::setw(4)<<o<<'|'<<std::setw(4)<<i<<": 1-"<<std::flush;

            // Phase 1: Set up objective function
            if (i==0) { // start new outer iteration
                if (o==0) {
                    make_initial_obj();
                    p_model->update();
                } else {
                    // give random costs in {-4,...,-1, +1,...,+4}
                    for (unsigned v=0; v<sz; ++v) {
                        int rnd_cost = std::rand()%8 - 4;
                        if (rnd_cost==0) ++rnd_cost;
                        var[v].set(GRB_DoubleAttr_Obj, (double)rnd_cost);
                    }
                }
            } else {
                make_next_obj__sparse();
            }

            // Phase 2: Solve the LP
            if (!warm_start_LPs) {
                p_model->reset();
                log<<"2+"<<std::flush;
            }
            else log<<"2-"<<std::flush;
            p_model->optimize();

            if (check_GRB_error(&grb_error_msg)) {
                if (p_return_info) *p_return_info = std::string("Gromov_Wasserstein::Memoli: Failed to solve LP #")+std::to_string(i+1)+" ("+grb_error_msg+"). Sorry.";
                return -1;
            }

            // Phase 3: Read out the solution from the LP variables
            log<<"3"<<std::flush;
            current_solution_value  =   make_current_LP_sol();

            // Phase 4: Check if there's a new best solution
            bool new_best = false;
            if (current_solution_value < best_solution_value) {
                log<<"-4"<<std::flush;
                best_solution_value = current_solution_value;
                store_current_LP_solution();
                new_best=true;
            } else log<<"  ";

            // Criteria for BREAKING out of inner loop:
            bool break_inner_loop = false;
            if (p_model->get(GRB_DoubleAttr_IterCount) == 0)  break_inner_loop=true;
            // More sophistic criteria here
            // ...

            // Output of LP stats:
            log<<" done"
               <<"; Value:"<<(new_best? "* " : "  ")<<std::setw(9)<<current_solution_value
               <<"; nnz: "<<std::setw(6)<<current_LP_sol.size()
               <<"; Simplex iters: "<<std::setw(9)<<p_model->get(GRB_DoubleAttr_IterCount)
               <<"; Simplex time: "<<std::setw(9)<<p_model->get(GRB_DoubleAttr_Runtime)<<"s"
               <<(break_inner_loop?  "; BREAK INNER LOOP" : "; continue" )
               <<std::endl;


            if (break_inner_loop)  break;
        } //^ for i (inner loop)

    } //^ for o (outer loop)
    log<<"Gromov_Wasserstein::Memoli: Done"<<std::endl;

    return 1;
} //^ compute_it()


template<class M1, class M2>
Gromov_Wasserstein::Memoli<M1,M2>::~Memoli()
{
    delete[] var;
    delete p_model;
    delete p_genv;
} //^ ~Memoli


namespace Gromov_Wasserstein {
    // template class Memoli<Simple_MSpace, Simple_MSpace>;
    // template class Memoli<Simple_MSpace, Simple_MMSpace>;
    // template class Memoli<Simple_MMSpace,Simple_MSpace>;
    // template class Memoli<Simple_MMSpace,Simple_MMSpace>;
    template class Memoli<MMSpace_from_Stream,MMSpace_from_Stream>;
}

#endif
// EOF
