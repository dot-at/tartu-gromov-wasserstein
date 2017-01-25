// estimate-gromov-wasserstein.cc

#include <fstream>
#include <ctime>

#include "gromov_wasserstein.hh"
#include "memoli.hh"

#include "gurobi_c++.h"

extern "C"
{
#include <unistd.h> // for getopt()
}

#include <cstdlib>

const char * const usgMsg =
    "Usage:\n"
    "estimate-gromov-wasserstein [-h]\n"
    "                            [-p#] [-d]  [-q | -v] [-l] [-O$]\n"
    "                            [-i#] [-j#] [-s#] [-w]\n"
    "                            file1 file2\n"
    "\n"
    "Compute the Gromov-Wasserstein p-distance\n"
    " where file1, file2 are names of mmspace files.\n"
    "General options:\n"
    "  -h      write longer help message to stdout\n"
    "  -p#     sets the p; # is float in ]0,infty[ (default 1)\n"
    "  -d      scale diameter of metric spaces to 1\n"
    "  -q      quiet, only print stats\n"
    "  -v      verbose output\n"
    "  -l      turn on logging of the LP solver (to log file)\n"
    "  -o$     ouput: write resulting PDF to file $\n"
    "Options for Memoli's heuristic:\n"
    "  -i#     max number of inner iterations (default 25)\n"
    "  -j#     number of outer iterations (default 3)\n"
    "  -s#     random seed\n"
    "  -w      warm-start LPs in inner iterations\n";

const char* more_help =
    "The mmspace File Format\n"
    "-----------------------\n"
    "The file defines a metric measure space on the point set {0,...,n-1}.\n"
    "The format for an explicit distribution q is this (lines matter):\n"
    "   'n' n\n"
    "   'g' q_0 q_1 p_2 ... q_{n-1}\n"
    "   'd' d(1,0)\n"
    "   'd' d(2,0)   d(2,1)\n"
    "   'd' d(3,0)   d(3,1) d(3,2)\n"
    "   ...\n"
    "   'd' d(n-1,0) ...    d(n-1,n-2)\n"
    "where '*' refers to the one-character string literal (in other words,\n"
    "omit the apostrophes).\n"
    "If the distribution is uniform, replace the 2nd line by\n"
    "   'u'\n"
    "The numbers are white space separated.\n"
    "\n"
    "The Output File\n"
    "---------------\n"
    "If -Ooutfile is given, the file \"outfile\" is created/truncated,\n"
    "and the probability density function P giving the distance estimate\n"
    "is written to \"outfile\" in one the following formats.\n"
    "Dense format:\n"
    "   'D' n_1 n_2\n"
    "   'p' P(0,0)     P(0,1)    ...  P(0,n_2-1)\n"
    "   'p' P(1,0)     P(1,1)    ...  P(1,n_2-1)\n"
    "   ... \n"
    "   'p' P(n_1-1,0) P(n_1-1,1) ... P(n_1-1,n_2-1)\n"
    "Sparse format:\n"
    "   'S' n_1 n_2 m\n"
    "   'q' x_1 y_1 P(x_1,y_1)\n"
    "   'q' x_2 y_2 P(x_2,y_2)\n"
    "   ... \n"
    "   'q' x_m y_m P(x_m,y_m)\n"
    "\n"
    "The Statistics\n"
    "--------------\n"
    "If -q is given a CSV-format row is printed to stdout.\n"
    "The cells are the following:\n"
    "      'Memoli'{str},success{'y','n'},time{float(sec)},\n"
    "      distance{float},marginal_violation{float},\n"
    "      nonnegativity_violation{float},nnz{int}\n"
    "If success is 'n' (i.e., failure), the row ends after time.\n";





typedef Gromov_Wasserstein::MMSpace_from_Stream          MMSpace1;
typedef Gromov_Wasserstein::MMSpace_from_Stream          MMSpace2;
typedef Gromov_Wasserstein::Memoli<MMSpace1,MMSpace2>    Memoli;

int main(int argc, char *argv[]) try
{
    if (argc == 1) {
        std::cerr<<usgMsg;
        return 1;
    }

    // Parameters
    double      p                 = 1.;     // -p:
    bool        d_scale           = false;  // -d
    int         output_level      = 0;      // -q -v
    bool        gurobi_log        = false;  // -l
    std::string outfilename       = "";     // -O:
    int         memoli_inner_iter = 25;     // -i:
    int         memoli_outer_iter = 3;      // -j:
    int         random_seed       = 12345;  // -s:
    bool        warm_start_LPs    = false;  // -w

    const char * options = "h" "p:" "d" "q" "v" "l" "O:" "i:" "j:" "s:" "w";
    for (char arg=getopt(argc, argv, options); arg!=-1; arg=getopt(argc, argv, options) ) {
        char * str_end;
        switch (arg) {
        case 'h': std::cout<<usgMsg<<'\n'<<more_help<<std::endl;                                                                                                  return 0;
        case 'p': p              = std::strtod(optarg,&str_end);        if(*str_end!='\0' || p<=0) return std::cerr<<"Invalid use of -p.\n", 1;                   break;
        case 'd': d_scale        = true;                                                                                                                          break;
        case 'q': output_level   = -1;                                                                                                                            break;
        case 'v': output_level   = +1;                                                                                                                            break;
        case 'l': gurobi_log     = true;                                                                                                                          break;
        case 'O': outfilename    = optarg;                                                                                                                        break;

        case 'i': memoli_inner_iter = std::strtol(optarg,&str_end,10);  if(*str_end!='\0' || memoli_inner_iter<1) return std::cerr<<"Invalid use of -i.\n", 1;    break;
        case 'j': memoli_outer_iter = std::strtol(optarg,&str_end,10);  if(*str_end!='\0' || memoli_outer_iter<1) return std::cerr<<"Invalid use of -j.\n", 1;    break;
        case 's': random_seed       = std::strtol(optarg,&str_end,10);  if(*str_end!='\0') return std::cerr<<"Invalid use of -s.\n", 1;                           break;
        case 'w': warm_start_LPs    = true;                                                                                                                       break;

        default:
            std::cerr<<"This is a bug in the parameter parsing part, sorry.  Write and email do dotheis@ut.ee";
            return -1;
        }
    }
    // check non-optional options :)
    if (p<=1.e-10 || p >1.e100)  return std::cerr<<"p out of range.\n"<<usgMsg, 2;
    if (optind != argc-2)        return std::cerr<<"Please give file names.\n"<<usgMsg, 2;

    const char * const fn1 = argv[argc-2];
    const char * const fn2 = argv[argc-1];

    std::ifstream file1 (fn1);
    std::ifstream file2 (fn2);

    if (! file1) {
        std::cerr<<"Could not open 1st input file: "<<fn1<<"\n";
        return 8;
    }
    if (! file2) {
        std::cerr<<"Could not open 2nd input file: "<<fn2<<"\n";
        return 8;
    }

    std::ofstream pdf_out;
    if (outfilename!="") {
        pdf_out.open(outfilename,std::ofstream::trunc);
        if (!pdf_out) {
            std::cerr<<"Could not open output file: "<<outfilename<<"\n";
            return 8;
        }
    }


    if (output_level>=0) {
        std::cout<<"estimate-gromov-wasserstein\n"
                 <<"Computing Memoli's heuristic upper bound on the Gromov-Wasserstein "<<p<<"-distance\n"
                 <<"between metric measure spaces in '"<<fn1<<"' and '"<<fn2<<"'"
                 <<""; // no newline
        if (outfilename!="") std::cout<<",\nwriting resulting PDF to "<<outfilename<<".\n";
        else                 std::cout<<".\n";
    }

    std::srand(random_seed);

    MMSpace1 X1 (file1, d_scale);
    MMSpace1 X2 (file2, d_scale);
    if (output_level>=0) {
        std::cout<<"Created metric measure spaces:\n"
                 <<"X1 with "<<X1.n()<<" points with "<<(X1.is_uniform()? "uniform":"explicit")<<" distribution;\n"
                 <<"X2 with "<<X2.n()<<" points with "<<(X2.is_uniform()? "uniform":"explicit")<<" distribution.\n";
    }

    int       stat_success;
    double    stat_distance           = -1.;
    double    stat_time;
    double    stat_marginal_violation = -1;
    double    stat_nonnegat_violation = -1;
    unsigned  stat_nnz                = -1;
    std::string return_string;

    const std::clock_t t_0 = std::clock();
    std::ostream devnull (0);
    Memoli  memoli (X1,X2,p, memoli_outer_iter, memoli_inner_iter,  warm_start_LPs,
                    (output_level > 0 ? std::cout : devnull), (gurobi_log? "estimate-gromov-wasserstein--LP-solver.log" : (const char *)NULL) );
    int     retval = memoli.compute_it(&return_string);
    const std::clock_t t_1 = std::clock();

    stat_time = (t_1-t_0)/(double)CLOCKS_PER_SEC;

    if (output_level>=0 && return_string.size()) std::cout<<"Solver return message:\n"<<return_string<<std::endl;

    if (retval >= 0) {
        stat_success  = 1;
        stat_distance = memoli.the_distance(&stat_marginal_violation,
                                            &stat_nonnegat_violation,
                                            &stat_nnz);
        if (output_level < 0) {
            std::cout<<"Memoli,"<<fn1<<','<<fn2
                     <<','<<(stat_success? 'y' : 'n')
                     <<','<<stat_time
                     <<','<<stat_distance
                     <<','<<stat_marginal_violation
                     <<','<<stat_nonnegat_violation
                     <<','<<stat_nnz
                     <<std::endl;
        } else {
            std::cout<<"Done!\n"
                     <<"Time:                         "<<stat_time<<"s\n"
                     <<"Distance:                     "<<stat_distance<<'\n'
                     <<"Total marginal violation:     "<<stat_marginal_violation<<'\n'
                     <<"Total nonnegatvity violation: "<<stat_nonnegat_violation<<'\n'
                     <<"Number of nonzeros:           "<<stat_nnz<<std::endl;
        } //^ if-else quiet

        if (outfilename!="") memoli.write_pdf(pdf_out);
        pdf_out<<std::flush;
        pdf_out.close();

    } else {
        stat_success = 0;
        if (output_level < 0) {
            std::cout<<"Memoli,"
                     <<(stat_success? 'y' : 'n')
                     <<std::endl;;
        } else {
            std::cout<<"Failure to compute GW distance.\n"
                     <<"You wasted "<<stat_time<<" seconds of your life.\n"
                     <<"Check if the reason is on your side (e.g., missing Gurobi license)\n"
                     <<"Otherwise, please send your input files and parameters to dotheis@ut.ee,\n"
                     <<"with subject line 'estimate-gromov-wasserstein problem'\n";
        }
    } //^ if-else success
    return 0;
}
catch(const Gromov_Wasserstein::GW_exception &gwe) {
    std::cout<<std::endl;
    std::cerr<<gwe.what();
    std::cerr<<std::endl;
}
catch(const std::exception &stde) {
    std::cout<<std::endl;
    std::cerr<<"The following std::exception was caught by main():\n";
    std::cerr<<stde.what();
    std::cerr<<std::endl;
}
catch(GRBException &grbe) {
    std::cout<<std::endl;
    std::cerr<<"The following GRBException was caught by main():\n";
    std::cerr<<grbe.getMessage();
    std::cerr<<std::endl;
}
catch(const std::string &stre) {
    std::cout<<std::endl;
    std::cerr<<"The following std::string was caught as an exception by main():\n";
    std::cerr<<stre;
    std::cerr<<std::endl;
}
catch(const char *cstr) {
    std::cout<<std::endl;
    std::cerr<<"The following const-char[] exception was caught as an exception by main():\n";
    std::cerr<<cstr;
    std::cerr<<std::endl;
}
catch(...) {
    std::cout<<std::endl;
    std::cerr<<"Unknown exception caught by main()\n";
}
//^ main()
