/**
 * @file dSFMTdc.cpp
 *
 * @brief The main function of parameter generator of dSFMT.
 *
 * The functions in this file are simple. They parse the command line
 * options and call all_in_one function which does almost all things.
 * Users can change this file so that it fits to their applications
 * and OS.
 *
 * @author Mutsuo Saito (Manieth corp.)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2015 Mutsuo Saito, Makoto Matsumoto,
 * Hiroshima University.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <sstream>
#include <fstream>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/AlgorithmCalculateParity.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/MersenneTwister64.hpp>
#include <NTL/GF2X.h>
#include <getopt.h>
#include "dSFMTsearch.hpp"
#include "AlgorithmDSFMTEquidistribution.hpp"
#include "Annihilate.h"
#include "calc_fixpoint.h"
#include "mpicontrol.hpp"

using namespace std;
using namespace MTToolBox;
using namespace NTL;

class options {
public:
    int mexp;
    bool verbose;
    int fixedSL1;
    int fixedPOS1;
    uint64_t seed;
    std::string filename;
    long count;
};

bool parse_opt(options& opt, int argc, char **argv);


int search(options& opt, ostream& os, int count);

/**
 * parse command line option, and search parameters
 * @param argc number of arguments
 * @param argv value of arguments
 * @return 0 if this ends normally
 */
int main(int argc, char** argv) {
    MPIControl mpi(&argc, &argv);
    options opt;
    bool parse = parse_opt(opt, argc, argv);
    if (!parse) {
        return -1;
    }
    ofstream ofs;
// MPI OUTPUT CHANGE
#if !defined(NO_MPI)
    // MPI では .txt を付けずに指定する。それが一番簡単
    opt.seed = opt.seed + mpi.getRank();
    if (!opt.filename.empty()) {
        char buff[200];
        sprintf(buff, ".s%04ld-%03d.txt", opt.seed, mpi.getRank());
        opt.filename += buff;
    }
#endif
    if (!opt.filename.empty()) {
        ofs.open(opt.filename.c_str());
        if (!ofs) {
            cerr << "can't open file:" << opt.filename << endl;
            return -1;
        }
        return search(opt, ofs, opt.count);
    } else {
        return search(opt, cout, opt.count);
    }
}

/**
 * search parameters using all_in_one function in the file search_all.hpp
 * @param opt command line options
 * @param count number of parameters user requested
 * @return 0 if this ends normally
 */
int search(options& opt, ostream& os, int count) {
    MersenneTwister64 mt(opt.seed);
    dSFMT g(opt.mexp);

    os << "seed = " << dec << opt.seed << endl;
    if (opt.verbose) {
        time_t t = time(NULL);
        os << "search start at " << ctime(&t);
    }
    if (opt.fixedSL1 > 0) {
        g.setFixedSL1(opt.fixedSL1);
    }
    if (opt.fixedPOS1 > 0) {
        g.setFixedPOS1(opt.fixedPOS1);
    }
    AlgorithmReducibleRecursionSearch<w128_t> ars(g, mt);
    int i = 0;
    AlgorithmCalculateParity<w128_t, dSFMT> cp;
    os << "# " << g.getHeaderString() << ", delta52"
         << endl;
    while (i < count) {
        if (ars.start(opt.mexp * 100)) {
            GF2X irreducible = ars.getIrreducibleFactor();
            GF2X characteristic = ars.getCharacteristicPolynomial();
            //os << "deg irreducible = " << dec << deg(irreducible) << endl;
            //os << "deg characteristic = " << dec << deg(characteristic)
            //     << endl;
            //os << "deg quotient = " << dec << deg(quotient) << endl;
            if (deg(irreducible) != opt.mexp) {
                os << "error" << endl;
                return -1;
            }
            getLCMPoly(characteristic, g);
            GF2X quotient = characteristic / irreducible;
            w128_t fixpoint = calc_fixpoint(g, irreducible, quotient);
            g.setFixPoint(fixpoint);
            cp.searchParity(g, irreducible);
            w128_t seed = {{1, 0, 0, 0}};
            g.seed(seed);
            if (!anni(g)) {
                return -1;
            }
            annihilate<w128_t>(&g, quotient);
            int veq52[52];
            DSFMTInfo info;
            info.bitSize = 128;
            info.elementNo = 2;
            int delta52
                = calc_dSFMT_equidistribution<w128_t, dSFMT>(g, veq52, 52, info,
                                                           opt.mexp);
            os << g.getParamString();
            os << dec << delta52 << endl;
            i++;
        } else {
            os << "search failed" << endl;
            break;
        }
    }
    if (opt.verbose) {
        time_t t = time(NULL);
        os << "search end at " << ctime(&t) << endl;
    }
    return 0;
}

static void output_help(string& pgm);

/**
 * command line option parser
 * @param opt a structure to keep the result of parsing
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @param start default start value
 * @return command line options have error, or not
 */
bool parse_opt(options& opt, int argc, char **argv) {
    opt.verbose = false;
    opt.mexp = 0;
    opt.count = 1;
    opt.seed = (uint64_t)clock();
    opt.filename = "";
    opt.fixedSL1 = -1;
    opt.fixedPOS1 = -1;
    int c;
    bool error = false;
    string pgm = argv[0];
    static struct option longopts[] = {
        {"verbose", no_argument, NULL, 'v'},
        {"file", required_argument, NULL, 'f'},
        {"count", required_argument, NULL, 'c'},
        {"seed", required_argument, NULL, 's'},
        {"fixed-sl1", required_argument, NULL, 'x'},
        {"fixed-pos1", required_argument, NULL, 'X'},
        {NULL, 0, NULL, 0}};
    errno = 0;
    for (;;) {
        c = getopt_long(argc, argv, "vs:f:c:x:X:", longopts, NULL);
        if (error) {
            break;
        }
        if (c == -1) {
            break;
        }
        switch (c) {
        case 's':
            opt.seed = strtoull(optarg, NULL, 0);
            if (errno) {
                error = true;
                cerr << "seed must be a number" << endl;
            }
            break;
        case 'x':
            opt.fixedSL1 = strtoull(optarg, NULL, 0);
            if (errno) {
                error = true;
                cerr << "fixed sl1 must be a number" << endl;
            }
            break;
        case 'X':
            opt.fixedPOS1 = strtoull(optarg, NULL, 0);
            if (errno) {
                error = true;
                cerr << "fixed pos1 must be a number" << endl;
            }
            break;
        case 'v':
            opt.verbose = true;
            break;
        case 'f':
            opt.filename = optarg;
            break;
        case 'c':
            opt.count = strtoll(optarg, NULL, 10);
            if (errno) {
                error = true;
                cerr << "count must be a number" << endl;
            }
            break;
        case '?':
        default:
            error = true;
            break;
        }
    }
    argc -= optind;
    argv += optind;
    if (argc < 1) {
        error = true;
    } else {
        long mexp = strtol(argv[0], NULL, 10);
        static const int allowed_mexp[] = {521, 1279,
                                           2203, 4253,
                                           11213, 19937,
                                           44497, -1};
        if (! errno) {
            bool found = false;
            for (int i = 0; allowed_mexp[i] > 0; i++) {
                if (mexp == allowed_mexp[i]) {
                    found = true;
                    break;
                }
            }
            if (! found) {
                error = true;
            }
        }
        if (errno || error){
            error = true;
            cerr << "mexp must be one of ";
            for (int i = 0; allowed_mexp[i] > 0; i++) {
                cerr << dec << allowed_mexp[i] << " ";
            }
            cerr << endl;
        }
        opt.mexp = mexp;
    }
#if 0
    if (!opt.filename.empty()) {
        ofstream ofs(opt.filename.c_str());
        if (ofs) {
            ofs.close();
        } else {
            error = true;
            cerr << "can't open file:" << opt.filename << endl;
        }
    }
#endif
    if (error) {
        output_help(pgm);
        return false;
    }
    return true;
}

/**
 * showing help message
 * @param pgm program name
 */
static void output_help(string& pgm) {
    cerr << "usage:" << endl;
    cerr << pgm
         << " [-s seed] [-v] [-c count]"
         << " [-f outputfile]"
         << " mexp"
         << endl;
    static string help_string1 = "\n"
"--verbose, -v        Verbose mode. Output parameters, calculation time, etc.\n"
"--file, -f filename  Parameters are outputted to this file. without this\n"
"                     option, parameters are outputted to standard output.\n"
"--count, -c count    Output count. The number of parameters to be outputted.\n"
"--seed, -s seed      seed of randomness.\n"
"--fixed-sl1          fix the parameter sl1 to given value.\n"
"--fixed-pos1         fix the parameter pos1 to given value.\n"
"mexp                 mersenne exponent.\n"
        ;
    cerr << help_string1 << endl;
}
