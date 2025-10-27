#include "src/DAGsampler.hpp"
#include "src/post_counts.cpp"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <filesystem>
#include <chrono>
#include <cmath>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <datafile or -I file> [options]\n";
        return 1;
    }

    bool shortcut_mode = false;
    std::string datafile;


    int iter = -1; 
    int burn_in_iter = -1; 
    int gibby_iter = 10000;
    int rev_iter = 200;
    int mbr_iter = 200;

    double ess = 1.0;
    int sig_bits = 10;
    int max_indegree = -1;
    int max_parents = 0;
    int structure_prior = 1;
    int pruning = 2; 
    int s_flag = 1;
    int M_param = 16; 

    std::string parent_scores_file = "";
    int seed_value = 0;

    if (std::string(argv[1]) == "-I") {
        shortcut_mode = true;
        if (argc < 3) {
            std::cerr << "Error: -I requires a file\n";
            return 1;
        }
        datafile = argv[2];

        for (int i = 3; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "-iter" && i + 1 < argc) iter = std::atoi(argv[++i]);
            else if (arg == "-burnin" && i + 1 < argc) burn_in_iter = std::atoi(argv[++i]);
            else if (arg == "-FBM" && i + 1 < argc) gibby_iter = std::atoi(argv[++i]);
            else if (arg == "-REV" && i + 1 < argc) rev_iter = std::atoi(argv[++i]);
            else if (arg == "-MBR" && i + 1 < argc) mbr_iter = std::atoi(argv[++i]);
            else if (arg == "-a" && i + 1 < argc) sig_bits = std::atoi(argv[++i]);
            else if (arg == "-M" && i + 1 < argc) M_param = std::atoi(argv[++i]);
            else if (arg == "-d" && i + 1 < argc) max_indegree = std::atoi(argv[++i]); 
            else if (arg[0] == '-') i++; 
        }

        if (iter <= 0) iter = 10000;  // default

    } else {

        datafile = argv[1];
        for (int i = 2; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "-iter" && i + 1 < argc) iter = std::atoi(argv[++i]);
            else if (arg == "-burnin" && i + 1 < argc) burn_in_iter = std::atoi(argv[++i]);
            else if (arg == "-FBM" && i + 1 < argc) gibby_iter = std::atoi(argv[++i]);
            else if (arg == "-REV" && i + 1 < argc) rev_iter = std::atoi(argv[++i]);
            else if (arg == "-MBR" && i + 1 < argc) mbr_iter = std::atoi(argv[++i]);
            else if (arg == "-O" && i + 1 < argc) parent_scores_file = argv[++i];
            else if (arg == "-e" && i + 1 < argc) ess = std::stod(argv[++i]);
            else if (arg == "-s" && i + 1 < argc) s_flag = std::atoi(argv[++i]);
            else if (arg == "-a" && i + 1 < argc) sig_bits = std::atoi(argv[++i]);
            else if (arg == "-d" && i + 1 < argc) max_indegree = std::atoi(argv[++i]);
            else if (arg == "-K" && i + 1 < argc) max_parents = std::atoi(argv[++i]);
            else if (arg == "-M" && i + 1 < argc) M_param = std::atoi(argv[++i]);
            else if (arg == "-p" && i + 1 < argc) structure_prior = std::atoi(argv[++i]);
            else if (arg == "-P" && i + 1 < argc) pruning = std::atoi(argv[++i]);
            else if (arg == "-R" && i + 1 < argc) seed_value = std::atoi(argv[++i]);
        }
        if (iter <= 0) iter = 10000;
    }

    if (burn_in_iter <= 0) {
        burn_in_iter = std::max(1, iter / 10);
    }


    if (seed_value == 0)
        seed_value = static_cast<int>(
            std::chrono::system_clock::now().time_since_epoch().count() % 1000000);


    std::string user_params;
    if (shortcut_mode) {
        user_params = "-I " + datafile + " ";
        user_params += "-a " + std::to_string(sig_bits) + " ";
        user_params += "-M " + std::to_string(M_param) + " ";
    } else {
        user_params = "-s " + std::to_string(s_flag) + " ";
        user_params += "-e " + std::to_string(ess) + " ";
        user_params += "-a " + std::to_string(sig_bits) + " ";
        user_params += "-P " + std::to_string(pruning) + " ";
        user_params += "-K " + std::to_string(max_parents) + " ";
        user_params += "-M " + std::to_string(M_param) + " ";
        user_params += "-R " + std::to_string(seed_value) + " ";
        if (!parent_scores_file.empty())
            user_params += "-O " + parent_scores_file + " ";
    }


    std::filesystem::create_directories("results");

    DAGsampler ds;
    if (shortcut_mode) {
        ds.init(user_params);
    } else {
        ds.init(datafile, true, user_params);
    }

    int n_nodes = ds.getn();
    std::vector<std::vector<int>> adj(n_nodes, std::vector<int>(n_nodes, 0));

    std::string base_name = datafile.substr(datafile.find_last_of("/\\") + 1);
    std::string score_output = "results/" + base_name + "_R" + std::to_string(seed_value) + "_scores.txt";
    std::string posterior_output = "results/" + base_name + "_R" + std::to_string(seed_value) + "_posterior.txt";


    if (shortcut_mode) {
        std::cout << " \n";
        std::cout << " \n";
        std::cout << " \n";
        std::cout << "█▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀█\n";
        std::cout << "█                              Gibby DAGs Sampler                               █\n";
        std::cout << "█▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄█\n";
        std::cout << " \n";
        std::cout << "╭───────────────────────────── Scoring parameters ───────────────────────────────\n";
        std::cout << "│ \n";
        std::cout << "│ Score file: " << datafile << "\n";
        if (max_indegree > 0) {
        std::cout << "│ Maximum in-degree: " << max_indegree << "\n";
        }
        std::cout << "│ Number of candidate parents per node: "
                  << (max_parents == 0 ? "unrestricted" : std::to_string(max_parents)) << "\n";
        std::cout << "│ Number of significant bits in approximations: " << sig_bits << "\n";
        std::cout << "│ Amount of RAM available (GiB): " << M_param << "\n";
        std::cout << "│ \n";
    } else {
        std::string prior_str = (structure_prior == 0) ? "uniform" :
                                (structure_prior == 1) ? "fair" :
                                (structure_prior == 2) ? "fair+" :
                                "edge(" + std::to_string(std::abs(structure_prior)) + ")";
        std::string pruning_str = (pruning == 0) ? "no pruning" :
                                  (pruning == 1) ? "top-down" :
                                  (pruning == 2) ? "bottom-up" : "unknown";
        std::cout << " \n";
        std::cout << " \n";
        std::cout << " \n";
        std::cout << "█▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀█\n";
        std::cout << "█                              Gibby DAGs Sampler                               █\n";
        std::cout << "█▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄█\n";
        std::cout << " \n";
        std::cout << "╭───────────────────────────── Scoring parameters ───────────────────────────────\n";
        std::cout << "│ \n";
        std::cout << "│ Dataset: " << datafile << "\n";
        std::cout << "│ Maximum in-degree: " << (max_indegree == -1 ? n_nodes - 1 : max_indegree) << "\n";
        std::cout << "│ Number of candidate parents per node: "
                  << (max_parents == 0 ? "unrestricted" : std::to_string(max_parents)) << "\n";
        std::cout << "│ Structure prior: " << prior_str << "\n";
        std::cout << "│ Pruning: " << pruning_str << "\n";
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "│ Score: BDeu (ESS=" << ess << ")\n";
        std::cout.unsetf(std::ios::floatfield);
        std::cout << "│ Number of significant bits in approximations: " << sig_bits << "\n";
        std::cout << "│ Amount of RAM available (GiB): " << M_param << "\n";
        std::cout << "│ \n";
    }

    std::cout << "╭────────────────────────────── MCMC parameters ─────────────────────────────────\n";
    std::cout << "│ \n";
    std::cout << "│ Burn-in iterations: " << burn_in_iter << "\n";
    std::cout << "│ Main iterations: " << iter << "\n";
    std::cout << "│ Fast basic moves per iteration: " << gibby_iter << "\n";
    std::cout << "│ REV moves per iteration: " << rev_iter << "\n";
    std::cout << "│ MBR moves per iteration: " << mbr_iter << "\n";
    std::cout << "│ Random seed (R): " << seed_value << "\n";
    std::cout << "│ \n";
    std::cout << "╭─────────────────────────────── Output files ───────────────────────────────────\n";
    std::cout << "│ \n";
    std::cout << "│ Sampled DAGs scores -> " << score_output << "\n";
    std::cout << "│ Edge probability matrix -> " << posterior_output << "\n";
    if (!parent_scores_file.empty())
        std::cout << "│ Parent sets scores -> " << parent_scores_file << "\n";
    std::cout << "│ \n";
    //std::cout << "---------------------------------------------------------------------------------\n";
    std::cout << "\n";
    std::cout << "♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ Sampling phase ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦ ♦\n";
    std::cout << "\n";

    std::cout << " Starting burn-in phase..." << std::endl;
    int progress_interval1 = 100;
    for (int it = 0; it < burn_in_iter; it++) {
        ds.sGib(gibby_iter);
        ds.sREV(rev_iter);
        ds.sMBR_alt(mbr_iter, false);

        if ((it + 1) % progress_interval1 == 0 || it == iter - 1)
            std::cout << "\r Burn-in iteration " << (it + 1) << "/" << burn_in_iter << "" << std::flush;
    }
    std::cout << std::endl;
    std::cout << std::endl;


    std::ofstream outfile(score_output);
    if (!outfile) {
        std::cerr << " Error: cannot open output file for scores: " << score_output << std::endl;
        return 1;
    }

    std::cout << " Starting main sampling phase..." << std::endl;
    int progress_interval = 100;
    for (int it = 0; it < iter; it++) {
        ds.sGib(gibby_iter);
        ds.sREV(rev_iter);
        ds.sMBR_alt(mbr_iter, false);

        double t = 0.0;
        std::vector<int> P;
        for (int i = 0; i < n_nodes; ++i) {
            P = ds.getp(i);
            t += ds.gets(i, P);
        }
        outfile << std::fixed << std::setprecision(4) << t << std::endl;
        add_post(ds, n_nodes, adj);

        if ((it + 1) % progress_interval == 0 || it == iter - 1)
            std::cout << "\r Iteration " << (it + 1) << "/" << iter << "" << std::flush;
    }

    std::cout << std::endl;
    outfile.close();
    save_matrix_to_file(adj, iter, posterior_output);

    std::cout << " Sampling complete.\n";
    std::cout << "\n";

    return 0;
}
