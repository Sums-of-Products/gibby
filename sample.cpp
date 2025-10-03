#include "src/DAGsampler.hpp"
#include "src/post_counts.cpp"
#include <iomanip>

/*. Sample from the posterior distribution of DAGs using Gibby algorithm.  
     
    The program initializes a DAGsampler object with a specified dataset and parameters.
    It then performs a burn-in phase followed by a sampling phase, where it repeatedly applies
    the fast Gibby, REV, and MBR moves to sample DAGs from the posterior distribution. The pruning method is integrated.
     
    The scores of the sampled DAGs are written to an output file, and the posterior probabilities
    of edges are computed and saved to another output file. 
     
    Make sure to set the correct file paths for input datasets and output files before running the program.

    Parameters to set in the init function:
        fname		    string		    File name (including the directory) 
        fmode		    bool			Data point per line (true) or data point per column (false)
        params	        string		    A string of parameters:
        -a accuracy		                The number of significant bits in approximations (default 15)
        -d maxind		                The maximum indegree parameter (default -1, i.e., infinity)
        -e ess		                    The ESS parameter of BDeu (default 1.0)
        -p prior		                The structure prior: 0 uniform, 1 fair, 2 fair+, -w edge(w) (default 1)
        -s model		                Scoring model (BDeu = 1, default)
        -K max		                    Maximum number of parents per node (default 0 = unrestricted)
        -M mem		                    Amount of RAM available in GiB (default 16)                                                                                     
        -P mode		                    Score pruning: 0 no pruning, 1 top-down, 2 bottom-up (default 0); 
                                        if preceded by digit k > 0, then only up to min{k, maxind} parents
        -R seed		                    Seed to the random number generator (default 0)  
        -I filen		                Input file name for local scores in the jkl format
        -O filen		                Output file name for local scores in the jkl format 
        
        
    Examples:	
        init (“data/mushroom.dat”, true, “-s 1 -e 10.0 -d 3 -a 16 -p 1 -P 22 -K 64 -M 8”)
        init (“-S scores/mushroom.jkl -a 16 -M 8”)
*/



int main(){                   // SELECT CORRECT PATH AND PARAMETERS                                 DATASETS USED IN THE EXPERIMENTS
    DAGsampler ds;
    //ds.init("data/pigs1k.dat", true, "-s 1 -e 1.0 -d 2 -K 64  -a 15 -p -441 -R 1 -P 2 ");         //Pigs  (441 nodes)
    //ds.init("data/andes1k.dat", true, "-s 1 -e 1.0 -d 6 -K 64  -a 15 -p -223 -R 1 -P 2 ");        //Andes  (223 nodes)
    //ds.init("data/path1k.dat", true, "-s 1 -e 1.0 -d 5 -K 64  -a 15 -p -109 -R 1 -P 2 ");         //Pathfinder  (109 nodes)
    //ds.init("data/hail1k.dat", true, "-s 1 -e 1.0 -d 4 -K 64  -a 15 -p -56 -R 1 -P 2 ");          //Hailfinder  (56 nodes)
    //ds.init("data/alarm1K.csv", true, "-s 1 -e 1.0 -d 4 -a 10 -p -37 -R 1 -P 2");                 //Alarm  (37 nodes)
    //ds.init("data/child1k.dat", true, "-s 1 -e 1.0 -d 9 -a 10 -p -20 -R 1 -P 2");                 //Child  (20 nodes)
    //ds.init("data/zoo.dat", true, "-s 1 -e 1.0 -d 13 -a 10 -p -17 -R 1 -P 2");                    //Zoo  (17 nodes)
    //ds.init("data/sachs1k.dat", true, "-s 1 -e 1.0 -d 10 -a 10 -p -11 -R 1 -P 2");                //Sachs  (11 nodes) 
    ds.init("data/asia1k.dat", true, "-s 1 -e 1.0 -d 7 -a 10 -p -8 -R 1 -P 2");                     //Asia  (8 nodes) 
    uint32_t infomask = 0b11111001;

    int n_nodes = ds.getn(); 
    int burn_in_iter=100; // Select the number of burn-in iterations
    int iter = 10000; //  Select the number of iterations
    std::vector<std::vector<int>> adj(n_nodes, std::vector<int>(n_nodes, 0));

    // Burn-in phase
    for (int it = 0; it < burn_in_iter; it++) {
        ds.sGib(10000);    // Select the number of Gibby iterations per main iteration  
        ds.sREV(200);      // Select the number of REV iterations per main iteration 
        ds.sMBR_alt(200, false);  // Select the number of MBR iterations per main iteration
    }
    
    // Sampling phase
    std::ofstream outfile("results/scores.txt"); // Output file for the scores of the sampled DAGs. SELECT CORRECT PATH
    for (int it = 0; it < iter; it++) {
        ds.sGib(10000);     
        ds.sREV(200);     
        ds.sMBR_alt(200, false);
        double t = 0.0;
        std::vector<int> P;
        for (int i = 0; i < n_nodes; ++i) {
            P = ds.getp(i);
            t += ds.gets(i, P);
            }
        outfile << std::fixed << std::setprecision(4) << t << std::endl;
        add_post(ds, n_nodes, adj);     
    }
    outfile.close();
    save_matrix_to_file(adj, iter, "results/posterior_prob_edges.txt"); // Output file for the posterior probability of the edges. SELECT CORRECT PATH
    std::cout << "sample complete" << std::endl;
    return 0;
}
