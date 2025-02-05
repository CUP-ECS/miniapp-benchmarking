#include <iostream>
#include <memory>
#include <type_traits>
#include <cstdlib>
#include <Cabana_Core.hpp>
#include <Kokkos_Core.hpp>

#include <string>

#include <algorithm>
#include <tclap/CmdLine.h>

#include <nlohmann/json.hpp>
#include <cstdlib>

#include <cmath>
#include <filesystem>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>




enum distribution {
    GAUSSIAN,
    EMPIRICAL
 };

typedef enum distribution distribution_t;

enum prefix {
    A,B,K,M,G
};

typedef enum prefix prefix_t;

static int typesize = 8;
static int numpes = 0;
static int nsamples = 25;
static int niterations = 100;
static int nneighbors = -1;
static int nneighbors_stdv = -1;
static int nowned = -1;
static int nowned_stdv = -1;
static int nremote = -1;
static int nremote_stdv = -1;
static int blocksz = -1;
static int blocksz_stdv = -1;
static int stride = -1;
static int stride_stdv = -1;
static int unit_div = 1;
static prefix unit_symbol = A;
static std::string filepath = "";
static distribution_t distribution_type = GAUSSIAN;
static bool irregularity = 1;
static bool irregularity_owned = 1;
static bool irregularity_neighbors = 1;
static bool irregularity_stride = 1;
static bool irregularity_blocksz = 1;
static bool irregularity_remote = 1;
static bool report_params = 0;
static int seed = -1;


int gauss_dist(double mean, double stdev) {
    // generates two random numbers that form the seeds
    // of the transform
    double u1, u2, r, theta;
    int generated = -1;

    while (generated <= 0) {
        u1 = (double)rand() / RAND_MAX;
        u2 = (double)rand() / RAND_MAX;

        // generates the R and Theta values from the above
        // documentation
        r = sqrt(-2.*log(u1));
        theta = (2*M_PI*u2);

        // an additional number can be generated in the
        // same distribution using the alternate form
        // ((r*sin(theta)) * stdev) + mean

        generated = round(((r*cos(theta)) * stdev ) + mean);
    }

    return generated;
}




void parse_config_file() {
  ifstream input_file("benchmark.json");
  json j;
  input_file >> j;

  // Accessing the data
  for (const auto& param : j["parameters"]) {
    cout << "PARAM: " << param["name"] << endl;
    cout << "BIN_COUNT: " << param["bin_count"] << endl;
    cout << "MIN: " << param["min"] << endl;
    cout << "MAX: " << param["max"] << endl;
    cout << "MEAN: " << param["mean"] << endl;
    cout << "STDEV: " << param["stdev"] << endl;

    // Iterate through bins
    for (const auto& bin : param["bins"]) {
      cout << "BIN_MIN: " << bin["bin_min"] << ", ";
      cout << "BIN_MAX: " << bin["bin_max"] << ", ";
      cout << "BIN_PROP: " << bin["bin_prop"] << ", ";
      cout << "BIN_MEAN: " << bin["bin_mean"] << ", ";
      cout << "BIN_STDEV: " << bin["bin_stdev"] << endl;
    }

  }



//    char* params[] = { "nowned",
//                       "nremote",
//                       "blocksize",
//                       "stride",
//                       "comm_partners" };
//
//    for (int index = 0; index < (sizeof(params) / sizeof(params[0])); index++) {
//        char* param = params[index];
//
//        char param_key[25] = "PARAM: ";
//        strcat(param_key, param);
//
//        FILE* fp = fopen(filepath, "r");
//        if (fp == NULL) {
//            fprintf(stderr, "Error: Unable to open file %s\n", filepath);
//            exit(1);
//        }
//
//        char* line = NULL;
//        size_t linecap = 0;
//        ssize_t linelen;
//
//        // iterates through lines in the file until we get to the line
//        // that contains the parameter we're trying to generate for
//        while ((linelen = getline(&line, &linecap, fp)) != -1) {
//            // remove trailing new line to better do comparison
//            if (linelen > 0 && line[linelen-1] == '\n') {
//                line[linelen-1] = '\0';
//            }
//
//            // if the correct parameter is found,
//            // stop iterating through the file
//            if (strcmp(param_key, line) == 0) {
//                break;
//            }
//        }
//
//        // iterates through non-bin data points
//        for (int i = 0; i < 5; i++) {
//            // gets the next line which contains the BIN_COUNT data
//            linelen = getline(&line, &linecap, fp);
//
//            // remove trailing new line to better do comparison
//            if (linelen > 0 && line[linelen-1] == '\n') {
//                line[linelen-1] = '\0';
//            }
//
//            // gets the name of the data point being read
//            char* token = strtok(line, ":");
//
//            // writes the correct token value to the correct variable from file input
//            if (strcmp(token, "MEAN") == 0) {
//                // puts the mean value in the correct parameter spot
//                // based on current parameter choice
//                if (strcmp(param, "nowned") == 0) {
//                    nowned = atoi(strtok(NULL, " "));
//                } else if (strcmp(param, "nremote") == 0) {
//                    nremote = atoi(strtok(NULL, " "));
//                } else if (strcmp(param, "blocksize") == 0) {
//                    blocksz = atoi(strtok(NULL, " "));
//                } else if (strcmp(param, "stride") == 0) {
//                    stride = atoi(strtok(NULL, " "));
//                } else if (strcmp(param, "comm_partners") == 0) {
//                    nneighbors = atoi(strtok(NULL, " "));
//                }
//            } else if (strcmp(token, "STDEV") == 0) {
//                // puts the stdev value in the correct parameter spot
//                // based on current parameter choice
//                if (strcmp(param, "nowned") == 0) {
//                    nowned_stdv = atoi(strtok(NULL, " "));
//                } else if (strcmp(param, "nremote") == 0) {
//                    nremote_stdv = atoi(strtok(NULL, " "));
//                } else if (strcmp(param, "blocksize") == 0) {
//                    blocksz_stdv = atoi(strtok(NULL, " "));
//                } else if (strcmp(param, "stride") == 0) {
//                    stride_stdv = atoi(strtok(NULL, " "));
//                } else if (strcmp(param, "comm_partners") == 0) {
//                    nneighbors_stdv = atoi(strtok(NULL, " "));
//                }
//            }
//        }
//    }
}

bool setValue(){
  return true;
}

void exitError(const std::string& error_message) {
    std::cerr << error_message << std::endl; // Use std::cerr for error messages
    std::exit(-1); // Exit the program with error code
}


void setAndCheckValue(int& value, TCLAP::ValueArg<int>& arg, const char* errorMessage, int minValue = 0, int maxValue = INT_MAX) {
  int tempValue = arg.getValue();

  if (tempValue != -1) {
    value = tempValue;
  }

  // General value check
  if (value < minValue || value > maxValue) {
    exitError(errorMessage);
  }
}


void parseArgs(int argc, char **argv){

  try {
        TCLAP::CmdLine cmd("\nNOTE: Setting parameters for the benchmark such as (neighbors, owned, remote, blocksize, and stride)"
                           "sets parameters to those values for the reference benchmark."
                           "Those parameters are then randomized for the irregular samples"
                           "where the user-set parameters become averages for the random generation."
                           "Use the `--disable-irregularity` flag to only run the reference benchmark.", ' ', "1.0");

        TCLAP::ValueArg<std::string> filepathArg("f", "filepath", "Path to the BENCHMARK_CONFIG file", false, "NOFILE", "string");
        TCLAP::ValueArg<int> typeSizeArg("t", "typesize", "Size of the variable being sent (in bytes)", false, 8, "int");
        TCLAP::ValueArg<int> samplesArg("I", "samples", "Number of random samples to generate", false, 25, "int");
        TCLAP::ValueArg<int> iterationsArg("i", "iterations", "Number of updates each sample performs", false, 100, "int");
        TCLAP::ValueArg<int> neighborsArg("n", "neighbors", "Average number of neighbors each process communicates with", false, -1, "int");
        TCLAP::ValueArg<int> neighborsStdvArg("N", "neighbors_stdv", "Standard deviation of the number of neighbors each process communicates with", false, -1, "int");
        TCLAP::ValueArg<int> ownedAvgArg("o", "owned_avg", "Average byte count for data owned per node", false, -1, "int");
        TCLAP::ValueArg<int> ownedStdvArg("O", "owned_stdv", "Standard deviation byte count for data owned per node", false, -1, "int");
        TCLAP::ValueArg<int> remoteAvgArg("r", "remote_avg", "Average amount of data each process receives", false, -1, "int");
        TCLAP::ValueArg<int> remoteStdvArg("R", "remote_stdv", "Standard deviation of the amount of data each process receives", false, -1, "int");
        TCLAP::ValueArg<int> blockSizeAvgArg("b", "blocksize_avg", "Average size of transmitted blocks", false, -1, "int");
        TCLAP::ValueArg<int> blockSizeStdvArg("B", "blocksize_stdv", "Standard deviation of transmitted block sizes", false, -1, "int");
        TCLAP::ValueArg<int> strideArg("s", "stride", "Average size of stride", false, -1, "int");
        TCLAP::ValueArg<int> strideStdvArg("T", "stride_stdv", "Standard deviation of stride", false, -1, "int");
        TCLAP::ValueArg<int> seedArg("S", "seed", "Positive integer to be used as seed for random number generation", false, -1, "int");
        TCLAP::ValueArg<std::string> distributionArg("d", "distribution", "Choose from: gaussian (default), empirical", false, "gaussian", "string");
        TCLAP::ValueArg<std::string> unitsArg("u", "units", "Choose from: a,b,k,m,g (auto, bytes, kilobytes, etc.)", false, "auto", "string");


        TCLAP::SwitchArg reportParamsArg("", "report-params", "Enables parameter reporting for use with analysis scripts", false);
        TCLAP::SwitchArg disableirregularityArg("", "disable-irregularity", "Use the `--disable-irregularity` flag to only run the reference benchmark.", false);
        cmd.add(filepathArg);
        cmd.add(typeSizeArg);
        cmd.add(samplesArg);
        cmd.add(iterationsArg);
        cmd.add(neighborsArg);
        cmd.add(neighborsStdvArg);
        cmd.add(ownedAvgArg);
        cmd.add(ownedStdvArg);
        cmd.add(remoteAvgArg);
        cmd.add(remoteStdvArg);
        cmd.add(blockSizeAvgArg);
        cmd.add(blockSizeStdvArg);
        cmd.add(strideArg);
        cmd.add(strideStdvArg);
        cmd.add(seedArg);
        cmd.add(distributionArg);
        cmd.add(unitsArg);
        cmd.add(reportParamsArg);
        cmd.add(disableirregularityArg);
        cmd.parse(argc, argv);


        filepath = filepathArg.getValue();
        bool config_file_used = false;



        if (filepath != "NOFILE") {
            if(std::filesystem::exists(filepath)) {
                config_file_used = true;
                parse_config_file();

            }else{
                exitError("ERROR: the specified filepath doesn't exist, exiting...");
            }
        }




    setAndCheckValue(typesize, typeSizeArg, "ERROR: Invalid typesize\n", 1, 8);

    // For nsamples, no specific range, only non-negative check
    setAndCheckValue(nsamples, samplesArg, "ERROR: Invalid number of samples\n", 0);

    // For niterations, same non-negative check
    setAndCheckValue(niterations, iterationsArg, "ERROR: Invalid number of iterations\n", 0);

    // For nneighbors, same non-negative check
    setAndCheckValue(nneighbors, neighborsArg, "ERROR: Invalid number of neighbors\n", 0);

    // For nneighbors_stdv, no specific range, only non-negative check
    setAndCheckValue(nneighbors_stdv, neighborsStdvArg, "ERROR: Invalid neighbors std\n", 0);

    // For nowned, no specific range, only non-negative check
    setAndCheckValue(nowned, ownedAvgArg, "ERROR: Invalid number of owned\n", 0);

    // For nowned_stdv, no specific range, only non-negative check
    setAndCheckValue(nowned_stdv, ownedStdvArg, "ERROR: Invalid owned std\n", 0);

    // For nremote, no specific range, only non-negative check
    setAndCheckValue(nremote, remoteAvgArg, "ERROR: Invalid number of remote\n", 0);

    // For nremote_stdv, no specific range, only non-negative check
    setAndCheckValue(nremote_stdv, remoteStdvArg, "ERROR: Invalid remote std\n", 0);

    // For blocksz, no specific range, only non-negative check
    setAndCheckValue(blocksz, blockSizeAvgArg, "ERROR: Invalid block size\n", 0);

    // For blocksz_stdv, no specific range, only non-negative check
    setAndCheckValue(blocksz_stdv, blockSizeStdvArg, "ERROR: Invalid block size std\n", 0);

    // For stride, no specific range, only non-negative check
    setAndCheckValue(stride, strideArg, "ERROR: Invalid stride\n", 0);

    // For stride_stdv, no specific range, only non-negative check
    setAndCheckValue(stride_stdv, strideStdvArg, "ERROR: Invalid stride std\n", 0);



        std::string unit = unitsArg.getValue();

        std::unordered_map<std::string, std::pair<char, int>> unit_map = {
            {"auto", {A, 1}},
            {"a", {A, 1}},
            {"bytes", {A, 1}},
            {"b", {A, 1}},
            {"kilobytes", {K, 1024}},
            {"k", {K, 1024}},
            {"megabytes", {M, 1024 * 1024}},
            {"m", {M, 1024 * 1024}},
            {"gigabytes", {G, 1024 * 1024 * 1024}},
            {"g", {G, 1024 * 1024 * 1024}}
        };




        auto it = unit_map.find(unit);
        if (it != unit_map.end()) {
            unit_symbol = it->second.first;
            unit_div = it->second.second;
        } else {
            exitError("ERROR: Invalid formatting choice [b, k, m, g]");
        }


        std::string distribution=distributionArg.getValue();

        if (distribution == "gaussian"  ||
            distribution == "g") {
          distribution_type = GAUSSIAN;

        }else if (distribution == "empirical"
                  || distribution == "e") {
          distribution_type=EMPIRICAL;
        }else {
          exitError("ERROR: Invalid distribution choice [empirical,gaussian]\n");
        }

       int seedholder=seedArg.getValue();
        if(seed!=-1&&seedholder==-1) {
          seed=time(NULL);
        }



        srand(seed);

        irregularity=disableirregularityArg.getValue();


        //        irregularity_owned
        //        irregularity_neighbors
        //        irregularity_stride
        //        irregularity_blocksz
        //        irregularity_remote
        //        report_params
















    }
    catch (TCLAP::ArgException &e) {
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
        exit(-1);
    }
}

int main(int argc, char** argv)
{

    parseArgs(argc, argv);



	printf("Hi\n");
	return 0;
}
