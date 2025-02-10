#include <iostream>
#include <memory>
#include <type_traits>
#include <cstdlib>
#include <Cabana_Core.hpp>
#include <Kokkos_Core.hpp>
#include <fstream>

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
using json = nlohmann::json;



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


int empirical_dist(const char* param) {
//    // set the parameter string to look for.
    char param_key[25] = "PARAM: ";
    strcat(param_key, param);
//
//    // attempts to open the file pointer
//    // and ensure that the file can be read
//    // exits if fails
//    FILE* fp = fopen(filepath, "r");
//    if (fp == NULL) {
//        fprintf(stderr, "Error: Unable to open file %s\n", filepath);
//        exit(1);
//    }
//
//    // stores the current read line and
//    // some auxiliary information that assists in parsing the file
//    char* line = NULL;
//    size_t linecap = 0;
//    ssize_t linelen;
//
//    // iterates through lines in the file until we get to the line
//    // that contains the parameter we're trying to generate for
//    while ((linelen = getline(&line, &linecap, fp)) != -1) {
//        // remove trailing new line to better do comparison
//        if (linelen > 0 && line[linelen-1] == '\n') {
//            line[linelen-1] = '\0';
//        }
//
//        // if the correct parameter is found,
//        // stop iterating through the file
//        if (strcmp(param_key, line) == 0) {
//            break;
//        }
//    }
//
//    int binCount;
//    int dataMin;
//    int dataMax;
//    int dataMean;
//    int dataStdev;
//
//    // Iterates through non-bin data points
//    // and collect information about the data and number of bins
//    for (int i = 0; i < 5; i++) {
//        // gets the next line which contains the BIN_COUNT data
//        linelen = getline(&line, &linecap, fp);
//
//        // remove trailing new line to better do comparison
//        if (linelen > 0 && line[linelen-1] == '\n') {
//            line[linelen-1] = '\0';
//        }
//
//        // gets the name of the data point being read
//        char* token = strtok(line, ":");
//
//        // writes the correct token value to the correct variable from file input
//        if (strcmp(token, "BIN_COUNT") == 0) {
//            binCount = atoi(strtok(NULL, " "));
//        } else if (strcmp(token, "MIN") == 0) {
//            dataMin = atoi(strtok(NULL, " "));
//        } else if (strcmp(token, "MAX") == 0) {
//            dataMax = atoi(strtok(NULL, " "));
//        } else if (strcmp(token, "MEAN") == 0) {
//            dataMean = atoi(strtok(NULL, " "));
//        } else if (strcmp(token, "STDEV") == 0) {
//            dataStdev = atoi(strtok(NULL, " "));
//        } else {
//            printf("Error: reading data file failed, invalid data point found: %s\n", token);
//            exit(1);
//        }
//    }
//
//    // do a quick sanity check to ensure that the data doesn't only
//    // fall into a single bin
//    if (dataMin == dataMax) {
//        // close file as reading is no longer required
//        fclose(fp);
//
//        return dataMin; // could also return data max, it is arbitrary
//    }
//
//    // chose a uniformly distributed value between (0,1)
//    double binSelection = (double)rand() / RAND_MAX;
//
//
//    // used to track if we have checked a bin
//    // if we have/are, the bin's proportion of the total dataset
//    // is added to the binSum.
//    // If the binSelection falls under the binSum, it is said
//    // to be in the bin and we can use that bin's data to generate values
//    double binSum = 0.0;
//
//    // iterates through the bins for a parameter to identity
//    // which one should be chosen based on the random binSelection
//    // value chosen above
//    for (int i = 0; i < binCount; i++) {
//        // gets the next line which contains the BIN_COUNT data
//        linelen = getline(&line, &linecap, fp);
//
//        // remove trailing new line to better do comparison
//        if (linelen > 0 && line[linelen-1] == '\n') {
//            line[linelen-1] = '\0';
//        }
//
//        // stores data from the current bin being checked
//        int binMin;
//        int binMax;
//        int binMean;
//        int binStdev;
//
//        // gets the bin minimum data point
//        char* token = strtok(line, ",");
//        binMin = atoi(token);
//
//        // gets the bin maximum data point
//        token = strtok(NULL, ",");
//        binMax = atoi(token);
//
//        // gets the binProp value and adds it to the binSum
//        // NOTE: I don't think this will ever chose a bin
//        //       that has a prop value of 0, as binSum would remain
//        //       the same and would not trigger the binSelection check
//        token = strtok(NULL, ",");
//        binSum += strtod(token, NULL);
//
//        // if so, then the selection is in the most
//        // recently searched bin, so we can proceed
//        // to generating a random value
//        // also triggers if on the last bin and no value has been selected
//        if (binSelection < binSum || i == (binCount-1)) {
//
//            // if the bin is the correct one,
//            // get the mean of the data points
//            // for this bin
//            token = strtok(NULL, ",");
//            binMean = atoi(token);
//
//            // if the bin is the correct one,
//            // get the standard deviation of the
//            // data points in this bin
//            token = strtok(NULL, ",");
//            binStdev = atoi(token);
//
//            // close file as reading is no longer required
//            fclose(fp);
//
//            // return a normally distributed value with the
//            // mean and standard deviation from the
//            // chosen bin
//            return gauss_dist(binMean, binStdev);
//        }
//
//    }
//
//    // if we have somehow made it to this point
//    // without generating a value, we should throw an error
//    printf("Error: empirical value could not be generated due to an unknown error.\n");
//    printf("Parameter type being generated: %s\n", param);
//    printf("Config file which could have caused this error: %s\n", filepath);
//    printf("If you encounter a bug which causes this for any reason, please open an issue on Github\n");
//    printf("Here is additional debugging information\n");
//    printf("\tbinSelection: %f\n", binSelection);
//    printf("\tbinSum: %f\n", binSum);
//    exit(1);
    return -1;
}


int benchmark(int penum) {
    // for benchmarks with irregularity disabled,
    // having more than 1 sample is not useless
    if (irregularity == 0) {
        if (penum == 0) {
            printf("Irregularity has been disabled, setting nsamples to 1...\n");
            printf("To increase iterations for benchmarks with irregularity disabled, please specify parameters via the \"-i [iterations]\"\n");
        }
        nsamples = 1;
    }

    // stores original parameter values if set by the user
    // this is because the original values change with irregularity enabled
    // and the new starting points for each iteration need to be set
    int nowned_orig = nowned;
    int nneighbors_orig = nneighbors;
    int nremote_orig = nremote;
    int blocksz_orig = blocksz;
    int stride_orig = stride;

    for(int sample_iter = 0; sample_iter < nsamples+1; sample_iter++) {
        // if irregularity is enabled, perform a reference benchmark first
        if (irregularity) {
            if (sample_iter == 0) {
                // if this is the first iteration, print the reference benchmark information
                if (penum == 0) {
                    printf("Non-Irregular Reference Benchmark\n");
                }
            } else {
                /*
                // if irregularity is enabled, set random parameters
                // the following code determines if the default values should
                // be used as random generation starting points
                // if a value is set by the user, it makes sense to use that
                // parameter as a starting point for all iterations values
                // if a value is unset, the benchmark's default values are used
                // as the starting point (mean, stdev) for the gaussian dist
                */
                if (irregularity_owned) {
                    switch (distribution_type) {
                        case GAUSSIAN:
                            nowned = gauss_dist(nowned_orig, nowned_stdv);
                            break;
                        case EMPIRICAL:
                            nowned = empirical_dist("nowned");
                            break;
                        default:
                            nowned = nowned_orig;
                            break;
                    }
                }
                if (irregularity_remote) {
                    switch (distribution_type) {
                        case GAUSSIAN:
                            nremote = gauss_dist(nremote_orig, nremote_orig);
                            break;
                        case EMPIRICAL:
                            nremote = empirical_dist("nremote");
                            break;
                        default:
                            nremote = nremote_orig;
                            break;
                    }
                }
                if (irregularity_blocksz) {
                    switch (distribution_type) {
                        case GAUSSIAN:
                            blocksz = gauss_dist(blocksz_orig, blocksz_stdv);
                            break;
                        case EMPIRICAL:
                            blocksz = empirical_dist("blocksize");
                            break;
                        default:
                            blocksz = blocksz_orig;
                            break;
                    }
                }
                if (irregularity_neighbors) {
                    switch (distribution_type) {
                        case GAUSSIAN:
                            nneighbors = gauss_dist(nneighbors_orig, nneighbors_stdv);
                            break;
                        case EMPIRICAL:
                            nneighbors = empirical_dist("comm_partners");
                            break;
                        default:
                            nneighbors = nneighbors_orig;
                            break;
                    }
                }
                if (irregularity_stride) {
                    switch (distribution_type) {
                        case GAUSSIAN:
                            stride = gauss_dist(stride_orig, stride_stdv);
                            break;
                        case EMPIRICAL:
                            stride = empirical_dist("stride");
                            break;
                        default:
                            stride = stride_orig;
                            break;
                    }
                }

                if (report_params) {
                    printf("PARAM: nowned - %d\n", nowned);
                    printf("PARAM: nremote - %d\n", nremote);
                    printf("PARAM: blocksize - %d\n", blocksz);
                    printf("PARAM: stride - %d\n", stride);
                    printf("PARAM: nneighbors - %d\n", nneighbors);
                }

                // print benchmark status each iteration
                if (penum == 0) {
                    printf("Benchmark Iteration: %d/%d\n", sample_iter, nsamples);
                }
            }
        } else {
            if (sample_iter == 0) {
                nsamples--;
            }
        }
    }








    return 0;
}

void migrationExample()
{
    /*
      The distributor is a communication plan allowing for the migration of
      data from one uniquely-owned distribution to another uniquely-owned
      distribution. Data migration may be applied to entire AoSoA data
      structures as well as slices.

      In this example we will demonstrate building a distributor communication
      plan and migrating data.

      Note: The distributor uses MPI for data migration. MPI is initialized
      and finalized in the main function below.

      Note: The distributor uses GPU-aware MPI communication. If AoSoA data is
      allocated in GPU memory, this feature will be used automatically.
    */

    std::cout << "Cabana Migration Example\n" << std::endl;

    /*
       Get parameters from the communicator. We will use MPI_COMM_WORLD for
       this example but any MPI communicator may be used.
    */
    int comm_rank = -1;
    MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );
    int comm_size = -1;
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size );

    /*
      Declare the AoSoA parameters.
    */
    using DataTypes = Cabana::MemberTypes<int, int>;
    const int VectorLength = 8;
    using MemorySpace = Kokkos::HostSpace;


    int num_tuple = 100;
    Cabana::AoSoA<DataTypes, MemorySpace, VectorLength> aosoa( "A", num_tuple );

    auto slice_ranks = Cabana::slice<0>( aosoa );
    auto slice_ids = Cabana::slice<1>( aosoa );
    for ( int i = 0; i < num_tuple; ++i )
    {
        slice_ranks( i ) = comm_rank;
        slice_ids( i ) = i;
    }

    if ( comm_rank == 0 )
    {
        std::cout << "BEFORE migration" << std::endl
                  << "(Rank " << comm_rank << ") ";
        for ( std::size_t i = 0; i < slice_ranks.size(); ++i )
            std::cout << slice_ranks( i ) << " ";
        std::cout << std::endl
                  << "(" << slice_ranks.size() << " ranks before migrate)"
                  << std::endl
                  << "(Rank " << comm_rank << ") ";
        for ( std::size_t i = 0; i < slice_ids.size(); ++i )
            std::cout << slice_ids( i ) << " ";
        std::cout << std::endl
                  << "(" << slice_ids.size() << " IDs before migrate)"
                  << std::endl
                  << std::endl;
    }


    Kokkos::View<int*, MemorySpace> export_ranks( "export_ranks", num_tuple );

    int previous_rank = ( comm_rank == 0 ) ? comm_size - 1 : comm_rank - 1;
    int next_rank = ( comm_rank == comm_size - 1 ) ? 0 : comm_rank + 1;
    for ( int i = 0; i < 10; ++i )
        export_ranks( i ) = next_rank;

    // Next 10 elements will be discarded. Use an export rank of -1 to
    // indicate this.
    for ( int i = 10; i < 20; ++i )
        export_ranks( i ) = -1;

    // The last 80 elements stay on this process.
    for ( int i = 20; i < num_tuple; ++i )
        export_ranks( i ) = comm_rank;

    std::vector<int> neighbors = { previous_rank, comm_rank, next_rank };
    std::sort( neighbors.begin(), neighbors.end() );
    auto unique_end = std::unique( neighbors.begin(), neighbors.end() );
    neighbors.resize( std::distance( neighbors.begin(), unique_end ) );
    Cabana::Distributor<MemorySpace> distributor( MPI_COMM_WORLD, export_ranks,
                                                  neighbors );



    Cabana::AoSoA<DataTypes, MemorySpace, VectorLength> destination(
        "destination", distributor.totalNumImport() );

    Cabana::migrate( distributor, aosoa, destination );

    auto slice_ranks_dst = Cabana::slice<0>( destination );
    auto slice_ids_dst = Cabana::slice<1>( destination );
    Cabana::migrate( distributor, slice_ranks, slice_ranks_dst );
    Cabana::migrate( distributor, slice_ids, slice_ids_dst );

    Cabana::migrate( distributor, aosoa );


    slice_ranks = Cabana::slice<0>( aosoa );
    slice_ids = Cabana::slice<1>( aosoa );

    if ( comm_rank == 0 )
    {
        std::cout << "AFTER migration" << std::endl
                  << "(Rank " << comm_rank << ") ";
        for ( std::size_t i = 0; i < slice_ranks.size(); ++i )
            std::cout << slice_ranks( i ) << " ";
        std::cout << std::endl
                  << "(" << slice_ranks.size() << " ranks after migrate)"
                  << std::endl
                  << "(Rank " << comm_rank << ") ";
        for ( std::size_t i = 0; i < slice_ids.size(); ++i )
            std::cout << slice_ids( i ) << " ";
        std::cout << std::endl
                  << "(" << slice_ids.size() << " IDs after migrate)"
                  << std::endl;
    }
}






void parse_config_file() {
  std::ifstream input_file("input.json");
  json j;
  input_file >> j;

  // Accessing the data
  for (const auto& param : j["parameters"]) {
//    printf("PARAM: %s\n", param["name"].get<std::string>().c_str());
//    printf("BIN_COUNT: %d\n", param["bin_count"].get<int>());
//    printf("MIN: %f\n", param["min"].get<double>());
//    printf("MAX: %f\n", param["max"].get<double>());
//    printf("MEAN: %f\n", param["mean"].get<double>());
//    printf("STDEV: %f\n", param["stdev"].get<double>());
//
//    // Iterate through bins
//    for (const auto& bin : param["bins"]) {
//      printf("BIN_MIN: %f, ", bin["bin_min"].get<double>());
//      printf("BIN_MAX: %f, ", bin["bin_max"].get<double>());
//      printf("BIN_PROP: %f, ", bin["bin_prop"].get<double>());
//      printf("BIN_MEAN: %f, ", bin["bin_mean"].get<double>());
//      printf("BIN_STDEV: %f\n", bin["bin_stdev"].get<double>());
//    }
      printf("test");
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
        bool config_file_used = false;//todo



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
      char unit_symbol = it->second.first;
      int unit_div = it->second.second;

      std::cout << "Unit Symbol: " << unit_symbol << std::endl;
      std::cout << "Unit Division: " << unit_div << std::endl;
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

  parse_config_file() ;
  MPI_Init( &argc, &argv );
  {
    Kokkos::ScopeGuard scope_guard( argc, argv );

    migrationExample();
  }
  MPI_Finalize();

//    parseArgs(argc, argv);




	printf("Hi\n");
	return 0;
}
