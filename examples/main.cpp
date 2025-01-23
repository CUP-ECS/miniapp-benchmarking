#include <iostream>
#include <memory>
#include <type_traits>

#include <Cabana_Core.hpp>
#include <Kokkos_Core.hpp>

#include <string>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>

#include <filesystem>


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
static int           unit_div = 1;
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

void parse_arguments(int argc, char **argv){
  try {
        // Create the command line parser
        TCLAP::CmdLine cmd("\nNOTE: Setting parameters for the benchmark such as (neighbors, owned, remote, blocksize, and stride)"
                           "sets parameters to those values for the reference benchmark."
                           "Those parameters are then randomized for the irregular samples"
                           "where the user-set parameters become averages for the random generation."
                           "Use the `--disable-irregularity` flag to only run the reference benchmark.", ' ', "1.0");

        // Arguments
        TCLAP::ValueArg<std::string> filepathArg("f", "filepath", "Path to the BENCHMARK_CONFIG file", true, "", "string");
        TCLAP::ValueArg<int> typeSizeArg("t", "typesize", "Size of the variable being sent (in bytes)", false, 8, "int");
        TCLAP::ValueArg<int> samplesArg("I", "samples", "Number of random samples to generate", false, 25, "int");
        TCLAP::ValueArg<int> iterationsArg("i", "iterations", "Number of updates each sample performs", false, 100, "int");
        TCLAP::ValueArg<int> neighborsArg("n", "neighbors", "Average number of neighbors each process communicates with", true, -1, "int");
        TCLAP::ValueArg<int> neighborsStdvArg("N", "neighbors_stdv", "Standard deviation of the number of neighbors each process communicates with", true, -1, "int");
        TCLAP::ValueArg<int> ownedAvgArg("o", "owned_avg", "Average byte count for data owned per node", true, -1, "int");
        TCLAP::ValueArg<int> ownedStdvArg("O", "owned_stdv", "Standard deviation byte count for data owned per node", true, -1, "int");
        TCLAP::ValueArg<int> remoteAvgArg("r", "remote_avg", "Average amount of data each process receives", true, -1, "int");
        TCLAP::ValueArg<int> remoteStdvArg("R", "remote_stdv", "Standard deviation of the amount of data each process receives", true, -1, "int");
        TCLAP::ValueArg<int> blockSizeAvgArg("b", "blocksize_avg", "Average size of transmitted blocks", true, -1, "int");
        TCLAP::ValueArg<int> blockSizeStdvArg("B", "blocksize_stdv", "Standard deviation of transmitted block sizes", true, -1, "int");
        TCLAP::ValueArg<int> strideArg("s", "stride", "Average size of stride", true, -1, "int");
        TCLAP::ValueArg<int> strideStdvArg("T", "stride_stdv", "Standard deviation of stride", true, -1, "int");
        TCLAP::ValueArg<int> seedArg("S", "seed", "Positive integer to be used as seed for random number generation", true, -1, "int");
//        TCLAP::ValueArg<std::string> memSpaceArg("m", "memspace", "Choose from: host, cuda, openmp, opencl", false, "host", "string");
        TCLAP::ValueArg<std::string> distributionArg("d", "distribution", "Choose from: gaussian (default), empirical", false, "gaussian", "string");
        TCLAP::ValueArg<std::string> unitsArg("u", "units", "Choose from: a,b,k,m,g (auto, bytes, kilobytes, etc.)", false, "auto", "string");

        // Optional flag argument for reporting
        TCLAP::SwitchArg reportParamsArg("", "report-params", "Enables parameter reporting for use with analysis scripts", false);
        TCLAP::SwitchArg disableirregularityArg("", "disable-irregularity", "Use the `--disable-irregularity` flag to only run the reference benchmark.", false);


        // Add all arguments to the command line parser
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
//        cmd.add(memSpaceArg);
        cmd.add(distributionArg);
        cmd.add(unitsArg);
        cmd.add(reportParamsArg);
        cmd.add(disableirregularityArg);
        cmd.parse(argc, argv);


        filepath = filepathArg.getValue();
        bool config_file_used = false;



        if (!std::filesystem::exists(filepath)) {
            printf("ERROR: the specified filepath doesn't exist, exiting...\n");
            exit(-1);
        }
        config_file_used = true;
        std::string dist = distributionArg.getValue();



         typesize = typeSizeArg.getValue();

         if (typesize < 1 || typesize > 8) {
             printf("ERROR: Invalid typesize\n");
             exit(-1);
         }


//numpes =00; \\todo
        nsamples=samplesArg.getValue();

        if (nsamples < 0){
            printf("ERROR: Invalid number of samples\n");
            exit(-1);
        }

        niterations=iterationsArg.getValue();
        if (niterations < 0){
          printf("ERROR: Invalid number of iterations\n");
          exit(-1);
        }




nneighbors=neighborsArg.getValue();
if (nneighbors < 0){
    printf("ERROR: Invalid number of neighbors\n");
    exit(-1);
}


nneighbors_stdv=neighborsStdvArg.getValue();

if (nneighbors_stdv < 0){
    printf("ERROR: Invalid nighbors std\n");
    exit(-1);
}




nowned=ownedAvgArg.getValue();
if (nneighbors_stdv < 0){
    printf("ERROR: Invalid number of owned\n");
    exit(-1);
}

nowned_stdv=ownedStdvArg.getValue();
if (nneighbors_stdv < 0){
    printf("ERROR: Invalid owned std\n");
    exit(-1);
}


nremote=remoteAvgArg.getValue();
if (nremote < 0){
    printf("ERROR: Invalid number of remote\n");
    exit(-1);
}
nremote_stdv=remoteStdvArg.getValue();
if (nremote_stdv < 0){
    printf("ERROR: Invalid remote std\n");
    exit(-1);
}

blocksz=blockSizeAvgArg.getValue();
if (blocksz < 0){
  printf("ERROR: Invalid block size\n");
  exit(-1);
}
blocksz_stdv=blockSizeStdvArg.getValue();
if (blocksz_stdv < 0){
  printf("ERROR: Invalid block size\n");
  exit(-1);
}
stride=strideArg.getValue();
if (stride < 0){
  printf("ERROR: Invalid stride\n");
}
stride_stdv=strideStdvArg.getValue();
if (stride_stdv < 0){
  printf("ERROR: Invalid stride std\n");
   exit(-1);

}


std::string unit = unitsArg.getValue();

if(unit=="auto" ||unit=="a" ) {
    unit_symbol=A;
    unit_div=0;
}else if(unit=="bytes" ||unit=="b") {
    unit_symbol=A;
    unit_div=0;
}else if(unit=="kilobytes" ||unit=="k") {
    unit_symbol=K;
    unit_div=1024 ;
}else if(unit=="megabytes" ||unit=="m") {
    unit_symbol=M;
    unit_div=1024*1024;
}else if(unit=="gigabytes" ||unit=="g") {
    unit_symbol=G;
    unit_div=1024*102*10244;
}else {
  printf("ERROR: Invalid formatting choice [b, k, m, g]\n");
  exit(-1);
}



std::string distribution=distributionArg.getValue();

if (distribution == "gaussian"  ||
    distribution == "g") {
  distribution_t = GAUSSIAN;

}else if (distribution == "empirical"
          || distribution == "e") {
  distribution_t=EMPIRICAL;
}else {
  printf("ERROR: Invalid distribution choice [empirical,gaussian]\n");
  exit(-1);

}



seed=seedArg.getValue();
if(seed==-1) {
  seed=time(NULL);
}
srand(seed);






irregularity=irregularityArg.getValue();

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
