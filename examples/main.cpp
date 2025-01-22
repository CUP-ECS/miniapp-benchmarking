#include <iostream>
#include <memory>
#include <type_traits>

#include <Cabana_Core.hpp>
#include <Kokkos_Core.hpp>

#include <string>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>










int main(int argc, char** argv)
{
    try {
        // Create the command line parser
        TCLAP::CmdLine cmd("\nNOTE: Setting parameters for the benchmark such as (neighbors, owned, remote, blocksize, and stride)"
                           "\tsets parameters to those values for the reference benchmark."
                           "\tThose parameters are then randomized for the irregular samples"
                           "\twhere the user-set parameters become averages for the random generation."
                           "\tUse the `--disable-irregularity` flag to only run the reference benchmark.", ' ', "1.0");

        // Arguments
        TCLAP::ValueArg<std::string> filepathArg("f", "filepath", "Path to the BENCHMARK_CONFIG file", false, "", "string");
        TCLAP::ValueArg<int> typeSizeArg("t", "typesize", "Size of the variable being sent (in bytes)", false, 0, "int");
        TCLAP::ValueArg<int> samplesArg("I", "samples", "Number of random samples to generate", false, 0, "int");
        TCLAP::ValueArg<int> iterationsArg("i", "iterations", "Number of updates each sample performs", false, 0, "int");
        TCLAP::ValueArg<int> neighborsArg("n", "neighbors", "Average number of neighbors each process communicates with", false, 0, "int");
        TCLAP::ValueArg<int> neighborsStdvArg("N", "neighbors_stdv", "Standard deviation of the number of neighbors each process communicates with", false, 0, "int");
        TCLAP::ValueArg<int> ownedAvgArg("o", "owned_avg", "Average byte count for data owned per node", false, 0, "int");
        TCLAP::ValueArg<int> ownedStdvArg("O", "owned_stdv", "Standard deviation byte count for data owned per node", false, 0, "int");
        TCLAP::ValueArg<int> remoteAvgArg("r", "remote_avg", "Average amount of data each process receives", false, 0, "int");
        TCLAP::ValueArg<int> remoteStdvArg("R", "remote_stdv", "Standard deviation of the amount of data each process receives", false, 0, "int");
        TCLAP::ValueArg<int> blockSizeAvgArg("b", "blocksize_avg", "Average size of transmitted blocks", false, 0, "int");
        TCLAP::ValueArg<int> blockSizeStdvArg("B", "blocksize_stdv", "Standard deviation of transmitted block sizes", false, 0, "int");
        TCLAP::ValueArg<int> strideArg("s", "stride", "Average size of stride", false, 0, "int");
        TCLAP::ValueArg<int> strideStdvArg("T", "stride_stdv", "Standard deviation of stride", false, 0, "int");
        TCLAP::ValueArg<int> seedArg("S", "seed", "Positive integer to be used as seed for random number generation", false, 0, "int");
        TCLAP::ValueArg<std::string> memSpaceArg("m", "memspace", "Choose from: host, cuda, openmp, opencl", false, "host", "string");
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
        cmd.add(memSpaceArg);
        cmd.add(distributionArg);
        cmd.add(unitsArg);
        cmd.add(reportParamsArg);
        cmd.add(disableirregularityArg);
//        cmd.setExtraUsage(
//            "\nNOTE: Setting parameters for the benchmark such as (neighbors, owned, remote, blocksize, and stride)\n"
//            "      sets parameters to those values for the reference benchmark.\n"
//            "      Those parameters are then randomized for the irregular samples\n"
//            "      where the user-set parameters become averages for the random generation.\n"
//            "      Use the `--disable-irregularity` flag to only run the reference benchmark.\n"
//        );

        // Parse the command line arguments
        cmd.parse(argc, argv);

        // Example of how you would retrieve and use the arguments
        std::string filepath = filepathArg.getValue();
        int typesize = typeSizeArg.getValue();
        int samples = samplesArg.getValue();

        printf("test1 %i\n",typesize);
        printf("test2 %i\n",samples);

    }
    catch (TCLAP::ArgException &e) {
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
        return -1;
    }

	printf("Hi\n");
	return 0;
}
