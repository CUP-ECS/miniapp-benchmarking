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


		TCLAP::CmdLine cmd("Command description message", ' ', "0.9");
		TCLAP::ValueArg<std::string> filepathArg("f","filepath","tspecify the path to the BENCHMARK_CONFIG file",true,"XXX","string");
		TCLAP::ValueArg<int>         typesizeArg("t","typesize","the size of the variable being sent (in bytes)",true,0,"int");
		TCLAP::ValueArg<int>         samplesArg("I","samples","the number of random samples to generate",true,0,"int");
		TCLAP::ValueArg<int>         neighborsArg("n","neighbors","specify the average number of neighbors each process communicates with",true,0,"int");
		TCLAP::ValueArg<int>         iterationsArg("i","iterations","specify the number of updates each sample performs",true,0,"int");
		TCLAP::ValueArg<int>         neighbors_stdvArg("N","neighbors_stdv","specify the stdev number of neighbors each process communicates with",true,0,"int");
		TCLAP::ValueArg<int>         owned_avgArg("o","owned_avg","x",true,0,"int");
		TCLAP::ValueArg<int>         owned_stdvArg("O","owned_stdv","x",true,0,"int");
		TCLAP::ValueArg<int>         remote_avgArg("r","remote_avg","x",true,0,"int");
		TCLAP::ValueArg<int>         remote_stdArg("R","remote_std","x",true,0,"int");
		TCLAP::ValueArg<int>         blocksize_avgArg("b","blocksize_avg","x",true,0,"int");
		TCLAP::ValueArg<int>         blocksize_stdvArg("B","blocksize_stdv","x",true,0,"int");
		TCLAP::ValueArg<int>         strideArg("s","stride","x",true,0,"int");
		TCLAP::ValueArg<int>         seedArg("S","seed","x",true,0,"int");
		TCLAP::ValueArg<int>         stride_stdvArg("T","stride_stdv","x",true,0,"int");
		TCLAP::ValueArg<int>         memspaceArg("m","memspace","x",true,0,"int");
		TCLAP::ValueArg<int>         distributionArg("d","distribution","x",true,0,"int");
		TCLAP::ValueArg<int>         unitsArg("u","units","x",true,0,"int");



        TCLAP::SwitchArg reverseSwitch("t","typesize","the size of the variable being sent (in bytes)", cmd, false);




        cmd.add(filepathArg);
        cmd.add(typesizeArg);
        cmd.add(samplesArg);
        cmd.add(neighborsArg);
        cmd.add(iterationsArg);
        cmd.add(neighbors_stdvArg);
        cmd.add(owned_avgArg);
        cmd.add(owned_stdvArg);
        cmd.add(remote_avgArg);
        cmd.add(remote_stdArg);
        cmd.add(blocksize_avgArg);
        cmd.add(blocksize_stdvArg);
        cmd.add(strideArg);
        cmd.add(seedArg);
        cmd.add(stride_stdvArg);
        cmd.add(memspaceArg);
        cmd.add(distributionArg);
        cmd.add(unitsArg);
;




		// Parse the argv array.
		cmd.parse( argc, argv );

		// Get the value parsed by each arg.
		std::string name = filepathArg.getValue();
		std::cout << "My name is: " << name << std::endl;





	} catch (TCLAP::ArgException &e)  // catch exceptions
	{
          std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;

        }


	printf("Hi\n");
	return 0;
}
