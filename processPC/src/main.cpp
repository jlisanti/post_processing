#include <iostream>
#include <vector>
#include <string>

#include "datafileinput.h"
#include "menu.h"

#include "readinputfiles.h"
#include "processfile.h"
#include "printoutput.h"

int main(int argc, char*argv[])
{
	/* Create menu */
	Menu cMenu("->");

	/* Store command line arguments */
	std::vector<std::string> arguments;
	for (int i = 0; i < argc; i++)
		arguments.push_back(argv[i]);

	/* Check input */
	if (argc <= 1)
	{
	    cMenu.display_message("Error : no input file");	
		exit(1);
	}
	else if (argc == 1)
	{
		/* No settings file */
		std::vector<std::string> type_menu;
		type_menu.push_back("pressureGain");
		type_menu.push_back("pulseCombustor");
    
		cMenu.add_menu("main",
			       type_menu);
		cMenu.display_menu("main");

		std::vector<std::string> file_names;
		read_input_files(argv[1],file_names);

		for (int i = 0; i < file_names.size(); i++)
		{
			std::cout << std::endl;
			std::cout << "Processing: " << file_names[i] << std::endl;
			//process_file(file_names[i],output);
			//print_output(output);
		}
	}
	else if ((argc > 1) && (arguments[2] == "-set"))
	{
		/* Settings file found */
		options optionsMenu;
		read_settings_file(argv[3],optionsMenu);
		std::vector<std::string> file_names;
		read_input_files(argv[1],file_names);

		for (int i = 0; i < file_names.size(); i++)
		{
			std::cout << std::endl;
			std::cout << "Processing: " << file_names[i] << std::endl;
			output_data output;
			process_file(file_names[i],output,optionsMenu);
			print_output_active(output,optionsMenu);
		}
	}
	else
	{
		/* Ambigious settings file declaration */
		cMenu.display_message("Error : use -set < settings file > ");
	}



	return 0;
}
