#include <iostream>

#include "menu.h"

int main()
{
	Menu cMenu("-> ");
	//cMenu.prompt();
	int check = cMenu.yes_no_question("Does this work");
    if(check==1)
		cMenu.display_message("answer yes");
	else if (check==0)
		cMenu.display_message("answer no");
	else
		cMenu.display_message("unrecognized response");
	return 0;
}
