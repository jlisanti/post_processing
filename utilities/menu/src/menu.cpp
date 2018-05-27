#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "menu.h"

Menu::Menu(std::string prompt)
	: m_prompt(prompt)
{
}

void Menu::display_message(std::string input)
{
	std::cout << std::endl;
	std::cout << " " << input << std::endl; 
	std::cout << std::endl;
}

void Menu::end_line()
{
	std::cout << std::endl;
}

int Menu::yes_no_question(std::string question)
{
	std::cout << question << " [Y/n]" << std::endl;
	prompt();
	std::string response;
	std::cin >> response;
	if(check_yes(response))
		return 1;
	else if(check_no(response))
		return 0;
	else
		return -1;
}

bool Menu::check_yes(std::string input)
{
	if ((input=="yes") ||
		(input=="Yes") ||
		(input=="yEs") ||
		(input=="yeS") ||
        (input=="YEs") ||
		(input=="YES") ||
		(input=="yES") ||
		(input=="y")   ||
		(input=="Y")   ||
		(input=="ye")  ||
		(input=="Ye")  ||
		(input=="YE")  ||
		(input=="yE")  ||
		(input=="es")  ||
		(input=="Es")  ||
		(input=="eS")  ||
		(input=="ES"))
		return true;
	else
		return false;
}

bool Menu::check_no(std::string input)
{
	if ((input=="no") ||
		(input=="No") ||
		(input=="nO") ||
		(input=="NO") ||
        (input=="n") ||
		(input=="N"))
		return true;
	else
		return false;
}

void Menu::add_menu(std::string menu_name,
		            std::vector<std::string> menu_choices)
{
	for (int i = 0; i < menu_choices.size(); i++)
		m_menus[menu_name][i] = menu_choices[i];
}

void Menu::delete_menu(std::string menu_name)
{
}

void Menu::display_menu(std::string menu)
{
	for (int i = 0; i < m_menus[menu].size(); i++)
		std::cout << "[" << i << "] " << m_menus[menu][i] << std::endl;
}
