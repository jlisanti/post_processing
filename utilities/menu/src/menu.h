#ifndef CODES_POST_PROCESSING_UTILITIES_MENU_H_
#define CODES_POST_PROCESSING_UTILITIES_MENU_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>

class Menu
{
protected:
    std::string m_prompt;
    std::map< std::string, std::map<int, std::string> > m_menus;
public:
	Menu(std::string prompt);

	void prompt() const {std::cout << m_prompt;};
    void display_message(std::string input);
	void end_line();

	int  yes_no_question(std::string question);

	bool check_yes(std::string input);
	bool check_no(std::string input);
	
	/* functions for creating and displaying menu */
	void add_menu(std::string menu_name,
			      std::vector<std::string>);
	void delete_menu(std::string menu_name);
	void display_menu(std::string menu);
};


#endif
