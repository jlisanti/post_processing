#include <iostream>
#include <vector>
#include <string>
#include <sstream>
bool is_number(const std::string& s)
{
    try
    {
        std::stod(s);
    }
    catch(...)
    {
        return false;
    }
    return true;
}

template <class value>
value table_value(int i, std::vector<std::string> table)
{
	if(is_number(table[i]))
		return std::stod(table[i]);
	else
	    return table[i];
}
int main()
{
	std::vector<std::string> table;
	table.push_back("stuff");
	table.push_back("4.0");
	table.push_back("3.0");
	table.push_back("2.0");
	table.push_back("1.0");

	std::string check = "3123";

	for (int i = 0; i < table.size(); i++)
		std::cout << i << '\t' << table[i] 
			      << '\t' << is_number(table[i]) 
				  << '\t' << table_value(i, table) << std::endl;
	
	return 0;
}
