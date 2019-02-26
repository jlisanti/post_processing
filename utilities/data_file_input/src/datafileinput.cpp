#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

#include "datafileinput.h"


DataFileInput::DataFileInput(std::string input_file, int number_rows_skip)
    : m_input_file(input_file), m_number_rows_skip(number_rows_skip)
{
    read_input_file(input_file);
}


void DataFileInput::read_input_file(std::string input_file)
{
    std::ifstream input_file_stream(input_file);
    std::string dummy_line;
    for (int i = 0; i < m_number_rows_skip; i++)
        getline(input_file_stream, dummy_line);

	int number_rows = 0;
	int number_cols = 0;

    while (true)
    {
        std::string line;
        double buffer;
		//std::string buffer;
        getline(input_file_stream, line);
        std::istringstream ss(line, std::ios_base::out|std::ios_base::in|
                std::ios_base::binary);

        if (!input_file_stream)
            break;
        if (line[0] == '#' || line.empty())
            continue;

		std::vector<double> row;
		//std::vector<std::string> row;
		while (ss >> buffer)
		{
			if(number_rows == 0)
			{
				row.push_back(buffer);
				number_cols++;
			}
			else
				    row.push_back(buffer);
		}

		m_table.push_back(row);

		number_rows++;
    }

	m_number_rows = number_rows;
	m_number_cols = number_cols;
}

bool DataFileInput::is_numeric(const std::string& s)
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

/*
bool DataFileInput::is_numeric (const std::string& str) {
    std::istringstream ss(str);
    double dbl;
    ss >> dbl;      // try to read the number
    ss >> std::ws;  // eat whitespace after number

    if (!ss.fail() && ss.eof()) {
        return true;  // is-a-number
    } else {
        return false; // not-a-number
    }
}
*/
