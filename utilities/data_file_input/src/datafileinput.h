#ifndef CODES_POST_PROCESSING_UTILITIES_DATAFILEINPUT_H_
#define CODES_POST_PROCESSING_UTILITIES_DATAFILEINPUT_H_

#include <vector>
#include <string>

class DataFileInput
{
protected:
    std::string m_input_file;
    int m_number_rows_skip;
	int m_number_rows;
	int m_number_cols;
	std::vector< std::vector<double> > m_table;
public:
    DataFileInput(std::string input_file, int number_rows_skip);
    void read_input_file(std::string input_file);
    double table_value(int i, int j) const { return m_table[i][j]; }
	std::vector<double> table_column(int j) const 
	{
		std::vector<double> tmp;
		for (int i = 0; i < m_number_rows; i++)
			tmp.push_back(m_table[i][j]);
		return tmp; 
	}
	std::vector<double> table_row(int i) const 
	{ 
		std::vector<double> tmp;
		for (int j = 0; j < m_number_cols; j++)
			tmp.push_back(m_table[i][j]);
		return tmp;
	}
	int file_length() const {return m_number_rows; }
	~DataFileInput()
	{
		m_table.clear();
	}
};

#endif
