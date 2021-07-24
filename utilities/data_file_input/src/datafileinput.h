#ifndef CODES_POST_PROCESSING_UTILITIES_DATAFILEINPUT_H_
#define CODES_POST_PROCESSING_UTILITIES_DATAFILEINPUT_H_

#include <vector>
#include <string>
#include <sstream>

class DataFileInput
{
protected:
    std::string m_input_file;
    int m_number_rows_skip;
	int m_number_rows;
	int m_number_cols;
	std::vector< std::vector<double> > m_table;
	//std::vector< std::vector<std::string> > m_table;
public:
    DataFileInput(std::string input_file, int number_rows_skip);
    bool is_numeric (const std::string& str); 
    void read_input_file(std::string input_file);
	/*
	template <class value>
		value table_value(int i, int j)
		{
			std::string tmp = m_table[i][j];
			if(is_numeric(tmp))
				return std::stod(tmp);
			else
				return tmp;
		}
		*/
    double table_value(int i, int j) const { return m_table[i][j]; }
	//double return_table_value(double value) const { return value; }
	//std::string return_table_value(std::string value) const { return value; }
	std::vector<double> table_column(int j) const { 
		std::vector<double> tmp;
		for (int i = 0; i < m_number_rows; i++)
			tmp.push_back(m_table[i][j]);
		return tmp;
	}
	//std::vector<std::string> table_column(int j) const
	//{
		//std::vector<double> tmp;
//		std::vector<std::string> tmp;
//		for (int i = 0; i < m_number_rows; i++)
//			tmp.push_back(m_table[i][j]);
//		return tmp; 
//	}
	std::vector<double> table_row(int i) const {
		std::vector<double> tmp;
		for (int j = 0; j < m_number_cols; j++)
			tmp.push_back(m_table[i][j]);
		return tmp;
	}
//	std::vector<std::string> table_row(int i) const
//	{ 
		//std::vector<double> tmp;
//		std::vector<std::string> tmp;
//		for (int j = 0; j < m_number_cols; j++)
//			tmp.push_back(m_table[i][j]);
//		return tmp;
//	}
	int file_length() const {return m_number_rows; }
	int length()      const {return m_number_rows; }
	int size()        const {return m_number_rows; }
	int number_cols() const {return m_number_cols; }
	int cols()        const {return m_number_cols; }
	int number_columns() const { return m_number_cols; }
    void delete_row(int i) 
    { 
		m_table[i].erase(m_table[i].begin(),m_table[i].begin() + m_number_cols); 
        m_number_rows--;
    }
	void clear() {
		for (int i = 0; i < m_number_rows; i++)
			m_table[i].clear();
	}
	void set_skip(int skip_rows) {m_number_rows_skip = skip_rows; };

	~DataFileInput()
	{
		m_table.clear();
	}
};

#endif
