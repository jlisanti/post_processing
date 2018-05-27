#include <iostream>

#include "datafileinput.h"

int main()
{
    DataFileInput cInput("test.dat", 3);
    std::cout << cInput.table_value(1,3) << std::endl;
	int row = 3;
	std::vector<double> test_rows = cInput.table_row(row);
	for (int i = 0; i < test_rows.size(); i++)
		std::cout << "Row " << row << " element " << i << " " << test_rows[i] << std::endl;

	int col = 3; 
	std::vector<double> test_cols = cInput.table_column(col);
	for (int i = 0; i < test_cols.size(); i++)
		std::cout << "Col " << col << " element " << i << " : " << test_cols[i] << std::endl;
    return 0;
}
