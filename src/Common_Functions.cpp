/*
 * Common_Functions.cpp
 *
 *  Created on: Aug 21, 2020
 *      Author: rodney
 */

#include "Common_Functions.hpp"

Common_Functions::Common_Functions() {

}

Common_Functions::~Common_Functions() {

}

void Common_Functions::import_csv(string filename, vector<double> &col1, vector<double> &col2)
{
	ifstream infile;
	infile.open(filename);
	if (!infile.good())
	{
		cout << "\"" << filename << "\" not found!\n";
		exit(1);
	}
	string line, word;
	vector<string> row;
	while (getline(infile, line))
	{
		row.clear();
		if (line[0] == '#' || line.empty() || std::all_of(line.begin(), line.end(), ::isspace))
		{
			continue;
		}
		else
		{
			stringstream str(line);
			while(getline(str, word, ','))
			{
				row.push_back(word);
			}
			col1.push_back(stod(row[0]));
			col2.push_back(stod(row[1]));
		}
	}
	infile.close();
}

void Common_Functions::import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3)
{
	ifstream infile;
	infile.open(filename);
	if (!infile.good())
	{
		cout << "\"" << filename << "\" not found!\n";
		exit(1);
	}
	string line, word;
	vector<string> row;
	while (getline(infile, line))
	{
		row.clear();
		if (line[0] == '#' || line.empty() || std::all_of(line.begin(), line.end(), ::isspace))
		{
			continue;
		}
		else
		{
			stringstream str(line);
			while(getline(str, word, ','))
			{
				row.push_back(word);
			}
			col1.push_back(stod(row[0]));
			col2.push_back(stod(row[1]));
			col3.push_back(stod(row[2]));
		}
	}
	infile.close();
}

void Common_Functions::import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4)
{
	ifstream infile;
	infile.open(filename);
	if (!infile.good())
	{
		cout << "\"" << filename << "\" not found!\n";
		exit(1);
	}
	string line, word;
	vector<string> row;
	while (getline(infile, line))
	{
		row.clear();
		if (line[0] == '#' || line.empty() || std::all_of(line.begin(), line.end(), ::isspace))
		{
			continue;
		}
		else
		{
			stringstream str(line);
			while(getline(str, word, ','))
			{
				row.push_back(word);
			}
			col1.push_back(stod(row[0]));
			col2.push_back(stod(row[1]));
			col3.push_back(stod(row[2]));
			col4.push_back(stod(row[3]));
		}
	}
	infile.close();
}

void Common_Functions::import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4, vector<double> &col5)
{
	ifstream infile;
	infile.open(filename);
	if (!infile.good())
	{
		cout << "\"" << filename << "\" not found!\n";
		exit(1);
	}
	string line, word;
	vector<string> row;
	while (getline(infile, line))
	{
		row.clear();
		if (line[0] == '#' || line.empty() || std::all_of(line.begin(), line.end(), ::isspace))
		{
			continue;
		}
		else
		{
			stringstream str(line);
			while(getline(str, word, ','))
			{
				row.push_back(word);
			}
			col1.push_back(stod(row[0]));
			col2.push_back(stod(row[1]));
			col3.push_back(stod(row[2]));
			col4.push_back(stod(row[3]));
			col5.push_back(stod(row[4]));
		}
	}
	infile.close();
}

void Common_Functions::import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4, vector<double> &col5, vector<double> &col6)
{
	ifstream infile;
	infile.open(filename);
	if (!infile.good())
	{
		cout << "\"" << filename << "\" not found!\n";
		exit(1);
	}
	string line, word;
	vector<string> row;
	while (getline(infile, line))
	{
		row.clear();
		if (line[0] == '#' || line.empty() || std::all_of(line.begin(), line.end(), ::isspace))
		{
			continue;
		}
		else
		{
			stringstream str(line);
			while(getline(str, word, ','))
			{
				row.push_back(word);
			}
			col1.push_back(stod(row[0]));
			col2.push_back(stod(row[1]));
			col3.push_back(stod(row[2]));
			col4.push_back(stod(row[3]));
			col5.push_back(stod(row[4]));
			col6.push_back(stod(row[5]));
		}
	}
	infile.close();
}

// returns interpolated value at x from parallel arrays (x_data, y_data)
// assumes that x_data has at least two elements, is sorted and is strictly monotonic increasing
double Common_Functions::interpolate(vector<double> &x_data, vector<double> &y_data, double x)
{
	int size = x_data.size();

	int i = 0;                    // find left end of interval for interpolation
	if (x >= x_data[size - 2])    // special case: beyond right end
	{
		i = size - 2;
	}
	else
	{
		while (x > x_data[i+1]) i++;
	}
	double x_left = x_data[i], y_left = y_data[i], x_right = x_data[i+1], y_right = y_data[i+1];  // points on either side (unless beyond ends)

	double dydx = (y_right - y_left) / (x_right - x_left);   // gradient

	return y_left + dydx * (x - x_left);   // linear interpolation
}

// returns uniformly distributed random number between 0 and 1
double Common_Functions::get_rand()
{
	return (double)rand() / (double)RAND_MAX;
}