

#include <bits/stdc++.h>
#include <fstream>

using namespace std;



struct point
{

	float x;
	float y;
	float z;


};

int main()
{

	ifstream file;
	file.open("vert.txt");
	map <int , point > map;

	int info;
	point point1;
	int greaterinfo = 0;

	while (!file.eof())
	{
		file >> info;
		file >> point1.x >> point1.y >> point1.z;
		map[info] = point1;



		if (info > greaterinfo)
		{

			greaterinfo = info;
		}



	}


	ofstream ofile;


	ofile.open("overt.ply", std::ofstream::out | std::ofstream::app);



	for (int i = 0; i <= greaterinfo; ++i)
	{

		ofile << map[i].x << " " << map[i].y << " " << map[i].z << endl;

	}



	ifstream infile("input_verti_info.txt");
	//ofs.open ("overt.txt", std::ofstream::out | std::ofstream::app);
	string content = "";
	int i;

	for (i = 0 ; infile.eof() != true ; i++) // get content of infile
		content += infile.get();

	i--;
	content.erase(content.end() - 1);   // erase last character


	infile.close();

	ofile << content;                 // output


	int countofnewline = 0;
	for (size_t i = 0; i < content.length(); i++)
		if (content[i] == '\n')
			countofnewline++;
	// count of newlines is no of facets;


	cout << greaterinfo + 1 << endl;


	cout << countofnewline << endl;

	ofile.close();



}




