#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

char atomNames[][5] = {"NA", "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar"};

int main()
{
    int nodeIndex = 0;

    char *filename = new char[10000];
    sprintf(filename,"../example/movie_files/movie%04d.bin",nodeIndex);
    ifstream movieFile(filename,std::ios::in | std::ios::binary);
    sprintf(filename,"movie.xyz");
    ofstream xyzFile(filename);

    delete filename;

    vector<double> dataArray;
    int timesteps = 0;

    char *tmpString = new char[10000];

    while(!movieFile.eof()) {
        unsigned long numberOfAtoms = 0;
        movieFile.read(reinterpret_cast<char*>(&numberOfAtoms), sizeof(unsigned long));
        if(numberOfAtoms == 0) break;
        cout << "Number of atoms: " << numberOfAtoms << endl;

        dataArray.resize(5*numberOfAtoms);
        movieFile.read(reinterpret_cast<char*>(&dataArray.front()), 5*numberOfAtoms*sizeof(double));

        sprintf(tmpString, "%ld\nstuff \n", numberOfAtoms);
        xyzFile << tmpString;
        for(unsigned long i=0; i<numberOfAtoms; i++) {
            double x = dataArray[5*i + 0];
            double y = dataArray[5*i + 1];
            double z = dataArray[5*i + 2];
            int atomType = int(dataArray[5*i + 3]);
            int atomID = int(dataArray[5*i + 4]);
            xyzFile << atomNames[atomType] << " " << x << " " << y << " " << z << " " << atomID << endl;
        }

        timesteps++;
    }

    xyzFile.close();
    movieFile.close();
    cout << "Created xyz with " << timesteps << " timesteps." << endl;
    return 0;
}

