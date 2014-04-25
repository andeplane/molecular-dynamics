#ifndef FILEMANAGER_H
#define FILEMANAGER_H
#include <string>
using std::string;

class System;

class FileManager
{
public:
    FileManager();
    void loadState(string stateFolder, System &system);
    void saveState(string stateFolder, System &system);
};

#endif // FILEMANAGER_H
