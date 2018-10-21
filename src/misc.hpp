#ifndef MISC_HPP
#define MISC_HPP
#include <string>
#include <vector>
#include <map>

bool IsFileExist(std::string fname);

std::string DecommentFile(std::string fname);

void SplitString(std::string str, std::vector<std::string>& result);

int GetNumCols(std::string fname, char c = ' ');

int GetNumRows(std::string fname);

typedef std::map<std::string, std::vector<double> > NamedArray;
NamedArray ReadNamedArray(std::string fname);

#endif
