#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include "misc.hpp"

bool IsFileExist(std::string fname)
{
  std::ifstream ifile(fname.c_str());
  return ifile.is_open();
}

std::string DecommentFile(std::string fname)
{
  std::stringstream msg;
  if (!IsFileExist(fname)) {
    msg << "### FATAL ERROR in DecommentFile. File " << fname << " doesn't exist." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  std::ifstream file(fname.c_str(), std::ios::in);
  std::string ss;
  char c;
  while (file) {
    file.get(c);
    if (c == '#') {
      while (c != '\n' && file) file.get(c);
      continue;
    }
    ss += c;
  }
  return ss;
}

void SplitString(std::string str, std::vector<std::string>& result)
{
  std::istringstream ss(str);
  result.clear();

  while (ss) {
    std::string sub;
    ss >> sub;
    result.push_back(sub);
  }

  // there is an empty space
  result.pop_back();
}

int GetNumCols(std::string fname, char c)
{
  std::ifstream inp(fname.c_str(), std::ios::in);
  std::string line;
  std::getline(inp, line);
  if (line.empty()) return 0;
  int cols = line[0] == c ? 0 : 1;

  for (int i = 1; i < line.length(); ++i)
    if (line[i-1] == c && line[i] != c) cols++;
  return cols;
}

int GetNumRows(std::string fname)
{
  std::ifstream inp(fname.c_str(), std::ios::in);
  std::string line;
  int rows = 0;

  while (std::getline(inp, line)) ++rows;
  return rows;
}

NamedArray ReadNamedArray(std::string fname)
{
  NamedArray amap;
  std::ifstream input(fname.c_str(), std::ios::in);
  std::stringstream ss, msg;
  std::string line, sbuffer; 
  std::vector<std::string> field;

  if (!input.is_open()) {
    msg << "### FATAL ERROR in ReadNamedArray. File " << fname << " doesn't exist." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  getline(input, line);
  ss.str(line);
  while (!ss.eof()) {
      ss >> sbuffer;
      field.push_back(sbuffer);
  }
  ss.clear();

  while (getline(input, line)) {
      if (line.empty()) continue;
      ss.str(line);
      for (std::vector<std::string>::iterator f = field.begin(); f != field.end(); ++f) {
          double value;
          ss >> value;
          amap[*f].push_back(value);
      }
      ss.clear();
  }

  return amap;
}
