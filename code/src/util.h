#include "main.h"

namespace util

{
  std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
  std::vector<std::string> split(const std::string &s, char delim);
  std::string toLower(const std::string &s);
}
