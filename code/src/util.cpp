#include "util.h"

namespace util
{
    std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            string temp_elem;
            int i = 0;
            while(isprint(item[i])){
                temp_elem.push_back(item[i]);
                i++;
            }
            if (temp_elem.size() > 0){
                elems.push_back(temp_elem);
            }
        }
        return elems;
    }

    std::vector<std::string> split(const std::string &s, char delim) {
        std::vector<std::string> elems;
        split(s, delim, elems);
        return elems;
    }

    std::string toLower(const std::string &s) {
      std::string new_string = s;
      std::transform(new_string.begin(), new_string.end(), new_string.begin(), (int (*)(int))std::tolower);
      return new_string;
    }

}
