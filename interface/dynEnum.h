#include <unordered_map>
#include <string>

class DynEnum 
{
 private:
    int count = 1;
    std::unordered_map<std::string, int> map;
 public:
    int addEnum(std::string name) {
	map[name] = count;
	return count++;
    }
    int getEnum(std::string name) {
	return map.at(name);
    }
    bool foundEnum(std::string name) {
	return map.find(name) != map.end();
    }
};
