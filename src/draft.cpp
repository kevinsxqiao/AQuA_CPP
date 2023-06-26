#include <iostream>
#include <vector>
#include <algorithm>
//#include <boost/math/distributions/normal.hpp>

int main() {
    std::vector<int> x = {4,2,3,4,5};
    int node3 = *std::min_element(x.begin(), x.end()) ;
    std::cout<< node3;

}
