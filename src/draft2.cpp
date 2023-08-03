#include <vector>
#include <unordered_map>
#include <algorithm>

std::vector<int> findMostFrequentElements(const std::vector<int>& pix0, int t_scl) {
    std::unordered_map<int, int> counts;
    for (const int& num : pix0) {
        ++counts[num];
    }

    std::vector<int> fgPix;
    for (const auto& pair : counts) {
        if (pair.second > t_scl / 2) {
            fgPix.push_back(pair.first);
        }
    }

    return fgPix;
}
