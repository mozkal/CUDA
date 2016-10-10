#include <algorithm>
#include <cassert>
#include <cstring>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "common.h"
#include "declarations.h"

_result_vector results;

void algorythm(const Fragment &frag1, const Fragment &frag2) {
    int max = needleman_wunsch(frag1.seq.c_str(), frag2.seq.c_str());
    if(max >= WEIGHT_LVL)
        results.emplace_back(frag1, frag2, max);
}

struct final_comparator {
    bool operator() (const _tuple &frag1, const _tuple &frag2)
    {
        long id11 = std::get<0>(frag1).id;
        long id12 = std::get<1>(frag1).id;
        long id21 = std::get<0>(frag2).id;
        long id22 = std::get<1>(frag2).id;
        return (id11 < id21 || (id11 == id21 && id12 < id22) );
    }
} final_comp;

//zwraca zakres sasiedztwa danej sekwencji zgodnie z heurystyką
unsigned evaluate_range(unsigned i, unsigned long len) {
    unsigned range;

    if(i < NEIGHBOURS)
        range = std::min((unsigned long) 2 * NEIGHBOURS - i, len - 1 - i);
    else if(i + 2 * NEIGHBOURS + 1 > len - 1)
        range = len - 1 - i;
    else
        range = NEIGHBOURS;

    return range;
}

void analyze(std::vector<Fragment> &seq) {
    unsigned i = 0, j, range;
    unsigned long len = seq.size();
    for(auto &s : seq) {
        range = evaluate_range(i, len);
        for(j = 1; j <= range; j++)
            algorythm(s, seq.at(i + j));

        ++i;
    }
}

void print(_result_vector &results) {
    std::sort(results.begin(), results.end(), final_comp);

    for(auto result : results)
        std::cout << std::get<0>(result).id << " " << std::get<1>(result).id
            << " " << std::get<0>(result).comp << " " << std::get<1>(result).comp
            << " " << std::get<2>(result)<< "\n";
}

int main() {
    std::vector<Fragment> seq = prepare_set();
    analyze(seq);
    std::sort(results.begin(), results.end(), final_comp);
    print(results);

    return 0;
}