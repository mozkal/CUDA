#include <cassert>
#include <cstring>
#include <set>
#include <sstream>
#include "declarations.h"
#include "common.h"

//Komparator zgodny z heurystyką
struct h_comparator {
    bool operator() (const Fragment& frag1, const Fragment& frag2)
    {
        return (frag1.gram_3 < frag2.gram_3);
    }
} comp;

void replace(std::string &s) {
    for(char &c : s)
        switch(c) {
            case 'A' : c = 'T'; break;
            case 'a' : c = 't'; break;
            case 'C' : c = 'G'; break;
            case 'c' : c = 'g'; break;
            case 'T' : c = 'A'; break;
            case 't' : c = 'a'; break;
            case 'G' : c = 'C'; break;
            case 'g' : c = 'c'; break;
        }
}

void add_complement(std::vector<Fragment> &frags) {
    std::string s_comp, str;
    std::vector<Fragment> comp_frags;

    for(auto &f : frags) {
        str = f.seq;
        s_comp = std::string((str).rbegin(), (str).rend());
        replace(s_comp);
        comp_frags.emplace_back(f.id, s_comp, true);
    }

    frags.insert(frags.end(), comp_frags.begin(), comp_frags.end());
}

std::string sort(const std::multiset<std::string> &grams) {
    std::string sorted_grams;
    std::set<std::string> set(grams.begin(), grams.end());
    std::vector<std::pair<int, std::string>> vector;
    
    vector.reserve(set.size());
    for(auto s : set)
        vector.emplace_back(grams.count(s), s);

    std::sort(vector.begin(), vector.end(),
              [&](const std::pair<int, std::string> &i1, const std::pair<int, std::string> &i2) {
                  return (i1.first > i2.first || (i1.first == i2.first && i1.second < i2.second));
              }
    );

    for(auto &v : vector)
        for(int i = 0; i < v.first; i++)
            sorted_grams += v.second;

    return sorted_grams;
}

std::string cut(const std::string s) {
    std::multiset<std::string> grams;

    assert(s.length() >= 3);
    for(unsigned i = 0; i <= s.length() - 3; i++) {
        std::string gram = s.substr(i, 3);
        grams.insert(gram);
    }

    return sort(grams);
}

//Zwraca najlepsze przyrównanie
int needleman_wunsch(const char* seq1, const char* seq2) {
    int i, j, current, prev, max = 0;
    int* prev_row = (int*) malloc(SINGLE_SEQ_LENGTH * sizeof(int));
    
    memset(prev_row, 0, SINGLE_SEQ_LENGTH * sizeof(int));

    for(i = 0; i < SINGLE_SEQ_LENGTH; i++) {
        for(j = 0; j < SINGLE_SEQ_LENGTH; j++) {
            current = (seq1[i] == seq2[j]) ? 1 : -2;

            if(j > 0) {
                current = std::max({current + prev_row[j - 1], prev_row[j] - GP, prev - GP});

                prev_row[j - 1] = prev;
                if(j == SINGLE_SEQ_LENGTH - 1)
                    prev_row[j] = current;
            } else
                current = std::max({current, prev_row[j] - GP, -GP});

            prev = current;
        }
        if(current >= max)
            max = current;
    }

    for(i = 0; i < SINGLE_SEQ_LENGTH; i++)
        if(prev_row[i] >= max)
            max = prev_row[i];

    free(prev_row);
    
    return max;
}

//Zwraca id z naglowka FASTA
long extract_id(std::string header) {
    std::stringstream ss(header);
    std::string str;

    std::getline(ss, str, '|');
    std::getline(ss, str, '|');

    return std::stol(str);
}

//Wczytuje sekwencje w formacie FASTA ze standardowego wejscia
std::vector<Fragment> fasta_input() {
    std::vector<Fragment> seq;
    std::string line, header, sequence;

    while(std::getline(std::cin, line)) {
        if(line[0] == '>') {
            if(!header.empty()) {

                seq.push_back(Fragment(extract_id(header), sequence));
                header.clear();
            }
            if(!line.empty()){
                header = line.substr(1);
            }
            sequence.clear();
        } else
            sequence += line;
    }
    if(!header.empty()) {

        seq.push_back(Fragment(extract_id(header), sequence));
    }

    return seq;
}

std::vector<Fragment> prepare_set() {
    std::vector<Fragment> seq = fasta_input();
   
    add_complement(seq);
    for(auto &s : seq) {
        s.gram_3 = cut(s.seq);
    }
    std::sort(seq.begin(), seq.end(), comp);

    return seq;
}

char* prepare_seq(const std::vector<Fragment> &seq) {
    std::string sequences;
    
    for(unsigned long i = 0; i < seq.size(); i++) {
        sequences += seq.at(i).seq;
    }
    char* return_seq = new char[sequences.length() + 1];
    strcpy(return_seq, sequences.c_str());
    
    return return_seq;
}

int* calculate_cpu(const std::vector<Fragment> &seq) {
    int* cpu_out = (int*) malloc(sizeof(int) * (NEIGHBOURS + 1) * NEIGHBOURS / 2);
    
    memset(cpu_out, -1, 129*64);
    
    int count = 0;
    for(int i = 0; i < NEIGHBOURS; i++)
        for(int j = i + NEIGHBOURS + 1; j <= 2 * NEIGHBOURS; j++) {
            cpu_out[count] = needleman_wunsch(seq.at(i).seq.c_str(), seq.at(j).seq.c_str());
            count++;
        }

    return cpu_out;
}

void final_print(_result_vector &results) {
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
    
    std::sort(results.begin(), results.end(), final_comp);

    for(auto result : results)
        std::cout << std::get<0>(result).id << " " << std::get<1>(result).id
            << " " << std::get<0>(result).comp << " " << std::get<1>(result).comp
            << " " << std::get<2>(result)<< "\n";
}

void print_results(int* cpu_out, int* gpu_out, std::vector<Fragment> const &seq) {
    _result_vector results;

    //Rekonstrukcja z gpu
    unsigned long idx, idx2;
    unsigned long size = seq.size();
    for(unsigned long i = 0; i < size; i++)
        for(int j = 0; j < 128; j++) {
            idx = i * NEIGHBOURS + j;
            idx2 = i + j + 1;
            if(idx2 >= size)
                idx2 = size - j - NEIGHBOURS - 1 - 1;

            if(gpu_out[idx] >= WEIGHT_LVL)
                results.emplace_back(seq.at(i), seq.at(idx2), gpu_out[idx]);
        }

    //Rekonstrukcja z cpu
    int count = 0;
    for(int i = 0; i < NEIGHBOURS; i++)
        for(int j = i + NEIGHBOURS + 1; j <= 2 * NEIGHBOURS; j++) {
            if(cpu_out[count] >= WEIGHT_LVL)
                results.emplace_back(seq.at(i), seq.at(j), cpu_out[count]);

            count++;
        }

    final_print(results);
}

