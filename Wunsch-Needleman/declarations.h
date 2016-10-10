#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <vector>

#define GP 2
#define NEIGHBOURS 128
#define WEIGHT_LVL 10
#define SINGLE_SEQ_LENGTH 60

struct Fragment {
    long id;
    std::string seq;
    std::string gram_3;
    bool comp;

    Fragment() {};
    Fragment(long id, std::string seq) : id(id), seq(seq), comp(0) {};
    Fragment(long id, std::string seq, bool comp) : id(id), seq(seq), comp(comp) {};
};

std::vector<Fragment> prepare_set();
char* prepare_seq(const std::vector<Fragment> &seq);
/* Zakres sąsiedztwa dla granicznych sekwencji jest nieregularny, dlatego część
 * obliczeń pierwszych 128 sekwencji jest doliczana na cpu (a dokładniej
 * 128 + 127 + ... + 1 wywołań needlman_wunsch */
int* calculate_cpu(const std::vector<Fragment> &seq);
int needleman_wunsch(const char* seq1, const char* seq2);
void print_results(int* cpu_out, int* gpu_out, std::vector<Fragment> const &seq);

#endif //DECLARATOINS_H