#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include "declarations.h"

#define MAX_BLOCKS 65535

/* Pretwarzamy macierz wiersz po wierszu. Dla ostatnich 128 sekwencji doklejamy
 * wcześniejsze, tak żeby w każdym bloku wszystkie 128 wątków pracowało */
__global__ void needlman(char* seq, int* d_out, int count, long size) {
    int prev_row[SINGLE_SEQ_LENGTH];
    int current, prev;
    int max = 0;
    int block = blockIdx.x + (MAX_BLOCKS * count);
    int idx = block * blockDim.x + threadIdx.x;
    int idx2 = block + threadIdx.x + 1;

    //doklejamy poprzednie sekwencje
    if(idx2 >= size) {
        idx2 = size - threadIdx.x - NEIGHBOURS - 1 - 1;
    }

    extern __shared__ char seq1[];
    if(threadIdx.x < SINGLE_SEQ_LENGTH)
        seq1[threadIdx.x] = seq[block * SINGLE_SEQ_LENGTH + threadIdx.x];
        

    char seq2[SINGLE_SEQ_LENGTH];

    for(int i = 0; i < SINGLE_SEQ_LENGTH; i++)
        seq2[i] = seq[(idx2) * SINGLE_SEQ_LENGTH + i];

    for(int i = 0; i < SINGLE_SEQ_LENGTH; i++)
        prev_row[i] = 0;

    for(int i = 0; i < SINGLE_SEQ_LENGTH; i++) {
        for(int j = 0; j < SINGLE_SEQ_LENGTH; j++) {
            current = (seq1[i] == seq2[j]) ? 1 : -2;

            if(j > 0) {
                current = (current + prev_row[j - 1] >= prev_row[j] - GP &&
                    current + prev_row[j - 1] >= prev - GP) ? current + prev_row[j - 1] :
                    prev_row[j] >= prev ? prev_row[j] - GP : prev - GP;

                prev_row[j - 1] = prev;
                if(j == SINGLE_SEQ_LENGTH - 1)
                    prev_row[j] = current;
            } else
                current = (current >= prev_row[j] - GP && current >= -GP) ? current :
                    prev_row[j] >= 0 ? prev_row[j] - GP :  -GP;

            prev = current;
        }
        if(current >= max)
            max = current;
    }

    for(int i = 0; i < SINGLE_SEQ_LENGTH; i++)
        if(prev_row[i] >= max)
            max = prev_row[i];

    d_out[idx] = max;
}

int* calculate_gpu(const std::vector<Fragment> &seq, char* sequences) {
    long  length     = 128 * seq.size();
    int*  gpu_out    = (int*) malloc(length * sizeof(int));
    long  seq_length = ((seq.size() * SINGLE_SEQ_LENGTH) + 1);
    char* d_seq;
    int*  d_out;

    memset(gpu_out, 0, length * sizeof(int));
    
    cudaMalloc(&d_seq, seq_length);
    cudaMalloc(&d_out, length * sizeof(int));
    cudaMemset(d_out, 0, length * sizeof(int));
    cudaMemcpy(d_seq, sequences, seq_length, cudaMemcpyHostToDevice);

    int count = 0;
    long size = seq.size();
    do {
        long grid_size = std::min((long) MAX_BLOCKS, size);
        needlman<<<grid_size, 128, SINGLE_SEQ_LENGTH>>>(d_seq, d_out, count, seq.size());
        cudaDeviceSynchronize();
        count++;
        size -= grid_size;
    } while(size > 0);

    cudaMemcpy(gpu_out, d_out, length * sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(d_seq);
    cudaFree(d_out);

    return gpu_out;
}



int main() {
    std::vector<Fragment> seq = prepare_set();
    char* sequences = prepare_seq(seq);
    int*  gpu_out = calculate_gpu(seq, sequences);
    int*  cpu_out = calculate_cpu(seq);

    print_results(cpu_out, gpu_out, seq);

    free(cpu_out);
    free(gpu_out);
    free(sequences);


    return 0;
}