#include <cuda_runtime.h>
#include <stdlib.h>
#include <stdio.h>
#include "mandelbrot.h"

#define BLOCK_SIZE 8

__global__ void mandel(double X0, double Y0, double X1, double Y1,
                       int POZ, int PION, int ITER, int *d_Iters) {
                           
    double dX = (X1 - X0) / (POZ - 1);
    double dY = (Y1 - Y0) / (PION-1);
    double x, y, Zx, Zy, tZx;
    int i;
    int pion, poz;

    pion = blockIdx.y * blockDim.y + threadIdx.y;
    poz  = blockIdx.x * blockDim.x + threadIdx.x;

    if(poz < POZ && pion < PION) {
    x = X0 + poz * dX;
    y = Y0 + pion * dY;
    Zx = x;
    Zy = y;
    i  = 0;

    while ( (i<ITER) && ((Zx*Zx+Zy*Zy)<4) ){
        tZx = Zx*Zx-Zy*Zy+x;
        Zy  = 2*Zx*Zy+y;
        Zx  = tZx;

        i++;
    }
    d_Iters[pion * POZ + poz] = i;
   }
}

int main(int argc, char **argv) {
    //Set computation area
    //{X0,Y0} - left bottom corner
    double X0 = atof(argv[1]);
    double Y0 = atof(argv[2]);

    //{X1,Y1} - right upper corner
    double X1=atof(argv[3]);
    double Y1=atof(argv[4]);

    //Set size in pixels
    int POZ  = atoi(argv[5]);
    int PION = atoi(argv[6]);

    //Set number of iterations
    int ITER = atoi(argv[7]);

    //Allocate memory for result array
    int *Iters;
    Iters = (int*) malloc(sizeof(int) * POZ * PION);


    //Do computations
    time_t start, end;

    printf("Computations for rectangle { (%lf %lf), (%lf %lf) }\n", X0, Y0, X1, Y1);
    int *d_Iters;
    start = clock();
    cudaMalloc((void **) &d_Iters, sizeof(int)*POZ*PION);
    end = clock();
    printf("cudaMalloc took %lf s\n\n", 1.0*(end - start)/CLOCKS_PER_SEC);

    start = clock();
    mandel<<<dim3(POZ / BLOCK_SIZE + 1, PION / BLOCK_SIZE + 1), dim3(BLOCK_SIZE, BLOCK_SIZE)>>>(X0, Y0, X1, Y1, POZ, PION, ITER, d_Iters);
    end = clock();
    printf("Computation took %lf s\n\n", 1.0 * (end - start) / CLOCKS_PER_SEC);

    start = clock();
    cudaMemcpy(Iters, d_Iters, sizeof(int) * POZ * PION, cudaMemcpyDeviceToHost);
    end = clock();
    printf("cudaMemcpy took %lf s\n\n", 1.0 * (end - start) / CLOCKS_PER_SEC);

    start = clock();
    makePicture(Iters, POZ, PION, ITER);
    end = clock();
    printf("Saving took %lf s\n\n", 1.0 * (end - start) / CLOCKS_PER_SEC);

    return 0;
}

void makePicture(int *Mandel, int width, int height, int MAX){
    int red_value, green_value, blue_value;

    int MyPalette[41][3]={
        {255,255,255}, //0
        {255,255,255}, //1  not used
        {255,255,255}, //2  not used
        {255,255,255}, //3  not used
        {255,255,255}, //4  not used
        {255,180,255}, //5
        {255,180,255}, //6  not used
        {255,180,255}, //7  not used
        {248,128,240}, //8
        {248,128,240}, //9  not used
        {240,64,224},  //10
        {240,64,224},  //11 not used
        {232,32,208},  //12
        {224,16,192},  //13
        {216,8,176},   //14
        {208,4,160},   //15
        {200,2,144},   //16
        {192,1,128},   //17
        {184,0,112},   //18
        {176,0,96},    //19
        {168,0,80},    //20
        {160,0,64},    //21
        {152,0,48},    //22
        {144,0,32},    //23
        {136,0,16},    //24
        {128,0,0},     //25
        {120,16,0},    //26
        {112,32,0},    //27
        {104,48,0},    //28
        {96,64,0},     //29
        {88,80,0},     //30
        {80,96,0},     //31
        {72,112,0},    //32
        {64,128,0},    //33
        {56,144,0},    //34
        {48,160,0},    //35
        {40,176,0},    //36
        {32,192,0},    //37
        {16,224,0},    //38
        {8,240,0},     //39
        {0,0,0}        //40
    };

    FILE *f = fopen("Mandel2_cu.ppm", "wb");
    fprintf(f, "P6\n%i %i 255\n", width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            //Compute index to the palette
            int indx = (int) floor(5.0 * log2f(1.0f * Mandel[j * width + i] + 1));
            red_value   = MyPalette[indx][0];
            green_value = MyPalette[indx][2];
            blue_value  = MyPalette[indx][1];

            fputc(red_value, f);   // 0 .. 255
            fputc(green_value, f); // 0 .. 255
            fputc(blue_value, f);  // 0 .. 255
        }
    }
    fclose(f);

}


void makePictureInt(int *Mandel,int width, int height, int MAX){
    double scale = 255.0 / MAX;
    int red_value, green_value, blue_value;

    int MyPalette[35][3]={
        {255,0,255},
        {248,0,240},
        {240,0,224},
        {232,0,208},
        {224,0,192},
        {216,0,176},
        {208,0,160},
        {200,0,144},
        {192,0,128},
        {184,0,112},
        {176,0,96},
        {168,0,80},
        {160,0,64},
        {152,0,48},
        {144,0,32},
        {136,0,16},
        {128,0,0},
        {120,16,0},
        {112,32,0},
        {104,48,0},
        {96,64,0},
        {88,80,0},
        {80,96,0},
        {72,112,0},
        {64,128,0},
        {56,144,0},
        {48,160,0},
        {40,176,0},
        {32,192,0},
        {16,224,0},
        {8,240,0},
        {0,0,0}
    };

    FILE *f = fopen("Mandel.ppm", "wb");

    fprintf(f, "P3\n%i %i 255\n", width, height);
    printf("MAX = %d, scale %lf\n", MAX, scale);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++)
        {
            int indx = (int) round(4 * log2((double) Mandel[j * width + i] + 1));
            red_value   = MyPalette[indx][0];
            green_value = MyPalette[indx][2];
            blue_value  = MyPalette[indx][1];

            fprintf(f,"%d ",red_value);   // 0 .. 255
            fprintf(f,"%d ",green_value); // 0 .. 255
            fprintf(f,"%d ",blue_value);  // 0 .. 255
        }
        fprintf(f,"\n");
    }
    fclose(f);
}
