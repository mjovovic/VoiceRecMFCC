#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>
#include <Python.h>
#include <string.h>
#include "matplotlib-cpp/matplotlibcpp.h"
#include "kissfft/kiss_fft.h"
#include "DTW_cpp/include/DTW.hpp"

namespace plt = matplotlibcpp;


typedef struct classP{

    std::vector<std::vector<double>> classPrototype;
    char name[20];

}classP;


int noiseLimit(SF_INFO* fileInfo,  short* sampleVal){
    //100ms for num of samples
    int sampleNum = fileInfo->samplerate/10;
    int midVal = 0;

    for (int i = 0; i < sampleNum; i++) {
        midVal += *(sampleVal+i);
    }
    
    midVal = midVal / sampleNum;
 
    int standardDev = 0;
 
    for (int i = 0; i < sampleNum; i++){
        standardDev += (*(sampleVal+i)-midVal) * (*(sampleVal+i)-midVal);
    }
    
    standardDev = sqrt( standardDev / sampleNum );
    
    return (midVal + 2 * standardDev);
    
}
// c=1 vrati start c=2 vrati kraj reci
int limitsOfWord(int limit, short* items,int size, int c){
    int start = 0, end = 0;
    int startFlag = 0;
    for (int i = 0; i < size; i++) {
        if (startFlag == 0 && limit < abs(*(items+i))) {
            start = i;
            startFlag = 1;
        }
        if (startFlag == 1 && limit < abs(*(items+i))) {
            end = i;
        }
    }
    if(c == 1) {
        return start;
    }
    return end;
}

//plotter
void plottingWord(int limit, short* items, int size) {

    int start = 0, end = 0;
    int startFlag = 0;
    for (int i = 0; i < size; i++) {
        if (startFlag == 0 && limit < abs(*(items+i))) {
            start = i;
            startFlag = 1;
        }
        if (startFlag == 1 && limit < abs(*(items+i))) {
            end = i;
        }
        
    }
    printf("start %d\nfinish %d\nlimit %d", start, end, limit);
    
    if (start < 5 ) {
        printf("JUST NOISE");
    }
    
    std::vector<int> pltVec(items, items + size);

    plt::plot( pltVec );

    plt::axvline(start);
    plt::axvline(end);
    plt::axhline(limit);
    plt::axhline(-limit);

    
    plt::show();

    
}

//konverzija u mono
short* conversionMono(short* items, SF_INFO* info){

    if(info->channels == 0)
        return items;
    short holder;
    int k = 0;
    short * monoItems = (short*)malloc(info->frames * sizeof(short));
    for (int i = 0; i < ( info->frames*info->channels ); i+=info->channels) {        
        holder = 0;
        for ( int j = i; j < i+info->channels; j++){
            holder += items[j];
        }
        monoItems[k] = holder/info->channels;
        k++;
        
    }
    return monoItems;
    free(items);     
}

void oneWindowSpec(SF_INFO* info, short* items){
     //samplerate/10 for 100ms
     
    kiss_fft_cfg cfg = kiss_fft_alloc( info->samplerate/10, 0, NULL, NULL );
    
    kiss_fft_cpx * pInput = (kiss_fft_cpx*)malloc(info->samplerate/10 * sizeof(kiss_fft_cpx));
    kiss_fft_cpx * pOutput = (kiss_fft_cpx*)malloc(info->samplerate/10 * sizeof(kiss_fft_cpx));

    for (int i = 0; i < info->samplerate/10; i++){
        pInput[i].r = *(items+i);
        pInput[i].i = 0;
    }
    

    kiss_fft(cfg, pInput, pOutput);

    std::vector<double> pltVec(items, items + info->samplerate/10);
    for (int i = 0; i < info->samplerate/10; i++){
        pltVec[i] = sqrt(pOutput[i].r*pOutput[i].r + pOutput[i].i*pOutput[i].i);//euklidska da bi se otarasili imaginarnog dela i delimo sa brojem framejvoa za normalizaciju
        
    }

    free(pInput);
    plt::plot( pltVec, 10 );
    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Magnitude");
    plt::show();
}

/*
    printf("%d - %d - %d - %d - %d - %d\n", pMaleInfo->channels, pMaleInfo->frames, pMaleInfo->format, pMaleInfo->samplerate, pMaleInfo->sections, pMaleInfo->seekable);
    printf("%d - %d - %d - %d - %d - %d\n", pFemaleInfo->channels, pFemaleInfo->frames, pFemaleInfo->format, pFemaleInfo->samplerate, pFemaleInfo->sections, pFemaleInfo->seekable);
    printf("%d - %d - %d - %d - %d - %d\n", pGuitarInfo->channels, pGuitarInfo->frames, pGuitarInfo->format, pGuitarInfo->samplerate, pGuitarInfo->sections, pGuitarInfo->seekable);
    printf("%d - %d - %d - %d - %d - %d\n", pWhitenoiseInfo->channels, pWhitenoiseInfo->frames, pWhitenoiseInfo->format, pWhitenoiseInfo->samplerate, pWhitenoiseInfo->sections, pWhitenoiseInfo->seekable);
    2 - 27958 - 65538 - 44100 - 1 - 1
    2 - 54455 - 65538 - 48000 - 1 - 1
    2 - 97374 - 65538 - 44100 - 1 - 1
    1 - 441001 - 65538 - 44100 - 1 - 1
  */

//n vrednost za koju racunamo vrednost, c je centar gde je vrednost 1 prevc je pocetak prozorne i centar prethodnog vrednost 0
double triangularWinFun(int n, int c,int prevc) {
    int start = 0;
    int centar = c - prevc;
    int size = 2*(centar);
    return( 1.0 - fabs( ((2.0* n) - size)/ size));
} 

//vidi gde ces skladistiti izvode da li velika matrica ili nes drugo
void derivativeMFCC(int q, int d, int numOfVec, double** matrixMFFC){
    int t = 2;
    if (numOfVec == 0)
        return;
    double firstDer[q][numOfVec];
    double secDer[q][numOfVec];
    double holder = 0;
    for ( int i = 0; i < numOfVec; i++ ) {
        for ( int j = 0; j < q; j++ ) {
            if(j < 2) {
                firstDer[j][i] = matrixMFFC[j+2][i];
            } else {
                firstDer[j][i] = matrixMFFC[j+2][i] - matrixMFFC[j-2][i];
            }
        }
    }
    if(numOfVec == 2){
        for ( int i = 0; i < numOfVec; i++ ) {
            for ( int j = 0; j < q; j++ ) {
                if(j < 2) {
                    secDer[j][i] = firstDer[j+2][i];
                } else {
                    secDer[j][i] = firstDer[j+2][i] - firstDer[j-2][i];
                }
            }
        }
    }


}


double** MFCC(SF_INFO* info, short* items, int startOfWord, int endOfWord,int q, int width, double** frameSpec,int derDegree){
    
    const int milisecNum = info->samplerate/50;//20 milisekundi
   /// printf("%d\n", milisecNum);
    kiss_fft_cfg cfg = kiss_fft_alloc( milisecNum, 0, NULL, NULL );
    kiss_fft_cpx * pInput = (kiss_fft_cpx*)malloc(milisecNum * sizeof(kiss_fft_cpx));
    kiss_fft_cpx * pOutput = (kiss_fft_cpx*)malloc(milisecNum * sizeof(kiss_fft_cpx));
    int filterNum = 25;
    int maxFreq = info->samplerate/2;
    double maxMelFreq = 1127 * log(1+ maxFreq/700.0);
    double deltaMel =  maxMelFreq / filterNum + 1;
 

    int centerMel[filterNum];
    int centerFreq[filterNum];


    double cepstar[q] = {0};
    /*/koji je kurac ovom mallocu
    double** frameSpec = (double**)malloc(sizeof(double*)*q);
    for (int i = 0; i < q; i++)
        frameSpec[i] = (double*)malloc(sizeof(double)*width);*/
  
    for ( int i = 0; i < filterNum; i++ ) {
        centerMel[i] = (i + 1) * deltaMel;
        centerFreq[i] = 700.0 * ( pow(M_E, ( centerMel[i]/1127.0 )) - 1);
    }

    int freqRes = info->samplerate/(milisecNum);
    int specIndex[filterNum];
    for ( int i = 0 ; i < filterNum; i++ ) { 
        specIndex[i] =  centerFreq[i]/freqRes;
    }
    
    for ( int i = startOfWord; i < endOfWord; i+=milisecNum/2 ) {

        for ( int j = i; j < i+milisecNum; j++ ) {
            pInput[j-i].r = *(items+j);
            pInput[j-i].i = 0;
        }

        kiss_fft(cfg, pInput, pOutput);

        std::vector<double> pltVec(items, items + milisecNum/2);
        for ( int l = 0; l < milisecNum/2; l++ ){
            pltVec[l] = /*sqrt(*/pOutput[l].r*pOutput[l].r + pOutput[l].i*pOutput[l].i/*)*/;//euklidska da bi se otarasili imaginarnog dela (isto radi sto i apsolutna)
           // pltVec[i] = pltVec[i] * pltVec[i];                                      //kvadriramo po formulama da bi dobili isti rez
            //printf("%d\n", milisecNum);
        }
        
        double melCoef[filterNum] = {0};
        double sum;
        for ( int j = 0; j < filterNum; j++ ) {
            sum = 0;
            double triWin;
            if (j == 0){
                for (int k = 0; k < specIndex[j]*2; k++) {
                    triWin=triangularWinFun(k, specIndex[j], 0);
                    sum += triWin * pltVec[k];
                }
            } else {
                for (int k = 0; k < (specIndex[j]-specIndex[j-1])*2; k++) {
                    triWin = triangularWinFun(k , specIndex[j], specIndex[j-1]);
                    sum += triWin * pltVec[k];                       
                    //printf("%f\n", sum);
                }
            }
            melCoef[j] = sum;
        }
        
      //  printf("\n-------------------\n");
        for (int j = 0; j < q; j++) {
            sum = 0;
            for (int k = 0; k < filterNum; k++) {
                sum += log(melCoef[k]) * cos( (M_PI * j *(k*2.0 + 1) )/(2.0*filterNum ));
            }
            cepstar[j] = sum;
          //  printf("%d  %d \n",j, (i-startOfWord)/(milisecNum/2));

            frameSpec[j][(i-startOfWord)/(milisecNum/2)] = sum;    
            
          //  printf("[%f] ",frameSpec[j][(i-startOfWord)/(milisecNum/2)]);

        }
        //IZVODI U ZAVISNOSTI OD INPUTA KORISNIKA
        if(derDegree > 1){
            for ( int j = q; j <q*2; j++ ) {
                if((j-q) < 2) {
                    frameSpec[j][(i-startOfWord)/(milisecNum/2)] = frameSpec[j+2-q][(i-startOfWord)/(milisecNum/2)];
                } else if((j-q) > q-2) {
                    frameSpec[j][(i-startOfWord)/(milisecNum/2)] = -frameSpec[j-2-q][(i-startOfWord)/(milisecNum/2)];
                } else {
                    frameSpec[j][(i-startOfWord)/(milisecNum/2)] = frameSpec[j+2-q][(i-startOfWord)/(milisecNum/2)] - frameSpec[j-2-q][(i-startOfWord)/(milisecNum/2)];
                }
            }    
        }
        if(derDegree == 3){
            for ( int j = q*2; j <q*3; j++ ) {
                if((j-(q*2)) < 2) {
                    frameSpec[j][(i-startOfWord)/(milisecNum/2)] = frameSpec[j+2-q][(i-startOfWord)/(milisecNum/2)];
                } else if((j-(q*2)) > q-2) {
                    frameSpec[j][(i-startOfWord)/(milisecNum/2)] = -frameSpec[j-2-q][(i-startOfWord)/(milisecNum/2)];
                } else {
                    frameSpec[j][(i-startOfWord)/(milisecNum/2)] = frameSpec[j+2-q][(i-startOfWord)/(milisecNum/2)] - frameSpec[j-2-q][(i-startOfWord)/(milisecNum/2)];
                }
            }
        }
    }
    return frameSpec;
       
}



std::vector<std::vector<double>> vectorMFCC(char* filename, int derDegree){
  
    SF_INFO * pInfo = (SF_INFO*)malloc( sizeof(SF_INFO) );
  //  SNDFILE * pFile = sf_open("./samples/v-please.wav", SFM_READ, pInfo );
    SNDFILE * pFile = sf_open(filename, SFM_READ, pInfo );
    short * pItems = (short*)malloc(pInfo->frames*pInfo->channels * sizeof(short));
    sf_readf_short( pFile, pItems, pInfo->frames);

    pItems = conversionMono(pItems, pInfo);
    int limit = noiseLimit(pInfo, pItems);
    int startOfWord = limitsOfWord(limit, pItems, pInfo->frames, 1);
    int endOfWord = limitsOfWord(limit, pItems, pInfo->frames, 2);
    int q = 13;
    int windowSize = pInfo->samplerate/50;       //20ms
    int width = (endOfWord-startOfWord)/(windowSize/2);
    
       
   // plottingWord(limit, pItems, pInfo->frames);

    double** matrixMFFC = (double**)malloc(sizeof(double)*q*derDegree);
    for (int i = 0; i < q*derDegree; i++)
        matrixMFFC[i] = (double*)malloc(sizeof(double)*width);
    matrixMFFC = MFCC(pInfo, pItems, startOfWord, endOfWord, q, width, matrixMFFC, derDegree);
    
    ///Creating vector for DTW
    std::vector<std::vector<double>> vectorMFFC;
    for (int i = 0; i < width; i++){
        std::vector<double> tmp ( q*derDegree );
        for (int j = 0; j < q*derDegree; j++){
           // printf("%d\n", j);
            tmp[j] = matrixMFFC[j][i];
        }
        vectorMFFC.push_back(tmp);
       // printf("%d - %d\n", i, width);
    }
    free(matrixMFFC);
    return vectorMFFC;
}

std::vector<std::vector<double>> DWTprototype( std::vector< std::vector<double> > a, std::vector< std::vector<double> > b, int height ){

    double p = 2;  
    
    DTW::DTW MyDtw (a, b, p);
    
    std::vector<std::vector<int> > path;
    path = MyDtw.path();
    std::vector<std::vector<double>> prototype;
    for ( int i = 0; i < path.size(); i++ ) {
        std::vector<double> tmp (height);
        for( int j = 0; j  < height; j++) { 
            tmp[j] = (a[path[i][0]][j] + b[path[i][1]][j])/2.0;
        //     printf("[%f] ", tmp[j]);
        }
     //  printf("\n");
       prototype.push_back(tmp);
    }

    return prototype;
}


int main(int argc, char**argv){
    
    int derDegree = 1;
    int q = 13;
    
    char name1[30];
    char name2[30];
    char name3[30];
    strcpy(name1, "srsly");
    strcpy(name2, "welcome");
    strcpy(name3, "please");

    std::vector<char*> classNames = {name1, name2, name3};
    
    // strcpy(name1, "./samples/v-srsly.wav");
    // strcpy(name2, "./samples/v-please.wav");
    // strcpy(name3, "./samples/v-welcome.wav");
    // std::vector<char*> testNames = {name1, name2, name3};

    
    // for (int i = 0; i < classNames.size(); i++) {
    //     //odvratno pravljenje imena za test fajlove
        
    //     char fFileName[30] = {0}; 
    //     char mFileName[30] = {0};
    //     strcpy(mFileName, "./samples/");
    //     strcat(mFileName, "m-");
    //     strcat(mFileName, classNames[i]);
    //     strcat(mFileName, ".wav");
    //     std::vector<std::vector<double>> male = vectorMFCC(mFileName, derDegree);

    //     strcpy(fFileName, "./samples/");
    //     strcat(fFileName, "f-");
    //     strcat(fFileName, classNames[i]);
    //     strcat(fFileName, ".wav");
    //     std::vector<std::vector<double>> female = vectorMFCC(fFileName, derDegree);

    //     std::vector<std::vector<double>> prototype = DWTprototype(male, female, (q*derDegree));
    //    // printf("%s\t%s\n", fFileName, mFileName);

    // }
    
    char fileName[30];
    strcpy(fileName,"./samples/m-welcome.wav" );
    std::vector<std::vector<double>> a = vectorMFCC(fileName, derDegree);
    strcpy(fileName,"./samples/f-welcome.wav" );
    std::vector<std::vector<double>> b = vectorMFCC(fileName, derDegree);
    strcpy(fileName,"./samples/v-welcome.wav" );
    std::vector<std::vector<double>> c = vectorMFCC(fileName, derDegree);
    
    std::vector<std::vector<double>> prototype = DWTprototype(a, b, (q*derDegree));

    std::cout << "DTW distance: " << DTW::dtw_distance_only(a, b, 2.0) << std::endl;
    
    // char fileName2[30];
    // strcpy(fileName2,"./samples/m-srsly.wav" );
    // std::vector<std::vector<double>> a1 = vectorMFCC(fileName2, derDegree);
    // strcpy(fileName2,"./samples/f-srsly.wav" );
    // std::vector<std::vector<double>> b1 = vectorMFCC(fileName2, derDegree);
    // strcpy(fileName2,"./samples/v-srsly.wav" );
    // std::vector<std::vector<double>> c1 = vectorMFCC(fileName2, derDegree);



    // std::cout << "DTW distance: " << DTW::dtw_distance_only(prototype, c, 2.0) << std::endl;
    // std::cout << "DTW distance: " << DTW::dtw_distance_only(b, c, 2.0) << std::endl;




    return 0;

}

