#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>  // For size_t.
#include <stdint.h>  // For uint8_t, uint32_t.

typedef uint32_t word;  // A word is an unsigned of 32-bit size.
typedef word block[4];  // A block is an array of 4 words.

// The standard Rijndael S-box
const static word SBOX[256] = {
        0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B,
        0xFE, 0xD7, 0xAB, 0x76, 0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0,
        0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0, 0xB7, 0xFD, 0x93, 0x26,
        0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
        0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07, 0x12, 0x80, 0xE2,
        0xEB, 0x27, 0xB2, 0x75, 0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0,
        0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84, 0x53, 0xD1, 0x00, 0xED,
        0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
        0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F,
        0x50, 0x3C, 0x9F, 0xA8, 0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5,
        0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2, 0xCD, 0x0C, 0x13, 0xEC,
        0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
        0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE, 0xB8, 0x14,
        0xDE, 0x5E, 0x0B, 0xDB, 0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C,
        0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79, 0xE7, 0xC8, 0x37, 0x6D,
        0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
        0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F,
        0x4B, 0xBD, 0x8B, 0x8A, 0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E,
        0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E, 0xE1, 0xF8, 0x98, 0x11,
        0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
        0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F,
        0xB0, 0x54, 0xBB, 0x16};

//Compute Hamming Weight of an 8-bit number
int weight(unsigned char x)
{
    //Compute using bitshift operations instead of loops
    return ((0x876543210 >>
            (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >>
            ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2)) & 0xf;
}

uint8_t T[600*55];
uint8_t H[600*256];

//Get value from H
//Line and entry have to be passed as indexes starting from 1 (like in lecture slides)
uint8_t getFromH(int line, int entry) {
    return H[(line-1)*256 + (entry - 1)];
}

//Get value from T
//Line and entry have to be passed as indexes starting from 1 (like in lecture slides)
uint8_t getFromT(int line, int entry) {
    return T[(line-1)*55 + (entry - 1)];
}

//Function to calculate correlation

double corr(int j,int l){
    int tempSum = 0;
    //Calculate the average h for column j
    for(int i=1;i<=55;i++){
        tempSum += getFromH(i,j);
    }
    double averageH = (double)tempSum/55.0;
    //Calculate the average t for l-th column
    tempSum = 0;
    for(int i=1;i<=55;i++){
        tempSum += getFromH(i,j);
    }
    double averageT = (double)tempSum/55.0;
    double numerator = 0;
    double denominator1 = 0;
    double denominator2 = 0;
    double temp1,temp2;
    //Calculate numerator and two denominator terms following the formula from slides
    for(int i=1;i<=600;i++){
        temp1 = (double)getFromH(i,j)-averageH;
        temp2 = (double)getFromT(i,l)-averageT;
        numerator += temp1*temp2;
        denominator1 += temp1 * temp1;
        denominator2 += temp2 * temp2;
    }

    return (numerator/sqrt(denominator1*denominator2));
}


#define MAX_LEN 1024
int main() {

    //Initialization of random number generator
    srand(time(NULL));

    //Load T table 600x55 from T2.dat file but store it as a single-dim array

    FILE* fp;
    fp = fopen("T0.dat", "r");
    if (fp == NULL) {
        perror("Failed: ");
        return 1;
    }

    char buffer[MAX_LEN];
    // -1 to allow room for NULL terminator for really long string
    int i=0;
    while (fgets(buffer, MAX_LEN - 1, fp))
    {
        // Remove trailing newline
        buffer[strcspn(buffer, "\n")] = 0;
        //Tokenize the line that was read
        char * pch;
        pch = strtok (buffer,",");
        while (pch != NULL)
        {
            //Add 0.5 to round float to nearest integer
            T[i] = (int)((strtof(pch,NULL) + 0.5f));
            //printf("%d%s%d\n",i,"th element is: ", T[i]);
            i++;
            pch = strtok (NULL, ",");
        }
    }
    fclose(fp);

    //Load P 600 entries array and store it as a single-dim array

    uint8_t P[600];

    fp = fopen("inputs0.dat", "r");
    if (fp == NULL) {
        perror("Failed: ");
        return 1;
    }
    char buffer2[3000];
    // -1 to allow room for NULL terminator for really long string
    i=0;
    while (fgets(buffer2, 3000, fp))
    {
        //Tokenize the line that was read
        char * pch;
        pch = strtok (buffer2,",");
        while (pch != NULL)
        {
            P[i] = atoi(pch);
            i++;
            //printf("%d%s%d\n",i,"th element is: ", P[i-1]);
            pch = strtok (NULL, ",");
        }
    }
    fclose(fp);

    //Construct H table of size N (600) by 256

    for(i=0;i<600;i++){
        for(int j=1;j<=256;j++) {
            //Calculate each H table entry
            H[(i*256)+(j-1)] = weight(SBOX[(P[i] ^ (j))]);
        }
    }

    //Calculate corr values for each key guess and each individual trace
    double corr_values[256];
    double max_prob = -10.0;    //-10 chosen arbitrarily to be lower than -1 (lowest probability)
    int closest_guess;
    for(int j=1;j<=256;j++){
        for(int l=1;l<=55;l++){
            corr_values[j-1] += corr(j,l);
            if(corr_values[j-1] > max_prob){
                max_prob = corr_values[j-1];
                closest_guess = j;
            }
        }
    }
    //Calculate the average correlation across all power entries
    max_prob /= 256;

    printf("With biggest correlation equal to: %f the most likely guess is %d",max_prob,closest_guess);


    printf("\nPress enter to continue...\n");
    getchar();

    return 0;
}
