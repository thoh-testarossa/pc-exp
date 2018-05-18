#include <iostream>
#include <cmath>
#include "mpi.h"

#define SCALE 2000000000
#define BLOCKSIZE SCALE/1000

char countZero(char x)
{
    x = (x & 0x55) + ((x >> 1) & 0x55);
    x = (x & 0x33) + ((x >> 2) & 0x33);
    x = (x & 0xf) + ((x >> 4) & 0xf);
    return 8 - x;
}

bool oddDivisibleJudge(int numerator, int denominator)
{
    while(numerator >= denominator)
    {
        while((numerator & 1) == 0) numerator >>= 1;
        if(numerator >= denominator) numerator -= denominator;
    }

    return numerator == 0;
}

bool checkBit(unsigned char byte, int posBit)
{
    char chkByte = (unsigned char) (1 << posBit);
    return (chkByte & byte) != 0;
}

int fromBitToNum(int posArray, int posBit)
{
    return (posArray << 4) + (posBit << 1) + 1;
}

void generateTable(unsigned char *bitArray, int upperBorder)
{
    for(int i = 3; i < upperBorder; i += 2)
    {
        int i_chk = (int)sqrt(i);
        int i_chk_posArray = i_chk >> 4;
        bool isPrime = true;
        for(int j = 0; j <= i_chk_posArray; j++)
        {
            if(fromBitToNum(j, 0) > i_chk) break;
            else if(checkBit(bitArray[j], 0) && oddDivisibleJudge(i, fromBitToNum(j, 0))) {isPrime = false; break;}
            if(fromBitToNum(j, 1) > i_chk) break;
            else if(checkBit(bitArray[j], 1) && oddDivisibleJudge(i, fromBitToNum(j, 1))) {isPrime = false; break;}
            if(fromBitToNum(j, 2) > i_chk) break;
            else if(checkBit(bitArray[j], 2) && oddDivisibleJudge(i, fromBitToNum(j, 2))) {isPrime = false; break;}
            if(fromBitToNum(j, 3) > i_chk) break;
            else if(checkBit(bitArray[j], 3) && oddDivisibleJudge(i, fromBitToNum(j, 3))) {isPrime = false; break;}
            if(fromBitToNum(j, 4) > i_chk) break;
            else if(checkBit(bitArray[j], 4) && oddDivisibleJudge(i, fromBitToNum(j, 4))) {isPrime = false; break;}
            if(fromBitToNum(j, 5) > i_chk) break;
            else if(checkBit(bitArray[j], 5) && oddDivisibleJudge(i, fromBitToNum(j, 5))) {isPrime = false; break;}
            if(fromBitToNum(j, 6) > i_chk) break;
            else if(checkBit(bitArray[j], 6) && oddDivisibleJudge(i, fromBitToNum(j, 6))) {isPrime = false; break;}
            if(fromBitToNum(j, 7) > i_chk) break;
            else if(checkBit(bitArray[j], 7) && oddDivisibleJudge(i, fromBitToNum(j, 7))) {isPrime = false; break;}
        }
        if(isPrime) bitArray[i >> 4] |= (unsigned char)(1 << ((i & 0xf) >> 1));
    }
}

int judgeAllFromPrimeTable(int from, int to, int *primeTable, int pCount)
{
    unsigned char *notPrimeBitArray = new unsigned char [(to - from) >> 4];
    for(int i = 0; i < (to - from) >> 4; i++) notPrimeBitArray[i] = 0;

    notPrimeBitArray[0] |= 1;

    for(int i = 0; i < pCount && primeTable[i] * primeTable[i] < to; i++)
    {
        int startNum = from - (from % primeTable[i]) + primeTable[i];
        if((startNum & 1) == 0) startNum += primeTable[i];
        if(startNum == primeTable[i]) startNum += 2 * primeTable[i];

        for(int j = startNum; j < to; j += 2 * primeTable[i])
        {
            int index = j - from;
            notPrimeBitArray[index >> 4] |= 1 << (((index & 0xf) >> 1));
        }
    }
    int ans = 0;
    for(int i = 0; i < (to - from) >> 4; i++)
        ans += countZero(notPrimeBitArray[i]);
    return ans;
}

int main(int argv, char *argc[])
{
    int thread_num = 4;
    if(argv > 1) thread_num = atoi(argc[1]);

    int factorUpperBorder = (int) sqrt(SCALE);
    int arraySize = (factorUpperBorder >> 4) + 1;

    unsigned char primeTable[arraySize];
    for (int i = 0; i < arraySize; i++) primeTable[i] = 0;

    int pCount = 0;

    generateTable(primeTable, factorUpperBorder);

    for (int i = 0; i < arraySize; i++)
    {
        for (int j = 0; j < 8; j++)
            if (checkBit(primeTable[i], j))
                pCount++;
    }

    int *primeArray = new int [pCount];

    for (int i = 0, k = 0; i < arraySize; i++)
    {
        for (int j = 0; j < 8; j++)
            if (checkBit(primeTable[i], j))
                primeArray[k++] = fromBitToNum(i, j);
    }

    int sum = 0, mysum = 0;

    int numprocs, myid;

    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    for(int from = BLOCKSIZE * myid; from < SCALE; from += BLOCKSIZE * numprocs)
    {
        int to = from + BLOCKSIZE;
        if(to > SCALE) to = SCALE;
        mysum += judgeAllFromPrimeTable(from, to, primeArray, pCount);
    }

    MPI_Reduce(&mysum, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    std::cout << myid << " " << mysum << std::endl;

    if(0 == myid)
        std::cout << sum + 1 << std::endl;

    MPI_Finalize();

    return 0;
}
