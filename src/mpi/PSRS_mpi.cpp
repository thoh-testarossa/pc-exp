//
// Created by Thoh Testarossa on 2018/6/4.
//

#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <cstdlib>
#include "mpi.h"

using namespace std;

#define MASTER_PROC 0

#define DATA_BOUND 1000000000
#define DEFAULT_INPUT_FILE "input.txt"
#define DEFAULT_OUTPUT_FILE "output.txt"

//#define TEST_PART

typedef struct minHeapNode
{
    int value;
    int source;
}minHeapNode;

/*
 * Notice: The minHeap begin at element 1
 */
void adjustMinHeap(minHeapNode *minHeap, int pos, int heapSize)
{
    int i, minpos = pos;
    minHeapNode tmpNode;
    do{
        i = minpos;
        if(2 * i <= heapSize && minHeap[minpos].value > minHeap[2 * i].value) minpos = 2 * i;
        if(2 * i + 1 <= heapSize && minHeap[minpos].value > minHeap[2 * i + 1].value) minpos = 2 * i + 1;
        if(i != minpos)
        {
            tmpNode = minHeap[i];
            minHeap[i] = minHeap[minpos];
            minHeap[minpos] = tmpNode;
        }
    }while(i != minpos);
}

/*
 * The function which build a minHeap from raw data
 */
void buildMinHeap(minHeapNode *minHeap, int heapSize)
{
    for(int i = heapSize / 2; i >= 1; i--) adjustMinHeap(minHeap, i, heapSize);
}

/*
 * The function which delete the minimum element from the minHeap (this operation changes the size of the heap)
 */
void deleteHeapElement(minHeapNode *minHeap, int &heapSize)
{
    if(heapSize <= 0);
    else
    {
        minHeap[1] = minHeap[heapSize--];
        adjustMinHeap(minHeap, 1, heapSize);
    }
}

/*
 * Check if the input parameters is legal
 * Only two parameters: the size of dataset, and the way to get the dataset needed.
*/
bool isParameterLegal(int argv, char *argc[])
{
    return (argv == 3) && atoi(argc[1]) > 0 && (string(argc[2]) == string("r") || string(argc[2]) == string("g"));
}

/*
 * Display help message
 */
void help()
{
    cout << "Usage:" << endl
         << "   mpirun -np <num_of_proc> <name_of_exec> <dataset_size> <r|g>" << endl
         << "   r: Read the dataset from a exist file" << endl
         << "   g: Generate a new dataset" << endl
         ;
}

/*
 * The value this function returned is the real size of the dataset
 * This function is used to generate a dataset and store it into a file
 */
int generateDataset(int *dataset, int size, const string &fileName)
{
    int realTotalSize = 0, tmp;

    std::random_device r;
    std::uniform_int_distribution<int> uniform_dest(0, DATA_BOUND);
    std::default_random_engine e1(r());

    ofstream fout(fileName);
    if(fout.is_open())
    {
        for (; realTotalSize < size; realTotalSize++)
        {
            tmp = uniform_dest(e1);
            fout << tmp << " ";
            dataset[realTotalSize] = tmp;
        }
    }

    return realTotalSize;
}

/*
 * The value this function returned is the real size of the dataset
 * This function is used to read a exist dataset from a available file
 */
int readDatasetFromFile(int *dataset, int size, const string &fileName)
{
    int realTotalSize = 0;

    ifstream fin(fileName);
    if(fin.is_open())
        for(; realTotalSize < size && !fin.eof(); realTotalSize++)
            fin >> dataset[realTotalSize];

    return realTotalSize;
}

/*
 * The value this function returned is the real size of the dataset
 */
int getDataset(int *dataset, int size, char method)
{
    int realTotalSize = 0;
    if(method == 'r') realTotalSize = readDatasetFromFile(dataset, size, string(DEFAULT_INPUT_FILE));
    else if(method == 'g') realTotalSize = generateDataset(dataset, size, string(DEFAULT_INPUT_FILE));
    return realTotalSize;
}

/*
 * This function is used for partial sort in each process
 */
void quicksort(int *array, int spos, int epos)
{
    if(spos >= epos);
    else
    {
        int c = spos, d = epos, tmp;
        while(c < d)
        {
            while (c < d && array[c] <= array[d]) c++;
            tmp = array[c];
            array[c] = array[d];
            array[d] = tmp;
            while (c < d && array[c] <= array[d]) d--;
            tmp = array[c];
            array[c] = array[d];
            array[d] = tmp;
        }
        quicksort(array, spos, c - 1);
        quicksort(array, d + 1, epos);
    }
}

int main(int argv, char *argc[])
{
    //Step 1: Initialization
    MPI_Init(nullptr, nullptr);

    int numprocs, myid;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    MPI_Status s;

    int *totalDataset = nullptr, inputTotalSize = 0, realTotalSize = 0;
    char datasetGetMethod;
    if(isParameterLegal(argv, argc))
    {
        if(MASTER_PROC == myid)
        {
            inputTotalSize = atoi(argc[1]);
            datasetGetMethod = argc[2][0];
            totalDataset = new int[inputTotalSize];
            realTotalSize = getDataset(totalDataset, inputTotalSize, datasetGetMethod);
        }
        MPI_Bcast(&realTotalSize, 1, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
    }
    else
    {
        if(MASTER_PROC == myid)
            help();
        MPI_Finalize();
        return -1;
    }

#ifdef TEST_PART
    cout << myid << ": Step 1 end" << endl;
#endif

    //Step 2: Scatter data, local sort and regular samples collected
    //Scatter part
    int mySize = ((myid + 1) * realTotalSize) / numprocs - (myid * realTotalSize) / numprocs;
    int *myDataset = new int [mySize];

    if(MASTER_PROC == myid)
    {
        for(int i = 0; i < mySize; i++) myDataset[i] = totalDataset[i];

        int sendSize, sendStartPos;

        for(int i = 1; i < numprocs; i++)
        {
            sendStartPos = (i * realTotalSize) / numprocs;
            sendSize = ((i + 1) * realTotalSize) / numprocs - (i * realTotalSize) / numprocs;

            MPI_Send(&totalDataset[sendStartPos], sendSize, MPI_INT, i, i, MPI_COMM_WORLD);
        }
    }
    else
        MPI_Recv(myDataset, mySize, MPI_INT, MASTER_PROC, myid, MPI_COMM_WORLD, &s);

    //Partial sort part
    quicksort(myDataset, 0, mySize - 1);

    //Regular sampling part
    int *mySamplingPoint = new int [numprocs];
    for(int i = 0; i < numprocs; i++)
        mySamplingPoint[i] = myDataset[(i * mySize) / numprocs];

#ifdef TEST_PART
    cout << myid << ": Step 2 end" << endl;
#endif

    //Step 3: Gather and merge samples, choose and broadcast p - 1 pivots
    int *totalSamplingPoint = nullptr;
    if(MASTER_PROC == myid)
    {
        totalSamplingPoint = new int [numprocs * numprocs];

        for(int i = 0; i < numprocs; i++)
            totalSamplingPoint[i] = mySamplingPoint[i];
        for(int i = 1; i < numprocs; i++)
            MPI_Recv(&totalSamplingPoint[i * numprocs], numprocs, MPI_INT, i, i, MPI_COMM_WORLD, &s);
    }
    else
        MPI_Send(mySamplingPoint, numprocs, MPI_INT, MASTER_PROC, myid, MPI_COMM_WORLD);

#ifdef TEST_PART
    cout << myid << ": Gather sample part end" << endl;
#endif

#ifdef TEST_PART
    if(MASTER_PROC == myid)
    {
        for(int i = 0; i < numprocs * numprocs; i++) cout << totalSamplingPoint[i] << " ";
        cout << endl << "==========================" << endl;
    }
#endif

    //Multimerge part using minHeap
    int *sortedTotalSamplingPoint = nullptr;
    if(MASTER_PROC == myid)
    {
        int *samplingPointArrayPointer = new int [numprocs];
        for(int i = 0; i < numprocs; i++) samplingPointArrayPointer[i] = i * numprocs;

        minHeapNode *minHeap = new minHeapNode [numprocs + 1];
        int heapSize = numprocs;
        for(int i = 1; i <= numprocs; i++)
        {
            minHeap[i].source = i - 1;
            minHeap[i].value = totalSamplingPoint[samplingPointArrayPointer[i - 1]];
        }

        buildMinHeap(minHeap, heapSize);

#ifdef TEST_PART
        for(int i = 1; i <= heapSize; i++) cout << minHeap[i].source << ": " << minHeap[i].value << endl;
#endif

        sortedTotalSamplingPoint = new int [numprocs * numprocs];
        int pt_sTSP = 0;
        while(heapSize > 0)
        {
            //Fetch the top element of the heap
            sortedTotalSamplingPoint[pt_sTSP++] = minHeap[1].value;
            //Replace the element from the same sorted table if possible
            samplingPointArrayPointer[minHeap[1].source]++;
            //If the corresponding sorted table has some elements left
            if(samplingPointArrayPointer[minHeap[1].source] < (minHeap[1].source + 1) * numprocs)
            {
#ifdef TEST_PART
                cout << "replace" << endl;
#endif
                minHeap[1].value = totalSamplingPoint[samplingPointArrayPointer[minHeap[1].source]];
                adjustMinHeap(minHeap, 1, heapSize);
            }
            //If the corresponding sorted table doesn't have any element left
            else
            {
#ifdef TEST_PART
                cout << "delete" << endl;
#endif
                deleteHeapElement(minHeap, heapSize);
            }

#ifdef TEST_PART
            cout << heapSize << endl;
            for(int i = 0; i < numprocs; i++) cout << samplingPointArrayPointer[i] << " ";
            cout << endl;
            for(int i = 0; i < numprocs * numprocs; i++) cout << sortedTotalSamplingPoint[i] << " ";
            cout << endl;
            for(int i = 1; i <= heapSize; i++) cout << minHeap[i].source << ": " << minHeap[i].value << endl;
#endif
        }
    }

    //Distribute p - 1 pivots
    int *selectedPivots = new int [numprocs];
    if(MASTER_PROC == myid)
    {
        for(int i = 1; i < numprocs; i++)
            selectedPivots[i] = sortedTotalSamplingPoint[i * numprocs];
    }
    MPI_Bcast(selectedPivots, numprocs, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);

#ifdef TEST_PART
    cout << myid << ": Step 3 end" << endl;
#endif

    //Step 4: Local data is partitioned
    //The num of ith class in this process generated
    int *posOfEachPartStartFrom = new int [numprocs];
    posOfEachPartStartFrom[0] = 0;
    int pivotId = 1;
    for(int i = 0; i < mySize; i++)
        while(myDataset[i] >= selectedPivots[pivotId] && pivotId < numprocs) posOfEachPartStartFrom[pivotId++] = i;

    int *numOfEachPart = new int [numprocs];
    for(int i = 0; i < numprocs; i++)
    {
        if(i == numprocs - 1) numOfEachPart[i] = mySize - posOfEachPartStartFrom[i];
        else numOfEachPart[i] = posOfEachPartStartFrom[i + 1] - posOfEachPartStartFrom[i];
    }

#ifdef TEST_PART
    cout << myid << ": Step 4 end" << endl;
#endif

    //Step 5: All *ith* classes are gathered and merged
    //The num of ith class in this process sent to other processes
    int *numOfCorrespondingPartFromEachProcess = new int [numprocs];
    for(int p_send = 0; p_send < numprocs; p_send++)
    {
        for(int p_recv = 0; p_recv < numprocs; p_recv++)
        {
            if(p_send == myid)
            {
                if(p_send != p_recv)
                    MPI_Send(&numOfEachPart[p_recv], 1, MPI_INT, p_recv, p_send, MPI_COMM_WORLD);
                else
                    numOfCorrespondingPartFromEachProcess[p_recv] = numOfEachPart[p_send];
            }
            if(p_recv == myid)
            {
                if(p_recv != p_send)
                    MPI_Recv(&numOfCorrespondingPartFromEachProcess[p_send], 1, MPI_INT, p_send, p_send, MPI_COMM_WORLD, &s);
                else;
            }
        }
    }
    int *startPosOfCorrespondingPartFromEachProcess = new int [numprocs];
    startPosOfCorrespondingPartFromEachProcess[0] = 0;
    for(int i = 1; i < numprocs; i++)
        startPosOfCorrespondingPartFromEachProcess[i] = startPosOfCorrespondingPartFromEachProcess[i - 1] + numOfCorrespondingPartFromEachProcess[i - 1];
    int correspondingPartSize = startPosOfCorrespondingPartFromEachProcess[numprocs - 1] + numOfCorrespondingPartFromEachProcess[numprocs - 1];

    //The corresponding part of data will be delivered to the corresponding process
    int *correspondingPartDataSet = new int [correspondingPartSize];
    for(int p_send = 0; p_send < numprocs; p_send++)
    {
        for(int p_recv = 0; p_recv < numprocs; p_recv++)
        {
            if(p_send == myid)
            {
                if (p_send != p_recv)
                    MPI_Send(&myDataset[posOfEachPartStartFrom[p_recv]], numOfEachPart[p_recv], MPI_INT, p_recv, p_send, MPI_COMM_WORLD);
                else
                {
                    for (int i = posOfEachPartStartFrom[p_recv]; i < posOfEachPartStartFrom[p_recv] + numOfEachPart[p_recv]; i++)
                    {
                        int pos = i - posOfEachPartStartFrom[p_recv];
                        correspondingPartDataSet[startPosOfCorrespondingPartFromEachProcess[p_send] + pos] = myDataset[i];
                    }
                }
            }
            if(p_recv == myid)
            {
                if(p_recv != p_send)
                    MPI_Recv(&correspondingPartDataSet[startPosOfCorrespondingPartFromEachProcess[p_send]], numOfCorrespondingPartFromEachProcess[p_send], MPI_INT, p_send, p_send, MPI_COMM_WORLD, &s);
                else;
            }
        }
    }

    //Every processes quicksort its corresponding part
    quicksort(correspondingPartDataSet, 0, correspondingPartSize - 1);

#ifdef TEST_PART
    cout << myid << ": Step 5 end" << endl;
#endif

    //Step 6: Root processor collects all the data
    //Collect the summary information of each part
    int *correspondingPartSizeArray = nullptr;
    int *startPosForSortedEachProcessPart = nullptr;
    int ret_totalSize = 0;
    if(MASTER_PROC == myid)
    {
        correspondingPartSizeArray = new int [numprocs];
        for(int i = 0; i < numprocs; i++)
        {
            if(i == MASTER_PROC) correspondingPartSizeArray[i] = correspondingPartSize;
            else
                MPI_Recv(&correspondingPartSizeArray[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &s);
        }

        startPosForSortedEachProcessPart = new int [numprocs];
        startPosForSortedEachProcessPart[0] = 0;
        for(int i = 1; i < numprocs; i++) startPosForSortedEachProcessPart[i] = startPosForSortedEachProcessPart[i - 1] + correspondingPartSizeArray[i - 1];
        ret_totalSize = startPosForSortedEachProcessPart[numprocs - 1] + correspondingPartSizeArray[numprocs - 1];
    }
    else
        MPI_Send(&correspondingPartSize, 1, MPI_INT, MASTER_PROC, myid, MPI_COMM_WORLD);

    int *ret_dataSet = nullptr;
    if(MASTER_PROC == myid)
    {
        ret_dataSet = new int [ret_totalSize];
        for(int i = 0; i < numprocs; i++)
        {
            if(i == MASTER_PROC)
            {
                for(int j = 0; j < correspondingPartSizeArray[i]; j++)
                    ret_dataSet[startPosForSortedEachProcessPart[i] + j] = correspondingPartDataSet[j];
            }
            else
                MPI_Recv(&ret_dataSet[startPosForSortedEachProcessPart[i]], correspondingPartSizeArray[i], MPI_INT, i, i, MPI_COMM_WORLD, &s);
        }
    }
    else
        MPI_Send(correspondingPartDataSet, correspondingPartSize, MPI_INT, MASTER_PROC, myid, MPI_COMM_WORLD);

#ifdef TEST_PART
    cout << myid << ": Step 6 end" << endl;
#endif

    //Test part
    if(MASTER_PROC == myid)
    {
        ofstream fout(string(DEFAULT_OUTPUT_FILE));
        if(fout.is_open())
        {
            for (int i = 0; i < ret_totalSize; i++) fout << ret_dataSet[i] << " ";
            fout << endl;
        }
    }

    MPI_Finalize();

    return 0;
}