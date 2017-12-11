#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <random>       // std::rand
#include <algorithm>    // std::random_shuffle

#include <mpi.h>

#include "InitGrid.hpp"
#include "ScaLapackSVD.hpp"


int main(int argc, char *argv[]){
  const int root = 0;
  clock_t t0, t1, t2, time;
  t0 = clock();

  MPI_Init(NULL, NULL);

  int Rows = 4000;//atoi(argv[1]);
  int Cols = 4000;//atoi(argv[2]);
  int BlockSizeRows = 2000;//atoi(argv[3]);
  int BlockSizeCols = 2000;//atoi(argv[4]);

  // Randomize a matrix

  t1 = clock();
  int NumProcs;
  int GridProcsRows;
  int GridProcsCols = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
  GridProcsRows = NumProcs;
  InitGrid Grid(Rows, Cols, BlockSizeRows, BlockSizeCols, GridProcsRows, GridProcsCols, root, root);
  ScaLapackSVD SLSVD(Grid.GetGridInfo(), Rows, Cols, BlockSizeRows, BlockSizeCols, GridProcsRows, GridProcsCols, root, root);
  t2 = clock();

  time = t2 - t1;
  MPI_Reduce(&time, &t2, 1, MPI_LONG, MPI_SUM, root, MPI_COMM_WORLD);
  if ( Grid.GetRank() == root )
    printf("rank: %d, init time: %f\n", Grid.GetRank(), (float)(t2)/(CLOCKS_PER_SEC*Grid.GetNumProcs()));

  t1 = clock();
  // Build Block cyclic data for each proc
  {
    // Need to find a way to load a matrix to test
    // CReadData readCSV(fileName, *delimiter.c_str());
    // readCSV.readAllLines();
    std::vector<double> mat(Rows*Cols, 0.0e0);
    std::random_shuffle( mat.begin(), mat.end() );

    SLSVD.BuildLocalBlockMatrix(mat);
    SLSVD.PrintLocalMatrix();
    MPI_Barrier(MPI_COMM_WORLD);
  }
  t2 = clock();
  time = t2 - t1;
  MPI_Reduce(&time, &t2, 1, MPI_LONG, MPI_SUM, root, MPI_COMM_WORLD);
  if ( Grid.GetRank() == root )
    printf("rank: %d, data dist. time: %f\n", Grid.GetRank(), (float)(t2)/(CLOCKS_PER_SEC*Grid.GetNumProcs()));

  // get the SVD
  t1 = clock();
  SLSVD.Compute();
  t2 = clock();


  time = t2 - t1;
  MPI_Reduce(&time, &t2, 1, MPI_LONG, MPI_SUM, root, MPI_COMM_WORLD);
  if ( Grid.GetRank() == root)
    printf("rank: %d, SVD time: %f\n", Grid.GetRank(), (float)(t2)/(CLOCKS_PER_SEC*Grid.GetNumProcs()));

  std::vector<double> S = SLSVD.GetS();
//    const std::vector<double> &leftSingularVectors = SLSVD.getLeftSingularVectors();
//    const std::vector<double> &rightSingularVectors = SLSVD.getRightSingularVectors();

  // root will print the singular values
  if ( Grid.GetRank() == root){
      SLSVD.PrintS();
  }
  MPI_Barrier(MPI_COMM_WORLD);


  t2 = clock();
  time = t2 - t0;
  MPI_Reduce(&time, &t2, 1, MPI_LONG, MPI_SUM, root, MPI_COMM_WORLD);
  if ( Grid.GetRank() == root )
    printf("rank: %d, total time: %f\n", Grid.GetRank(), (float)(t2)/(CLOCKS_PER_SEC*Grid.GetNumProcs()));

  MPI_Finalize();
  return 0;
}
