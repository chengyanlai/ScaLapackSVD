#ifndef __SCALAPACKSVD_HPP__
#define __SCALAPACKSVD_HPP__

#include <vector>
#include <algorithm>
#include <mkl_scalapack.h>

#include "CBLASInterface.hpp"

#define DTYPE_ 0
#define CTXT_  1
#define M_     2
#define N_     3
#define MB_    4
#define NB_    5
#define RSRC_  6
#define CSRC_  7
#define LLD_   8

class ScaLapackSVD{
private:
  int Rank;
  int NumProcs;
  int Context;
  int GridNumProcsRows, GridNumProcsCols;
  int RankRow, RankCol;
  int TotalRows, TotalCols;
  int Rows, Cols;
  int BlockSizeRows, BlockSizeCols;
  int ProcsWithFirstRow, ProcsWithFirstCol;
  std::vector<double> *Data;

  // ScaLAPACK variables
  int size, sizep, sizeq;
  char jobu, jobvt;
  MKL_INT info, ia, ja, iu, ju, ivt, jvt, lwork;
  std::vector<MKL_INT> descA;
  std::vector<MKL_INT> descU;
  std::vector<MKL_INT> descVT;

  // result variables
  std::vector<double> S;
  std::vector<double> VT;
  std::vector<double> U;
  std::vector<double> work;

  void CreateArrayDescriptor(std::vector<MKL_INT> &descVec, int dtype, int ctxt, int m, int n, int mb, int nb, int rsrc, int csrc, int lld);
  void InitSVDVariables();
public:
  ScaLapackSVD(GridParameters GridInfo, int totalRows, int totalCols, int blockSizeRows, int blockSizeCols, int gridNumProcRows, int gridNumProcCols, int procWithFirstRow, int procWithFirstCol);
  virtual ~ScaLapackSVD(){};

  void BuildLocalBlockMatrix(const std::vector<double> &CoordData);
  void Compute();
  inline std::vector<double> GetS()const{return S;};
  inline std::vector<double> GetU()const{return U;};
  inline std::vector<double> GetVT()const{return VT;};
  inline void PrintLocalMatrix()const{
    printf("Rank: %d, local A:\n", Rank);
    for (int j = 0; j < Rows; j++){
      for (int i = 0; i < Cols; i++){
        printf("  %d", (int)Data->operator[](j + i * Rows));
      }
      printf("\n");
    }
    printf("\n");
  };
  inline void PrintS()const{
    printf("Rank: %d, lwork: %d, singular values: \n", Rank, (int)work[0]);
    for (int i = 0; i < size; i++){
      printf("%f\n", S[i]);
    }
    printf("\n");
  };
  inline void PrintU()const{
    printf("Rank: %d, U: \n", Rank);
    for (int i = 0; i < Rows ; i++){
      for (int j = 0; j < Cols; j++){
        printf("  %f", U[i + j*Rows]);
      }
      printf("\n");
    }
  };
  inline void PrintVT()const{
    printf("Rank: %d, VT: \n", Rank);
    for (int i = 0; i < Rows ; i++){
      for (int j = 0; j < Cols; j++){
        printf("  %f", VT[i + j*Rows]);
      }
      printf("\n");
    }
  };
};
#endif/* end of include guard: __SCALAPACKSVD_HPP__ */
