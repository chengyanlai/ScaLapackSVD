#include <stdio.h>
#include "InitGrid.hpp"
#include "ScaLapackSVD.hpp"

ScaLapackSVD::ScaLapackSVD(GridParameters GridInfo, int totalRows, int totalCols, int blockSizeRows, int blockSizeCols, int gridNumProcRows, int gridNumProcCols, int procWithFirstRow, int procWithFirstCol):
Rank(GridInfo.Rank), NumProcs(GridInfo.NumProcs), Context(GridInfo.Context),
RankRow(GridInfo.Row), RankCol(GridInfo.Col),
Rows(GridInfo.NumRows), Cols(GridInfo.NumCols),
TotalRows(totalRows), TotalCols(totalCols),
BlockSizeRows(blockSizeRows), BlockSizeCols(blockSizeCols),
GridNumProcsRows(gridNumProcRows), GridNumProcsCols(gridNumProcCols),
ProcsWithFirstRow(procWithFirstRow), ProcsWithFirstCol(procWithFirstCol),
descA(9), descU(9), descVT(9),
S(std::min(totalRows, totalCols)),U(Rows * Cols), VT(Rows * Cols),
work(Rows * blockSizeRows){
  // initialize the variables for the SVD routine
  ia = 1;
  ja = 1;
  iu = 1;
  ju = 1;
  ivt = 1;
  jvt = 1;
  lwork = -1; // will first determine the size
  size = std::min(totalRows, totalCols);
  sizeq = Rows;
  sizep = Cols;
  jobu = 'V';
  jobvt = 'V';

  // initialize the local matrices
  Data = new std::vector<double>(Rows * Cols, 0);
}

void ScaLapackSVD::BuildLocalBlockMatrix(const std::vector<double> &CoordData)
{
  int procRow, procCol, procRank;
  int blockRow, blockCol, linearDisp;
  int localRow, localCol;

  int row, col;
  double value;
  for (size_t i = 0; i < CoordData.size(); i+=3){
    row = (int)CoordData[i + 0] - 1;
    col = (int)CoordData[i + 1] - 1;
    value = CoordData[i + 2];

    procRow = (row/BlockSizeRows) % GridNumProcsRows;
    procCol = (col/BlockSizeCols) % GridNumProcsCols;
    procRank = procCol + procRow * GridNumProcsCols;

    if ( procRank == Rank ){
      // block coordinate and the coordinates of "value" in the block
      blockRow = row/(GridNumProcsRows * BlockSizeRows);
      blockCol = col/(GridNumProcsCols * BlockSizeCols);
      localRow = row % (BlockSizeRows + 0);
      localCol = col % (BlockSizeCols + 0);
      linearDisp = localRow + localCol * Rows + blockCol * BlockSizeCols * Rows + blockRow * BlockSizeRows;
      Data->operator[](linearDisp) = value;
    }
  }
}

void ScaLapackSVD::CreateArrayDescriptor(std::vector<MKL_INT> &descVec, int dtype, int ctxt, int m, int n, int mb, int nb, int rsrc, int csrc, int lld){
  // array descriptors
  descVec[DTYPE_] = dtype;
  descVec[CTXT_] = ctxt;
  descVec[M_] = m;
  descVec[N_] = n;
  descVec[MB_] = mb;
  descVec[NB_] = nb;
  descVec[RSRC_] = rsrc;
  descVec[CSRC_] = csrc;
  descVec[LLD_] = lld;
}

void ScaLapackSVD::InitSVDVariables(){
  CreateArrayDescriptor(descA, 1, Context, TotalRows, TotalCols, BlockSizeRows, BlockSizeCols, 0, 0, Rows);
  CreateArrayDescriptor(descU, 1, Context, TotalRows, TotalRows, BlockSizeRows, BlockSizeCols, 0, 0, Rows);
  CreateArrayDescriptor(descVT, 1, Context, TotalCols, TotalCols, BlockSizeRows, BlockSizeCols, 0, 0, Rows);
}

void ScaLapackSVD::Compute(){
  InitSVDVariables();
  // To determine lwork
  pdgesvd(&jobu, &jobvt, &TotalRows, &TotalCols, Data->data(), &ia, &ja, descA.data(), S.data(),
          U.data(), &iu, &ju, descU.data(),
          VT.data(), &ivt, &jvt, descVT.data(),
          work.data(), &lwork, &info);
  if (info != 0)
    printf("rank: %d, info: %d\n", Rank, info);

  if (Rank == 0)
    printf("rank: %d, about to start SVD comp. ...\n", Rank);

  // re-allocate work using returned lwork and run SVD again
  lwork = work[0];
  work.resize(lwork);
  pdgesvd(&jobu, &jobvt, &TotalRows, &TotalCols, Data->data(), &ia, &ja, descA.data(), S.data(),
          U.data(), &iu, &ju, descU.data(),
          VT.data(), &ivt, &jvt, descVT.data(),
          work.data(), &lwork, &info);
  if (info != 0) printf("rank: %d, info: %d\n", Rank, info);
}



