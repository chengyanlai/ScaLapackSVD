#ifndef __INIT_GRID_HPP__
#define __INIT_GRID_HPP__

#include <stdio.h>
#include "CBLASInterface.hpp"

struct GridParameters{
  int Context;
  int Rank;
  int NumProcs;
  int Row;
  int Col;
  int NumRows;
  int NumCols;
};

class InitGrid
{
private:
// protected:
  int Rank;
  int NumProcs;
  int Context;
  int GridNumProcsRows, GridNumProcsCols;
  int RankRow, RankCol;
  int TotalRows, TotalCols;
  int Rows, Cols;
  int BlockSizeRows, BlockSizeCols;
  int ProcsWithFirstRow, ProcsWithFirstCol;
  GridParameters GridInfo;

  inline void SetProcsGrid(){
    Cblacs_pinfo(&Rank, &NumProcs);
    Cblacs_get(0, 0, &Context);
    Cblacs_gridinit(&Context, "Row-major", GridNumProcsRows, GridNumProcsCols);
    Cblacs_pcoord(Context, Rank, &RankRow, &RankCol);
    Cblacs_barrier(Context, "All");
  };
  inline int NumRoC(int NumRowsCols, int RowColBlockSize, int ProcsRowColCoord, int ProcsWithFirstRowCol, int NumRowColProcs){
    int ExtraBlocks, Dist, NBlocks, NumLocalRowsCols;

    Dist = (NumRowColProcs + ProcsRowColCoord - ProcsWithFirstRowCol) % NumRowColProcs;
    NBlocks = NumRowsCols / RowColBlockSize;

    NumLocalRowsCols = (NBlocks/NumRowColProcs) * RowColBlockSize;
    ExtraBlocks = NBlocks % NumRowColProcs;

    if ( Dist < ExtraBlocks ) NumLocalRowsCols += RowColBlockSize;
    else if ( Dist == ExtraBlocks) NumLocalRowsCols += (NumRowsCols % RowColBlockSize);

    return NumLocalRowsCols;
  };
  inline void SetProcsRowsOrCols(){
    Rows = NumRoC(TotalRows, BlockSizeRows, RankRow, ProcsWithFirstRow, GridNumProcsRows);
    Cols = NumRoC(TotalCols, BlockSizeCols, RankCol, ProcsWithFirstCol, GridNumProcsCols);
    GridInfo.Context = Context;
    GridInfo.Rank = Rank;
    GridInfo.NumProcs = NumProcs;
    GridInfo.Row = RankRow;
    GridInfo.Col = RankCol;
    GridInfo.NumRows = Rows;
    GridInfo.NumCols = Cols;
  };
public:
  InitGrid(int totalRows, int totalCols, int blockSizeRows, int blockSizeCols, int gridNumProcRows, int gridNumProcCols, int procWithFirstRow, int procWithFirstCol):
  TotalRows(totalRows), TotalCols(totalCols),
  BlockSizeRows(blockSizeRows), BlockSizeCols(blockSizeCols),
  GridNumProcsRows(gridNumProcRows), GridNumProcsCols(gridNumProcCols),
  ProcsWithFirstRow(procWithFirstRow),ProcsWithFirstCol(procWithFirstCol){
    // initiaize the rectangular grid
    SetProcsGrid();
    // allocate number of rows and cols for each proc
    SetProcsRowsOrCols();
  };
  virtual ~InitGrid(){};

  GridParameters GetGridInfo()const{return GridInfo;};
  inline int GetRank()const{return Rank;};
  inline int GetNumProcs()const{return NumProcs;};
  inline int GetContext()const{return Context;};
  inline int GetRow()const{return RankRow;};
  inline int GetCol()const{return RankCol;};
  inline int GetNumRows()const{return Rows;};
  inline int GetNumCols()const{return Cols;};
};

#endif /* end of include guard: __INIT_GRID_HPP__ */
