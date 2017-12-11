#ifndef __CBLASInterface_HPP__
#define __CBLASInterface_HPP__

// Cblacs declarations not declared in MKL
extern "C" {
  void Cblacs_pinfo(int *rank, int *nprocs);
  void Cblacs_get(int context, int what, int *val);
  void Cblacs_gridinit(int *context, const char *layout, int proc_rows, int proc_cols);
  void Cblacs_pcoord(int context, int rank, int *row, int *col);
  void Cblacs_gridexit(int);
  void Cblacs_barrier(int context, const char *scope);
  void Cdgerv2d(int, int, int, double*, int, int, int);
  void Cdgesd2d(int, int, int, double*, int, int, int);
}

#endif /* end of include guard: __CBLASInterface_HPP__ */
