#ifndef LCP_H
#define LCP_H

#include <iostream>
#include <vector>
#include <string>
#include <cassert>

typedef std::vector<double> Array1D;
typedef std::vector<Array1D> Array2D;

typedef std::vector<unsigned int> Array1UI;
typedef std::vector<Array1UI> Array2UI;

template<typename T> void print(T& v);
void print(Array2D& m);
void print(Array2UI& b);

struct Tableau {
  Array2UI basis;
  Array2D data;
  int q_col, z0_col, z0_row;
  int n_rows, n_cols;

  Array1UI get_elem_from_basis(unsigned int i);
  Array1UI elem_from_col(unsigned int i);
  int col_from_elem(const Array1UI& b);
  Array1UI change_in(unsigned int i, const Array1UI& b);
  void print_elem(const Array1UI& b);
  void print_sol();
  void print_status();
  void print_status(int ctr, const Array1UI& bi, const Array1UI& bo);
  int find_initial_pivot();
  int find_pivot(int col_p);
  void reduce_from_with_pivot(int i, int row_p, int col_p);
  void pivot_reduction(int row_p, int col_p);
};

struct LCP {
  Array1D q;
  Array2D M;

  void read_LCP(void);
  void print_LCP(void);
  void make_tableau(Tableau& tableau);
  void lemke_algorithm(void);
};

template<typename T>void print(T& v) {
  for (int i=0; i<v.size(); i++)
    std::cout << v[i] << " ";
  std::cout << std::endl;
}

void print(Array2D& m) {
  for (int i=0; i<m.size(); i++)
    print(m[i]);
}

void print(Array2UI& b) {
  for (int i=0; i<b[0].size(); i++) {
    std::cout << ((b[0][i]==0)? "w": "z");
    std::cout << b[1][i];
    std::cout << " ";
  }
  std::cout << std::endl;
}

void LCP::read_LCP() {
  int n;
  // problem size
  std::cin >> n;
  q.resize(n,0.0);
  M.resize(n,Array1D(n,0.0));
  // set q vector
  for (int i=0; i<q.size(); i++)
    std::cin >> q[i];
  // set M vector
  for (int i=0; i<M.size(); i++)
    for (int j=0; j<M[i].size(); j++)
      std::cin >> M[i][j];
}

void LCP::print_LCP() {
  // print problem statement: size, q and M
  std::cout << q.size() << std::endl;
  print(q);
  print(M);
}

void LCP::make_tableau(Tableau& tableau) {
  // make lemke's tableau and basis
  int n = q.size();
  tableau.data.resize(n, Array1D(2*n+2,0.0));
  tableau.basis.resize( 2, Array1UI( n, 0.0) );
  for (int i=0; i<n; i++) {
    // data
    tableau.n_rows = n;
    tableau.n_cols = 2*n+2;
    tableau.z0_col = 2*n;
    tableau.q_col = 2*n+1;
    tableau.data[i][i] = 1.0;
    tableau.data[i][2*n] = -1.0;
    tableau.data[i][2*n+1] = q[i];
    // basis
    tableau.basis[0][i]=0.0; tableau.basis[1][i]=i+1.0;
  }
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      tableau.data[i][j+n] = M[i][j];
}

void LCP::lemke_algorithm() {

  Tableau tableau;

  make_tableau(tableau);

  std::cout << "INITIAL TABLEAU: " << std::endl;
  tableau.print_status();

  // step 1.
  // find initial row pivot
  int row_p = tableau.find_initial_pivot();
  int col_p = tableau.z0_col;
  if (row_p==-1) {
    std::cout << "Find Solution!!!\nz = "; print(q);
    return;
  }
  // step 2. row operations in idx-row to enter basis
  tableau.pivot_reduction( row_p, col_p );
  int counter=1, max_iter=6;
  do {
    // step 3. with row_p_old exit find complement enter the basis
    // b_i -> Basis -> b_o
    Array1UI b_i=tableau.elem_from_col( col_p );
    Array1UI b_o = tableau.change_in( row_p, b_i );
    tableau.print_status(counter, b_i, b_o);
    // find complement
    b_o[0] = (b_o[0]==0 && b_o[1]!=0)? 1: 0;
    col_p = tableau.col_from_elem( b_o );
    // exit condition
    if (row_p==tableau.z0_row && 1<counter) {
      std::cout << "Find Solution!!!\nz = "; tableau.print_sol();
      return;
    }
    // step 4.
    row_p = tableau.find_pivot( col_p );
    if (row_p==-1) {
      std::cout << "Ray Solution!!!\n";
      return;
    }
    // step 5.
    tableau.pivot_reduction(row_p,col_p);
  } while(counter++<max_iter);
  std::cout << "Max iteration done!!!\n";
}

Array1UI Tableau::get_elem_from_basis(unsigned int i) {
  assert( i<n_rows );
  Array1UI e(2,0);
  e[0] = basis[0][i];
  assert((e[0]==0)||(e[0]==1));
  e[1] = basis[1][i];
  return e;
}

Array1UI Tableau::elem_from_col(unsigned int i) {
  assert( i<n_cols-1 );
  Array1UI e(2,0);
  e[0] = (i<n_rows)? 0: 1;
  i = ( i<n_rows )? i+1:
    ( i<2*n_rows )? i+1-n_rows: 0;
  assert((e[0]==0)||(e[0]==1));
  assert( i<=n_rows );
  e[1] = i;
  return e;
}

int Tableau::col_from_elem(const Array1UI& b) {
  return ( (b[0]==0)? b[1]-1:
	   (b[0]==1 && 0<b[1])?
	   b[1]+n_rows-1:
	   z0_col );
}

Array1UI Tableau::change_in(unsigned int i, const Array1UI& b) {
  Array1UI b_o = get_elem_from_basis( i );
  basis[0][i] = b[0];
  basis[1][i] = b[1];
  return b_o;
}

void Tableau::print_elem(const Array1UI& b) {
  std::cout << ( (b[0]==0)? "w": "z" ) << b[1];
}

void Tableau::print_status() {
  std::cout << "basis:\n"; print(basis);
  print(data);
}

void Tableau::print_status( int counter,
			    const Array1UI& b_i,
			    const Array1UI& b_o ) {
  std::cout << "TABLEAU " << counter;
  std::cout << ": ("; print_elem(b_i);
  std::cout << "->B->"; print_elem(b_o);
  std::cout << ")" << std::endl;
  std::cout << "basis:\n"; print(basis);
  print(data);
}

void Tableau::print_sol() {
  for (int i=0; i<n_rows; i++)
	std::cout << data[i][q_col] << " ";
  std::cout << std::endl;
}

void Tableau::reduce_from_with_pivot(int row_i, int row_p, int col_p) {
  int n = data.size();
  int m = 2*n+2;
  double pivot = data[row_p][col_p];
  if (row_i==row_p)
    for (int j=0; j<n_cols; j++)
      data[row_p][j]/=pivot;
  else {
    double bij=data[row_i][col_p]/pivot;
    for (int j=0; j<n_cols; j++)
      data[row_i][j]-=data[row_p][j]*bij;
  }
}

void Tableau::pivot_reduction(int row_p, int col_p) {
  for (int i=0; i<n_rows; i++) {
    if (i==row_p)
      continue;
    reduce_from_with_pivot(i,row_p,col_p);
  }
  reduce_from_with_pivot(row_p,row_p,col_p);
}

int Tableau::find_pivot(int col_p) {
  int row_p = -1;
  double rmin = 1e18;
  for (int i=0; i<n_rows; i++)
    if (data[i][q_col]/data[i][col_p]<=rmin &&
	0.0<data[i][col_p]) {
      row_p=i;
      rmin=data[i][q_col]/data[i][col_p];
    }
  return row_p;
}

int Tableau::find_initial_pivot() {
  int row_p = -1;
  double qmin = 1e18;
  for (int i=0; i<n_rows; i++) // exist(qi<0) and argmin_i qi
    if (data[i][q_col]<=qmin
	&& data[i][q_col]<0.0) {
      row_p = i;
      qmin = data[i][q_col];
    }
  return z0_row=row_p;
}

#endif // #ifndef LITTLEMATH_H
