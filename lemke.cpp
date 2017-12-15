#include <iostream>
#include <vector>

typedef std::vector<double> Array1D;
typedef std::vector<Array1D> Array2D;

template<typename T> void print(T& v);
void print(Array2D& m);

struct LCP {
  Array1D q;
  Array2D M;
  Array2D tableau;

  void read_LCP(void);
  void print_LCP(void);
  void make_tableau();
  void lemke_algorithm();

  void reduce_from_with_pivot(int i, int idx, int p_idx);
  void pivot_reduction(int& idx, int& p_idx);
};

int main() {

  // lemke's algorithm
  LCP problem;
  
  problem.read_LCP();
  
  problem.print_LCP();

  problem.lemke_algorithm();

  return 0;
}

template<typename T>void print(T& v) {
  for (int i=0; i<v.size(); i++)
    std::cout << v[i] << " ";
  std::cout << std::endl;
}

void print(Array2D& m) {
  for (int i=0; i<m.size(); i++)
    print(m[i]);
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

void LCP::make_tableau() {
  // make lemke's tableau
  int n = q.size();
  tableau.resize(n, Array1D(2*n+2,0.0));
  for (int i=0; i<n; i++) {
    tableau[i][i] = 1.0;
    tableau[i][2*n] = -1.0;
    tableau[i][2*n+1] = q[i];
  }
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      tableau[i][j+n] = M[i][j];
}

void LCP::lemke_algorithm() {
  make_tableau();

  std::cout << "INITIAL TABLEAU: " << std::endl;
  print(tableau);

  int n = tableau.size();
  // create and init basis indicator
  std::vector<unsigned int> basis(2*n+1,0.0);
  for (int i=0; i<n; i++)
    basis[i]=i+1.0;
  std::cout << "basis:\n";  print(basis);

  // step 1. z0 -> Basis:
  int row_p = -1;
  double qmin = 1e18;
  // find initial row pivot
  for (int i=0; i<n; i++) // exist(qi<0) and argmin_i qi
    if (q[i]<=qmin
	&& q[i]<0.0) {
      row_p = i;
      qmin = q[i];
    }
  // std::cout << qmin << std::endl;
  if (row_p==-1) {
    std::cout << "Find Solution!!!\nz = ";
    print(q);
    return;
  }
  // z0 -> Basis and basis -> wi
  int idx, p_idx, col_p = 2*n;
  // step 2. row operations in idx-row to enter basis
  pivot_reduction( row_p, col_p );
  basis[col_p] = basis[row_p]; basis[row_p] = 0;
  std::cout << "COUNTER: " << 1 << std::endl;
  print(tableau);
  std::cout << "basis:\n";  print(basis);
  int counter=2, max_iter=6;
  do {
    // step 3. i-row exit and j-column complement enter the basis
    std::cout << "  ***************************" << std::endl;
    std::cout << row_p << " exits the basis" << std::endl;
    if (basis[col_p]==0) {
      col_p=col_p-4;
      std::cout << "por aca" << std::endl;
    }
    else
      col_p=row_p+n;
    std::cout << col_p << " enters the basis" << std::endl;
    // step 4.
    row_p = -1;
    double rmin = 1e18;
    for (int i=0; i<n; i++)
      if (tableau[i][2*n+1]/tableau[i][col_p]<=rmin &&
	  0.0<tableau[i][col_p]) {
	row_p=i;
	rmin=tableau[i][2*n+1]/tableau[i][col_p];
      }
    // std::cout << rmin << std::endl;
    std::cout << "out: " << row_p << std::endl;
    if (row_p==-1) {
      std::cout << "Ray Solution!!!\n";
      return;
    }
    // step 5.
    std::cout << col_p << " " << row_p << std::endl;
    pivot_reduction(row_p,col_p);
    basis[col_p] = basis[row_p]; basis[row_p] = 0;
    std::cout << "COUNTER: " << counter << std::endl;
    print(tableau);
    std::cout << "basis:\n";  print(basis);
  } while(counter++<max_iter);
}

void LCP::reduce_from_with_pivot(int row_i, int row_p, int col_p) {
  int n = tableau.size();
  int m = 2*n+2;
  double pivot = tableau[row_p][col_p];
  if (row_i==row_p)
    for (int j=0; j<m; j++)
      tableau[row_p][j]/=pivot;
  else {
    double bij=tableau[row_i][col_p]/pivot;
    for (int j=0; j<m; j++)
      tableau[row_i][j]-=tableau[row_p][j]*bij;
  }
}

void LCP::pivot_reduction(int& row_p, int& col_p) {
  int n = tableau.size();
  for (int i=0; i<n; i++) {
    if (i==row_p)
      continue;
    reduce_from_with_pivot(i,row_p,col_p);
  }
  reduce_from_with_pivot(row_p,row_p,col_p);
}
