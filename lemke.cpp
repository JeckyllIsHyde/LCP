#include <iostream>
#include <vector>
#include <string>

typedef std::vector<double> Array1D;
typedef std::vector<Array1D> Array2D;

typedef std::vector<unsigned int> Array1UI;
typedef std::vector<Array1UI> Array2UI;

template<typename T> void print(T& v);
void print(Array2D& m);
void print(Array2UI& b);

struct LCP {
  Array1D q;
  Array2D M;
  Array2D tableau;

  void read_LCP(void);
  void print_LCP(void);
  void make_tableau();
  void lemke_algorithm();

  void reduce_from_with_pivot(int i, int row_p, int col_p);
  void pivot_reduction(int& row_p, int& col_p);
  void find_initial_pivot(int& row_p);
  void find_pivot(int& row_p, int& col_p);
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

  int n = tableau.size();
  // create and init basis indicator
  Array2UI basis( 2, Array1UI( n, 0.0) );
  for (int i=0; i<n; i++) {
    basis[0][i]=0.0; basis[1][i]=i+1.0;
  }
  std::cout << "basis:\n";  print(basis);
  print(tableau);
  int col_z0=2*n, row_z0, row_p;

  // step 1.
  // find initial row pivot
  find_initial_pivot(row_p);
  // std::cout << qmin << std::endl;
  if (row_p==-1) {
    std::cout << "Find Solution!!!\nz = ";
    print(q);
    return;
  }
  // z0 -> Basis and basis -> wi
  int row_p_old;
  basis[0][row_p] = 1;
  basis[1][row_p] = row_z0 = row_p;
  row_p_old = row_p;
  int col_p = col_z0;

  // step 2. row operations in idx-row to enter basis
  pivot_reduction( row_p, col_z0 );
  row_z0 = row_p_old = row_p;

  std::cout << "TABLEAU " << 1 << std::endl;
  std::cout << "basis:\n";  print(basis);
  print(tableau);

  int counter=2, max_iter=3;
  do {
    // step 3. with row_p_old exit find complement enter the basis

    // step 4.
    //    row_p_old = row_p;
    find_pivot( row_p, col_p );
    if (row_p==-1) {
      std::cout << "Ray Solution!!!\n";
      return;
    }

    // step 5.
    pivot_reduction(row_p,col_p);
    std::cout << "TABLEAU " << counter << std::endl;
    std::cout << "basis:\n"; print(basis);
    print(tableau);
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

void LCP::find_pivot(int& row_p, int& col_p) {
  int n = tableau.size();
  row_p = -1;
  double rmin = 1e18;
  for (int i=0; i<n; i++)
    if (tableau[i][2*n+1]/tableau[i][col_p]<=rmin &&
	0.0<tableau[i][col_p]) {
      row_p=i;
      rmin=tableau[i][2*n+1]/tableau[i][col_p];
    }
}

void LCP::find_initial_pivot(int& row_p) {
  row_p = -1;
  double qmin = 1e18;
  for (int i=0; i<q.size(); i++) // exist(qi<0) and argmin_i qi
    if (q[i]<=qmin
	&& q[i]<0.0) {
      row_p = i;
      qmin = q[i];
    }
}
