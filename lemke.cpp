#include <iostream>
#include <vector>

typedef std::vector<double> Array1D;
typedef std::vector<Array1D> Array2D;

void print(Array1D& v);
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

void print(Array1D& v) {
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
    basis[i]=1.0;

  // step 1. z0 -> Basis:
  int idx = -1;
  double qmin = 1e18;
  //    exist(qi<0) and argmin_i qi
  for (int i=0; i<n; i++)
    if (tableau[i][2*n+1]<qmin && tableau[i][2*n+1]<0.0) {
      idx = i;
      qmin = tableau[i][2*n+1];
    }
  std::cout << qmin << std::endl;
  if (idx==-1) {
    std::cout << "Find Solution!!!\nz = \n";
    for (int i=0; i<n; i++)
      std::cout << tableau[i][2*n+1] << " ";
    std::cout << std::endl;
    return;
  }
  // z0 -> Basis and basis -> wi
  int p_idx = 2*n;
  // step 2. row operations in idx-row to enter basis
  pivot_reduction(idx,p_idx);
  basis[idx] = 0; basis[p_idx] = 1;
  std::cout << "COUNTER: " << 1 << std::endl;
  print(tableau);
  int counter=2, max_iter=6;
  do {
    // step 3. i-row exit and j-column complement enter the basis
    (idx<n)? idx+=n: idx-=n;
    std::cout << idx << " exits the basis" << std::endl;
    std::cout << idx << " enters the basis" << std::endl;
    // step 4.
    p_idx = idx;
    idx = -1;
    double rmin = 1e18;
    for (int i=0; i<n; i++)
      if (tableau[i][2*n+1]/tableau[i][p_idx]<=rmin &&
	  0.0<tableau[i][p_idx]) {
	idx=i;
	rmin=tableau[i][2*n+1]/tableau[i][p_idx];
      }
    std::cout << rmin << std::endl;
    std::cout << "out: " << idx << std::endl;
    if (idx==-1) {
      std::cout << "Ray Solution!!!\n";
      return;
    }
    // step 5.
    std::cout << p_idx << " " << idx << std::endl;
    pivot_reduction(idx,p_idx);
    basis[idx] = 0; basis[p_idx] = 1;
    std::cout << "COUNTER: " << counter << std::endl;
    print(tableau);
  } while(counter++<max_iter);
}

void LCP::reduce_from_with_pivot(int i, int idx, int p_idx) {
  int n = q.size();
  int m = 2*n+2;
  double pivot = tableau[idx][p_idx];
  if (i==idx)
    for (int j=0; j<m; j++)
      tableau[idx][j]/=pivot;
  else {
    double bij=tableau[i][p_idx]/pivot;
    for (int j=0; j<m; j++)
      tableau[i][j]-=tableau[idx][j]*bij;
  }
}

void LCP::pivot_reduction(int& idx, int& p_idx) {
  int n=q.size();
  for (int i=0; i<n; i++) {
    if (i==idx)
      continue;
    reduce_from_with_pivot(i,idx,p_idx);
  }
  reduce_from_with_pivot(idx,idx,p_idx);
}
