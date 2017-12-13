#include <iostream>
#include <vector>

struct LCP {
  typedef std::vector<double> Array1D;
  typedef std::vector< Array1D > Array2D;
  Array1D q;
  Array2D M;

  void read(void);
  void print(void);
};

int main() {

  // lemke's algorithm
  LCP problem;
  
  problem.read();
  
  problem.print();

  // make lemke's tableau
  int n = problem.q.size();
  LCP::Array2D tableau(n, LCP::Array1D(2*n+2,0.0));
  for (int i=0; i<n; i++) {
    tableau[i][i] = 1.0;
    tableau[i][2*n] = -1.0;
    tableau[i][2*n+1] = problem.q[i];
  }
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      tableau[i][j+n] = problem.M[i][j];

  return 0;
}

void LCP::read() {
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

void LCP::print() {
  // print problem statement: size, q and M
  std::cout << q.size() << std::endl;
  for (int i=0; i<q.size(); i++)
    std::cout << q[i] << " ";
  std::cout << std::endl;
  for (int i=0; i<M.size(); i++) {
    for (int j=0; j<M[i].size(); j++)
      std::cout << M[i][j] << " ";
    std::cout << std::endl;
  }
}
