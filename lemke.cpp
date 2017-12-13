#include <iostream>
#include <vector>

typedef std::vector<double> Array1D;
typedef std::vector< Array1D > Array2D;

void print(Array1D& v);
void print(Array2D& m);

struct LCP {
  Array1D q;
  Array2D M;

  void read_LCP(void);
  void print_LCP(void);
};

int main() {

  // lemke's algorithm
  LCP problem;
  
  problem.read();
  
  problem.print();

  // make lemke's tableau
  int n = problem.q.size();
  Array2D tableau(n, Array1D(2*n+2,0.0));
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
