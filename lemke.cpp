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
