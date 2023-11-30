#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <vector>

class FlowGraph {

public:
  explicit FlowGraph(int n) : _capacity(n), _size(0) {
    edges = new int[_capacity * _capacity]();
    dfn = new int[_capacity]();
    rpon = new int[_capacity]();
    visited = new int[_capacity]();
    FATHER = new int[_capacity]();
    ND = new int[_capacity]();
    HIGHPT = new int[_capacity];
    for (int i = 0; i < _capacity; i++) {
      HIGHPT[i] = -1;
    }
  }

  ~FlowGraph() {
    delete[] edges;
    delete[] dfn;
    delete[] rpon;
    delete[] visited;
    delete[] FATHER;
    delete[] ND;
    delete[] HIGHPT;
  }

  void addVertex(std::string str);

  void addEdge(std::string str1, std::string str2);

  void dfs();

  int checkReducibility();

  // dot -T png -o zsy_test1.png zsy_test1.dot; sxiv zsy_test1.png
  void dumpDot(std::string filename);

private:
  int findIdByName(std::string str);
  int dfs(int n);
  int isAncestor(int u, int v);

  int _capacity;
  int _size;
  std::vector<std::string> vertex_names{};
  std::map<std::string, int> vertex_name_id_map{};
  // 0: no edge
  // 1: has edge, but yet handled
  // 2: tree edge
  // 3: fronds
  // 4: reverse fronds
  // 5: cross links
  int *edges;
  int *dfn;
  int *rpon;
  int *visited;
  int dfs_i, dfs_c;
  int *FATHER;
  int *ND;
  int *HIGHPT;
};

void FlowGraph::addVertex(std::string str) {
  assert(_size + 1 <= _capacity);
  vertex_names.push_back(str);
  vertex_name_id_map.emplace(str, _size++);
}

int FlowGraph::findIdByName(std::string str) {
  auto it = vertex_name_id_map.find(str);
  if (it == vertex_name_id_map.end())
    return -1;
  return it->second;
}

void FlowGraph::addEdge(std::string str1, std::string str2) {
  int id1 = findIdByName(str1);
  int id2 = findIdByName(str2);
  assert(id1 != -1 && id2 != -1);
  edges[id1 * _capacity + id2] = 1;
}

void FlowGraph::dumpDot(std::string filename) {
  std::ofstream outfile;
  outfile.open(filename, std::ios::out | std::ios::trunc);

  outfile << "digraph G {\n";
  outfile << "label=\"[Vertex]: [DFN], [RPO], [ND], [HIGHPT]\"";
  // Vertex: dfn, rpon
  for (int i = 0; i < _size; i++) {
    outfile << "\t" << vertex_names[i] << " [label =\"" << vertex_names[i]
            << ": " << dfn[i] << ", " << rpon[i] << ", " << ND[i] << ", "
            << (HIGHPT[i] == -1 ? "NULL" : vertex_names[HIGHPT[i]])
            << "\", shape=circle];\n";
  }
  outfile << "\n";

  for (int i = 0; i < _size; i++) {
    for (int j = 0; j < _size; j++) {
      if (edges[i * _capacity + j] > 0) {
        outfile << "\t" << vertex_names[i] << " -> " << vertex_names[j];
        if (edges[i * _capacity + j] == 2) {
          outfile << " [label=\"T\"];\n";
        } else if (edges[i * _capacity + j] == 3) {
          outfile << " [label=\"F\", style=dashed, color=red];\n";
        } else if (edges[i * _capacity + j] == 4) {
          outfile << " [label=\"R\", style=dashed, color=blue];\n";
        } else if (edges[i * _capacity + j] == 5) {
          outfile << " [label=\"C\", style=dashed, color=darkgreen];\n";
        } else {
          assert(1);
        }
      }
    }
  }

  outfile << "}\n";

  outfile.close();
}

int FlowGraph::dfs(int n) {
  visited[n] = 1;
  dfn[n] = dfs_i++;
  for (int s = 0; s < _size; s++) {
    if (edges[n * _capacity + s] > 0 && visited[s] == 0) {
      assert(edges[n * _capacity + s] == 1);
      edges[n * _capacity + s] = 2;
      FATHER[s] = n;
      ND[n] += (dfs(s) + 1);
    }
  }
  rpon[n] = dfs_c--;
  return ND[n];
}

void FlowGraph::dfs() {
  FATHER[0] = -1;
  dfs_i = 0;
  dfs_c = _capacity - 1;
  dfs(0);
  for (int i = 0; i < _capacity; i++) {
    for (int j = 0; j < _capacity; j++) {
      if (edges[i * _capacity + j] == 1) {
        if ((dfn[i] > dfn[j]) && (rpon[i] > rpon[j])) {
          // fronds
          edges[i * _capacity + j] = 3;
        } else if ((dfn[i] < dfn[j]) && (rpon[i] < rpon[j])) {
          // reverse fronds
          edges[i * _capacity + j] = 4;
        } else if (((dfn[i] < dfn[j]) && (rpon[i] > rpon[j])) ||
                   ((dfn[i] > dfn[j]) && (rpon[i] < rpon[j]))) {
          // cross links
          edges[i * _capacity + j] = 5;
        } else {
          assert(1);
        }
      }
    }
  }
}

int FlowGraph::checkReducibility() {
  for (int w = _capacity - 1; w >= 0; w--) {
    for (int u = 0; u < _capacity; u++) {
      if (edges[u * _capacity + w] == 3) {
        std::cout << "Check frond(" << vertex_names[u] << ", "
                  << vertex_names[w] << ")\n";
        std::queue<int> CHECK{};
        CHECK.push(u);
        while (!CHECK.empty()) {
          int tmp_u = CHECK.front();
          CHECK.pop();
          if (!isAncestor(w, tmp_u)) {
            std::cout << "[ERROR] " << vertex_names[w]
                      << " is not the ancestor of " << vertex_names[tmp_u]
                      << " on DFST, NOT Reduciable! STOP!\n";
            return 1;
          }
          std::cout << "\t" << vertex_names[w] << " is ancestor of "
                    << vertex_names[tmp_u] << " on DFST\n";
          while (tmp_u != w) {
            if (HIGHPT[tmp_u] == -1) {
              HIGHPT[tmp_u] = w;
              std::cout << "\tSet HIGHPT[" << vertex_names[tmp_u]
                        << "] := " << vertex_names[w] << "\n";
              for (int tmp_v = 0; tmp_v < _capacity; tmp_v++) {
                if (edges[tmp_v * _capacity + tmp_u] == 5) {
                  std::cout << "\tFound Cross-link(" << vertex_names[tmp_v]
                            << ", " << vertex_names[tmp_u] << ") and add "
                            << vertex_names[tmp_v] << " to CHECK\n";
                  CHECK.push(tmp_v);
                }
              }
            } else {
              std::cout << "\tHIGHPT[" << vertex_names[tmp_u]
                        << "] has been set, SKIP!\n";
            }
            tmp_u = FATHER[tmp_u];
          }
        }
      }
    }
  }
  for (int u = 0; u < _capacity; u++) {
    for (int v = 0; v < _capacity; v++) {
      if (edges[u * _capacity + v] == 4) {
        std::cout << "Check reverse-frond(" << vertex_names[u] << ", "
                  << vertex_names[v] << ")\n";
        if (dfn[u] < dfn[HIGHPT[v]]) {
          std::cout << "[ERROR] dfn[" << vertex_names[u] << "]: " << dfn[u]
                    << " < dfn[HIGHPT[" << vertex_names[v]
                    << "]]: " << dfn[HIGHPT[v]] << ", NOT Reduciable! STOP!\n";
          return 2;
        } else {
          std::cout << "\tdfn[" << vertex_names[u] << "]: " << dfn[u]
                    << " >= dfn[HIGHPT[" << vertex_names[v]
                    << "]]: " << dfn[HIGHPT[v]] << "\n";
        }
      }
    }
  }
  std::cout << "Test Successful! This flowgraph is Reduciable!\n";
  return 0;
}

int FlowGraph::isAncestor(int u, int v) {
  return (dfn[u] < dfn[v]) && (dfn[v] <= dfn[u] + ND[u]);
}

void zsy_test1() {
  std::cout << "\n[Test 1]:\n";
  FlowGraph g(12);
  g.addVertex("S");
  g.addVertex("C");
  g.addVertex("B");
  g.addVertex("E");
  g.addVertex("A");
  g.addVertex("D");
  g.addVertex("H");
  g.addVertex("K");
  g.addVertex("F");
  g.addVertex("G");
  g.addVertex("J");
  g.addVertex("I");

  g.addEdge("S", "B");
  g.addEdge("S", "C");
  g.addEdge("B", "A");
  g.addEdge("B", "D");
  g.addEdge("B", "E");
  g.addEdge("A", "D");
  g.addEdge("E", "H");
  g.addEdge("D", "H");
  g.addEdge("H", "K");
  g.addEdge("C", "F");
  g.addEdge("C", "G");
  g.addEdge("F", "I");
  g.addEdge("G", "I");
  g.addEdge("G", "J");
  g.addEdge("J", "I");
  g.addEdge("I", "K");
  g.addEdge("I", "S");
  g.addEdge("K", "S");
  // which matters
  g.addEdge("H", "B");

  g.dfs();
  g.checkReducibility();
  g.dumpDot("zsy_test1.dot");
}

void zsy_test2() {
  std::cout << "\n[Test 2]:\n";
  FlowGraph g(12);
  g.addVertex("S");
  g.addVertex("C");
  g.addVertex("B");
  g.addVertex("E");
  g.addVertex("A");
  g.addVertex("D");
  g.addVertex("H");
  g.addVertex("K");
  g.addVertex("F");
  g.addVertex("G");
  g.addVertex("J");
  g.addVertex("I");

  g.addEdge("S", "B");
  g.addEdge("S", "C");
  g.addEdge("B", "A");
  g.addEdge("B", "D");
  g.addEdge("B", "E");
  g.addEdge("A", "D");
  g.addEdge("E", "H");
  g.addEdge("D", "H");
  g.addEdge("H", "K");
  g.addEdge("C", "F");
  g.addEdge("C", "G");
  g.addEdge("F", "I");
  g.addEdge("G", "I");
  g.addEdge("G", "J");
  g.addEdge("J", "I");
  g.addEdge("I", "K");
  g.addEdge("I", "S");
  g.addEdge("K", "S");
  // which matters
  g.addEdge("H", "E");

  g.dfs();
  g.checkReducibility();
  g.dumpDot("zsy_test2.dot");
}

void zsy_test3() {
  std::cout << "\n[Test 3]:\n";
  FlowGraph g(8);
  g.addVertex("B1");
  g.addVertex("B2");
  g.addVertex("B3");
  g.addVertex("B4");
  g.addVertex("B5");
  g.addVertex("B6");
  g.addVertex("B7");
  g.addVertex("B8");

  g.addEdge("B1", "B2");
  g.addEdge("B2", "B3");
  g.addEdge("B3", "B4");
  g.addEdge("B3", "B5");
  g.addEdge("B4", "B6");
  g.addEdge("B5", "B6");
  g.addEdge("B5", "B7");
  g.addEdge("B6", "B5");
  g.addEdge("B6", "B8");
  g.addEdge("B8", "B1");

  g.dfs();
  g.checkReducibility();
  g.dumpDot("zsy_test3.dot");
}

int main() {
  std::cout << "Test Reduciability by zhaosiying12138 from LiuYueCity Academy "
               "of Sciences\n";
  zsy_test1();
  zsy_test2();
  zsy_test3();

  return 0;
}
