#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

class FlowGraph {

public:
  explicit FlowGraph(int n) : _capacity(n), _size(0) {
    edges = new int[_capacity * _capacity]();
    dfn = new int[_capacity]();
    rpon = new int[_capacity]();
    visited = new int[_capacity]();
  }

  ~FlowGraph() { delete[] edges; }

  void addVertex(std::string str);

  void addEdge(std::string str1, std::string str2);

  void dfs();

  // dot -T png -o zsy_test1.png zsy_test1.dot; sxiv zsy_test1.png
  void dumpDot(std::string filename);

private:
  int findIdByName(std::string str);
  void dfs(int n);

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
  // Vertex: dfn, rpon
  for (int i = 0; i < _size; i++) {
    outfile << "\t" << vertex_names[i] << " [label =\"" << vertex_names[i]
            << ": " << dfn[i] << ", " << rpon[i] << "\", shape=circle];\n";
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

void FlowGraph::dfs(int n) {
  visited[n] = 1;
  dfn[n] = dfs_i++;
  for (int s = 0; s < _size; s++) {
    if (edges[n * _capacity + s] > 0 && visited[s] == 0) {
      assert(edges[n * _capacity + s] == 1);
      edges[n * _capacity + s] = 2;
      dfs(s);
    }
  }
  rpon[n] = dfs_c--;
}

void FlowGraph::dfs() {
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

int main() {
  std::cout << "Hello Reduciability\n";
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
  g.dumpDot("zsy_test1.dot");
}
