#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>

class FlowGraph {

public:
  class AdjacentView {
  public:
    explicit AdjacentView(const FlowGraph &graph)
        : _graph(graph), _size(graph._size),
          _coalesced_vertex_names(graph.vertex_names) {
      coalesced = new int[_size]();
      for (int i = 0; i < _size; i++) {
        std::set<int> tmp_set{};
        precessors.push_back(tmp_set);
      }
      for (int i = 0; i < _size; i++) {
        std::set<int> tmp_set{};
        successors.push_back(tmp_set);
      }
      for (int u = 0; u < _size; u++) {
        for (int v = 0; v < _size; v++) {
          if (graph.edges[u * _size + v] > 0) {
            precessors[v].insert(u);
            successors[u].insert(v);
          }
        }
      }
    }

    ~AdjacentView() { delete[] coalesced; }

    void T1(int u) {
      std::cout << "T1(" << _coalesced_vertex_names[u] << ")\n";
      assert(precessors[u].find(u) != precessors[u].end() &&
             successors[u].find(u) != successors[u].end());
      precessors[u].erase(u);
      successors[u].erase(u);
    }

    int checkSelfCycle(int u) {
      return precessors[u].find(u) != precessors[u].end();
    }

    void T2(int v) {
      // Do not coalesce entry node
      assert(v != 0);
      // only edge entering v from u
      assert(precessors[v].size() == 1);
      int u = *(precessors[v].begin());
      // u != v, otherwise call T1() transformation
      assert(u != v);
      // this also means v don't have a v->v self-cycle
      assert(successors[v].find(v) == successors[v].end());
      // trival test
      assert(successors[u].find(v) != successors[u].end());
      std::cout << "T2(" << _coalesced_vertex_names[u] << ", "
                << _coalesced_vertex_names[v] << ")\n";

      // coalesce v to u
      coalesced[v] = 1;
      successors[u].erase(v);
      for (int v_succ : successors[v]) {
        precessors[v_succ].erase(v);
        precessors[v_succ].insert(u);
      }
      successors[u].insert(successors[v].begin(), successors[v].end());
      _coalesced_vertex_names[u] += "|" + _coalesced_vertex_names[v];
    }

    void dump() {
      std::cout << "### Dump Collapse Result ###\n";
      std::cout << "vextex_names: ";
      for (int i = 0; i < _size; i++) {
        if (coalesced[i] == 1)
          continue;
        std::cout << _coalesced_vertex_names[i] << ", ";
      }
      std::cout << "\n";

      for (int u = 0; u < _size; u++) {
        if (coalesced[u] == 1)
          continue;
        std::cout << "precessor[" << _coalesced_vertex_names[u] << "]: ";
        for (auto v : precessors[u]) {
          std::cout << _coalesced_vertex_names[v] << ", ";
        }
        std::cout << "\n";

        std::cout << "successors[" << _coalesced_vertex_names[u] << "]: ";
        for (auto v : successors[u]) {
          std::cout << _coalesced_vertex_names[v] << ", ";
        }
        std::cout << "\n";
      }
    }

    int getUncoalescedCount() {
      int cnt = 0;
      for (int i = 0; i < _size; i++) {
        if (coalesced[i] == 0) {
          cnt++;
        }
      }
      return cnt;
    }

    int checkCollapsibilitySlowly() {
      std::cout << "### Check Collapsibility Slowly in O(V^2) ###\n";
      int changed = 1;
      for (int round = 1; changed == 1; round++) {
        std::cout << "Round " << round << ":\n";
        changed = 0;
        for (int u = 0; u < _size; u++) {
          if (coalesced[u] == 1)
            continue;
          if (successors[u].find(u) != successors[u].end()) {
            T1(u);
            changed = 1;
          }
        }
        // T2 transformation need skip the entry node!
        for (int v = 1; v < _size; v++) {
          if (coalesced[v] == 1)
            continue;
          // no need to check self-cycle v->v since it already erased by T1
          // transformation just now
          if (precessors[v].size() == 1) {
            T2(v);
            changed = 1;
          }
        }
        if (!changed)
          std::cout << "Done!\n";
      }
      int unCoalescedCount = getUncoalescedCount();
      std::cout << "#UnCollased Node is " << unCoalescedCount << "\n";
      return unCoalescedCount == 1;
    }

  private:
    const FlowGraph &_graph;
    int _size;
    std::vector<std::string> _coalesced_vertex_names;
    int *coalesced;
    std::vector<std::set<int>> precessors{};
    std::vector<std::set<int>> successors{};
  };

  class ReductionID {
  public:
    ReductionID(int highpt, int snumber) : _highpt(highpt), _snumber(snumber) {}

    bool operator<(const ReductionID &rhs) const {
      return ((_highpt > rhs._highpt) ||
              ((_highpt == rhs._highpt) && (_snumber < rhs._snumber)));
    }

    int _highpt;
    int _snumber;
  };

  explicit FlowGraph(int n) : _capacity(n), _size(0) {
    edges = new int[_capacity * _capacity]();
    dfn = new int[_capacity]();
    dfn_reverse = new int[_capacity]();
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
    delete[] dfn_reverse;
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

  void dumpReductionList() {
    std::cout << "### Dump Reduction List ###\n";
    for (auto it : reduction_list) {
      std::cout << vertex_names[it.second] << "(" << it.first._highpt << ", "
                << it.first._snumber << "), ";
    }
    std::cout << "\n";
  }

  void doReduceAsList(AdjacentView &adjview) {
    std::cout << "### Do Reduction As List Illustrates ###\n";
    for (auto it : reduction_list) {
      int u = it.second;
      if (adjview.checkSelfCycle(u)) {
        adjview.T1(u);
      }
      adjview.T2(u);
    }
    assert(adjview.getUncoalescedCount() == 1);
  }

  std::map<ReductionID, int> reduction_list;

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
  int *dfn_reverse;
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
  dfn_reverse[dfn[n]] = n;
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
  std::cout << "### Check Reducibility ###\n";
  reduction_list.clear();
  for (int w_id = _capacity - 1; w_id >= 0; w_id--) {
    int w = dfn_reverse[w_id];
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
              reduction_list.emplace(ReductionID{dfn[w], rpon[tmp_u]}, tmp_u);
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

  FlowGraph::AdjacentView g_adjview1{g};
  FlowGraph::AdjacentView g_adjview2{g};
  g_adjview1.checkCollapsibilitySlowly();
  std::cout << "\n";
  g_adjview1.dump();
  std::cout << "\n";

  g.dfs();
  g.checkReducibility();
  g.dumpDot("zsy_test1.dot");
  std::cout << "\n";
  g.dumpReductionList();
  std::cout << "\n";
  g.doReduceAsList(g_adjview2);
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

  FlowGraph::AdjacentView g_adjview{g};
  g_adjview.checkCollapsibilitySlowly();
  std::cout << "\n";
  g_adjview.dump();
  std::cout << "\n";

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
