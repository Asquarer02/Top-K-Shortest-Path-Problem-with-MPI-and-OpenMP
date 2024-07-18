#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <limits>
#include <string>
#include <queue>
#include<algorithm>
#include <ctime>
#include<mpi.h>
#include<omp.h>

using namespace std;

const double INF = numeric_limits<double>::max();

void printDistanceMatrix(const unordered_map<int, unordered_map<int, double>>& distance_matrix) {
  for (const auto& src_pair : distance_matrix) {
    int from_node = src_pair.first;
    cout << "Distances from node " << from_node << ":" << endl;

    bool has_reachable = false;
    for (const auto& dest_pair : distance_matrix.at(from_node)) {
      int to_node = dest_pair.first;
      double distance = dest_pair.second;

      if (distance < INF) {
        cout << "  To node " << to_node << ": " << distance << endl;
        has_reachable = true;
      }
    }

    if (!has_reachable) {
      cout << "  No reachable nodes from " << from_node << endl;
    }
  }
}

bool readHeader(ifstream& file, int& N, int& M) {
  string line;
  while (getline(file, line)) {
    if (line.find("Nodes:") != string::npos && line.find("Edges:") != string::npos) {
      istringstream iss(line);
      string temp;
      iss >> temp >> temp >> N >> temp >> M; 
      return true;
    }
  }
  return false; 
}

void readGraph(ifstream& file, vector<vector<int>>& matrix) {
    string line;
    while (getline(file, line)) {
        
        if (!line.empty() && line[0] != '#') {
            istringstream ss(line);
            int from, to, weight;
            ss >> from >> to >> weight;  

            if (ss.fail()) {
                weight = 1;  
                ss.clear();  
            }
            matrix.push_back({from, to, weight});
        }
    }
}



struct Path {
    int cost;
    vector<int> nodes;

    bool operator>(const Path& other) const {
        return cost > other.cost;
    }
};

vector<pair<int, int>> extractEdges(const vector<int>& path, const vector<vector<pair<int, int>>>& g) {
    vector<pair<int, int>> edges;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        int u = path[i];
        int v = path[i + 1];
        auto it = find_if(g[u].begin(), g[u].end(), [v](const pair<int, int>& edge) { return edge.first == v; });
        if (it != g[u].end()) {
            edges.emplace_back(u, v);
        }
    }
    return edges;
}

void yenAlgorithm(const vector<vector<int>>& matrix, int n, int src, int dest, int k) {
    
    vector<vector<pair<int, int>>> g(n + 1);

    #pragma omp parallel for
    for (int index = 0; index < matrix.size(); ++index) {
        const auto& edge = matrix[index];
        int from = edge[0];
        int to = edge[1];
        int weight = edge[2];
        #pragma omp critical
        g[from].push_back({to, weight});
    }

    auto dijkstra = [&g, n](int start, int end) -> pair<int, vector<int>> {
        vector<int> dist(n + 1, numeric_limits<int>::max());
        vector<int> prev(n + 1, -1);
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
        dist[start] = 0;
        pq.push({0, start});

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();
            if (!pq.empty() && dist[u] < pq.top().first) continue;

            for (const auto& adj : g[u]) {
                int v = adj.first;
                int weight = adj.second;
                if (dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    prev[v] = u;
                    pq.push({dist[v], v});
                }
            }
        }

        vector<int> path;
        for (int at = end; at != -1; at = prev[at]) {
            path.push_back(at);
        }
        reverse(path.begin(), path.end());
        return {dist[end], path};
    };

    vector<Path> A;
    priority_queue<Path, vector<Path>, greater<Path>> B;

    auto [initialCost, initialPath] = dijkstra(src, dest);
    A.push_back({initialCost, initialPath});

    for (int t = 1; t < k; ++t) {
        const vector<int>& lastPath = A[t-1].nodes;

        #pragma omp parallel for
        for (size_t i = 0; i < lastPath.size() - 1; ++i) {
            int spurNode = lastPath[i];
            vector<pair<int, pair<int, int>>> removedEdges;

            #pragma omp critical
            {
                for (size_t j = i + 1; j < lastPath.size(); ++j) {
                    int u = lastPath[j-1];
                    int v = lastPath[j];
                    auto& edges = g[u];
                    for (auto it = edges.begin(); it != edges.end();) {
                        if (it->first == v) {
                            removedEdges.push_back({u, *it});
                            it = edges.erase(it);
                        } else {
                            ++it;
                        }
                    }
                }
            }

            auto [spurCost, spurPath] = dijkstra(spurNode, dest);
            if (!spurPath.empty() && spurCost < numeric_limits<int>::max()) {
                vector<int> totalPath(lastPath.begin(), lastPath.begin() + i + 1);
                totalPath.insert(totalPath.end(), spurPath.begin() + 1, spurPath.end());

                int rootCost = 0;
                for (size_t j = 1; j <= i; ++j) {
                    int u = lastPath[j-1];
                    int v = lastPath[j];
                    auto it = find_if(g[u].begin(), g[u].end(), [v](const pair<int, int>& edge) { return edge.first == v; });
                    if (it != g[u].end()) {
                        rootCost += it->second;
                    }
                }

                int totalCost = rootCost + spurCost;

                Path newPath{totalCost, totalPath};

                #pragma omp critical
                {
                    if (find_if(A.begin(), A.end(), [&newPath](const Path& p) {
                        return p.nodes == newPath.nodes;
                    }) == A.end()) {
                        B.push(newPath);
                    }
                }
            }

            
            #pragma omp critical
            {
                for (const auto& edge : removedEdges) {
                    g[edge.first].push_back(edge.second);
                }
            }
        }

        if (B.empty()) {
            cout << "No more paths exist." << endl;
            break;
        }

        #pragma omp critical
        {
            A.push_back(B.top());
            B.pop();
        }
    }

    sort(A.begin(), A.end(), [](const Path& a, const Path& b) {
        return a.cost < b.cost;
    });

    cout << "The " << k << " shortest paths from node " << src << " to node " << dest << " are:" << endl;
    for (const auto& path : A) {
        if(path.cost==2147483647){
            cout<<"No path"<<endl;
        }
        else
        {cout << "Path Cost: " << path.cost << " Route: ";
        for (int node : path.nodes) {
            cout << node << " ";
        }
        cout << endl;}
    }
}

int main(int argc, char* argv[]) {
    
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double startTime = MPI_Wtime();
    ifstream file("Email-EuAll.txt");
    if (!file.is_open()) {
        cout << "Failed to open file." << endl;
        return -1;
    }

    int N = 0, M = 0; 

    if (!readHeader(file, N, M)) {
        cout << "Invalid or missing header." << endl;
        return -1;
    }
    if(rank==0){
    cout << "Number of nodes: " << N << endl;
    cout << "Number of edges: " << M << endl;}
     vector<vector<int>> matrix;
    //cout<<"Hi"<<endl;
    readGraph(file, matrix);
    file.close();
    
    int K = 3; 
    int numPairsPerProcess = 10 / size; 
    int remainder = 10 % size; 

    srand(time(0) + rank); 

    for (int i = 0; i < numPairsPerProcess + (rank < remainder ? 1 : 0); ++i) {
        int src, dest;
        do {
            src = rand() % N + 1; 
            dest = rand() % N + 1; 
        } while (src == dest); 

        cout << "Process: " << rank << " explores src: " << src << " dest: " << dest << endl;
        yenAlgorithm(matrix, N, src, dest, K);
    }
    
    double endTime = MPI_Wtime();
    MPI_Finalize();
    
    if (rank == 0) { 
        cout << "Total execution time: " << (endTime - startTime) << " seconds." << endl;
    }

    return 0;
}
