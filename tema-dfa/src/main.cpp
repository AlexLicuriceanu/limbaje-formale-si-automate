#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <limits>
#include <queue>

using namespace std;

auto readInputFile(const string& inputFile) {

    ifstream file(inputFile);

    int subtask;
    int n, sigmaLength, m, k;
    int qx;
    vector<int> finalStates;
    vector<vector<int>> transitionMatrix;
    
    file >> subtask >> n >> sigmaLength >> m >> k >> qx;
    
    qx--;

    finalStates.resize(m);
    transitionMatrix.resize(n, vector<int>(sigmaLength));

    for (int i = 0; i < m; i++) {
        int aux;
        file >> aux;
        aux--;

        finalStates[i] = aux;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < sigmaLength; j++) {

            int aux;
            file >> aux;
            aux--;

            transitionMatrix[i][j] = aux;
        }
    }

    file.close();

    return make_tuple(subtask, n, sigmaLength, m, k, qx, finalStates, transitionMatrix);
}

void writeOutputFile(const string& outputFile, int wordLength, string word) {
    
    ofstream file;
    file.open(outputFile);

    file << wordLength << "\n" << word;
    file.close();
}

vector<vector<int>> transpose(vector<vector<int>>& matrix) {

    int n = matrix.size();
    vector<vector<int>> transposed(n, vector<int>(n, 0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}

void dfs(int node, const vector<vector<int>>& matrix, vector<bool>& visited) {
    visited[node] = true;

    for (int neighbor = 0; neighbor < matrix[node].size(); neighbor++) {
        if (matrix[node][neighbor] == 1 && !visited[neighbor]) {
            dfs(neighbor, matrix, visited);
        }
    }
}

vector<bool> getReachableNodes(const vector<vector<int>>& matrix, const vector<int>& finalStates) {

    int n = matrix.size();
    vector<bool> reachable(n, false);

    // Mark all reachable nodes from final states.
    for (int state : finalStates) {
        dfs(state, matrix, reachable);
    }

    return reachable;
}

pair<int, string> findWord(int subtask, int n, int sigmaLength, int m, int k, int qx, vector<int>& finalStates, vector<vector<int>>& transitionMatrix) {

    // 0. Implementation only for task 1.
    if (subtask != 1) {
        return make_pair(-1, "");
    }

    // 1. Set up adjacency matrix.
    vector<vector<int>> matrix(n, vector<int>(n, 0));

    for (int i = 0; i < transitionMatrix.size(); i++) {
        for (int j = 0; j < transitionMatrix[i].size(); j++) {
            matrix[i][transitionMatrix[i][j]] = 1; 
        }
    }

    // 1.1. Compute transposed adjacency matrix.
    vector<vector<int>> transposedMatrix = transpose(matrix);

    // 2. Get reachable nodes from final states.
    vector<bool> reachable = getReachableNodes(transposedMatrix, finalStates);

    // 2.1. Remove unreachable nodes and their transitions from matrix and transposedMatrix.
    vector<vector<int>> prunedMatrix(n, vector<int>(n, 0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (reachable[i] && reachable[j]) {
                prunedMatrix[i][j] = matrix[i][j];
            }
        }
    }

    // 2.2. Transpose the pruned matrix.
    vector<vector<int>> prunedTransposedMatrix = transpose(prunedMatrix);

    // 3. Compute distances vector. Run BFS and compute minimum distances.
    vector<vector<int>> distances(finalStates.size(), vector<int>(n, numeric_limits<int>::max()));

    for (int i = 0; i < finalStates.size(); i++) {

        int finalState = finalStates[i];

        vector<int> distance(n, numeric_limits<int>::max());
        distance[finalState] = 0;

        queue<int> bfsQueue;
        bfsQueue.push(finalState);

        while (!bfsQueue.empty()) {
            int current = bfsQueue.front();
            bfsQueue.pop();

            for (int neighbor = 0; neighbor < n; neighbor++) {
                if (prunedTransposedMatrix[current][neighbor] == 1 && distance[neighbor] == numeric_limits<int>::max()) {
                    distance[neighbor] = distance[current] + 1;
                    bfsQueue.push(neighbor);
                }
            }
        }

        distances[i] = distance;
    }

    // 3.1. Compute minimum distance vector.
    vector<int> d(n);

    for (int j = 0; j < n; j++) {
        
        int minimum = distances[0][j];

        for (int i = 0; i < finalStates.size(); i++) {
            if (minimum > distances[i][j]) {
                minimum = distances[i][j];
            }
        }

        d[j] = minimum;
    }


    // 4. Compute the DP table.
    vector<vector<int>> dp(k + 1, vector<int>(n, -1));

    dp[0][qx] = qx;

    // O(k*n*sigmaLength)

    for (int h = 1; h <= k; h++) {
        for (int i = 0; i < n; i++) {

            if (dp[h-1][i] != -1) {

                for (int j : transitionMatrix[i]) {
                    if (reachable[i] && reachable[j]) {
                        dp[h][j] = i;
                    }
                }
            }
        }
    }
 
    // 4.1. Find node at level k with shortest distance to a final node.
    int node1 = -1;
    int node2 = -1;
    int minimumDistance = numeric_limits<int>::max();

    for (int i = 0; i < n; i++) {
        if (dp[k][i] != -1 && d[i] < minimumDistance) {
            node1 = i;
            node2 = dp[k][i];

            minimumDistance = d[i];
        }
    }
    
    // 5. Check if there is no solution.
    if (node1 == -1) {
        return make_pair(-1, "");
    }

    // 5.1. Compute the prefix, walk the DP matrix backwards.
    int selectedNode = node1;

    string prefix = "";

    for (int h = k; h > 0; h--) {
        
        int t = 0;

        // Find the transition index.
        for (int i = 0; i < sigmaLength; i++) {
            if (transitionMatrix[node2][i] == node1) {
                t = i;
                break;
            }
        }

        prefix = (char)('a' + t) + prefix;

        node1 = node2;
        node2 = dp[h-1][node2];
    }

    // If the path to a final node is exactly k, exit.
    if (d[selectedNode] == 0) {
        return make_pair(prefix.size(), prefix);
    }


    // 5.2. Compute the suffix, BFS until a final node is reached.
    string suffix = "";

    int remainingDistance = d[selectedNode];
    
    while (remainingDistance > 0) {

        // Find neighbor of node1 with minimum distance.
        int node2 = transitionMatrix[selectedNode][0];
        int t = 0;

        for (int i = 0; i < sigmaLength; i++) {
            if (d[transitionMatrix[selectedNode][i]] < d[node2]) {
                node2 = transitionMatrix[selectedNode][i];
                t = i;
            }
        }

        suffix = suffix + (char)('a' + t);

        remainingDistance--;
        selectedNode = node2;
    }

    return make_pair(prefix.size() + suffix.size(), prefix + suffix);
}

int main() {

    // Read input data.
    auto [subtask, n, sigmaLength, m, k, qx, finalStates, transitionMatrix] = readInputFile("input.txt");

    // Run the algorithm.
    auto [wordLength, word] = findWord(subtask, n, sigmaLength, m, k, qx, finalStates, transitionMatrix);

    // Write output data.
    writeOutputFile("output.txt", wordLength, word);
    
    return 0;
}