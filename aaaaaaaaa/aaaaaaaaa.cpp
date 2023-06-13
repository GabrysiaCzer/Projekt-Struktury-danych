#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <climits>
#include <queue>
#include <fstream>

struct Edge {
    int src;
    int dest;
    int weight;
    Edge* next;

    Edge(int s, int d, int w) : src(s), dest(d), weight(w), next(nullptr) {}
};

struct ListNode {
    Edge data;
    ListNode* next;

    ListNode(const Edge& e) : data(e), next(nullptr) {}
};

struct List {
    ListNode* head;
    ListNode* tail;

    List() : head(nullptr), tail(nullptr) {}

    bool isEmpty() {
        return head == nullptr;
    }

    void insert(const Edge& e) {
        ListNode* newNode = new ListNode(e);
        if (isEmpty()) {
            head = newNode;
            tail = newNode;
        }
        else if (e.weight < head->data.weight) {
            newNode->next = head;
            head = newNode;
        }
        else {
            ListNode* curr = head;
            while (curr->next != nullptr && curr->next->data.weight < e.weight) {
                curr = curr->next;
            }

            newNode->next = curr->next;
            curr->next = newNode;

            if (newNode->next == nullptr) {
                tail = newNode;
            }
        }
    }

    Edge removeFront() {
        if (isEmpty()) {
            return Edge(-1, -1, -1);
        }
        ListNode* temp = head;
        Edge e = temp->data;
        head = head->next;
        if (head == nullptr) {
            tail = nullptr;
        }
        delete temp;
        return e;
    }
};

struct Graph {
    int numVertices;
    std::vector<List> adjList;

    Graph(int vertices) : numVertices(vertices), adjList(vertices) {}

    void addEdge(int src, int dest, int weight) {
        Edge newEdge(src, dest, weight);
        adjList[src].insert(newEdge);

        newEdge = Edge(dest, src, weight);
        adjList[dest].insert(newEdge);
    }

    std::vector<Edge> primMSTQueue() {
        std::vector<bool> visited(numVertices, false);
        std::vector<Edge> minimumSpanningTree;
        std::priority_queue<Edge*, std::vector<Edge*>, std::greater<Edge*>> pq;

        visited[0] = true;

        for (ListNode* node = adjList[0].head; node != nullptr; node = node->next) {
            pq.push(&node->data);
        }

        while (!pq.empty()) {
            Edge* current = pq.top();
            pq.pop();

            int nextWeight = current->weight;
            int nextVertex = current->dest;

            if (!visited[nextVertex]) {
                visited[nextVertex] = true;
                minimumSpanningTree.push_back(*current);

                for (ListNode* node = adjList[nextVertex].head; node != nullptr; node = node->next) {
                    if (!visited[node->data.dest]) {
                        pq.push(&node->data);
                    }
                }
            }
        }

        return minimumSpanningTree;
    }

    std::vector<Edge> primMSTim() {
        std::vector<bool> visited(numVertices, false);
        std::vector<int> parent(numVertices, -1);
        std::vector<int> key(numVertices, INT_MAX);

        List pq;

        int startVertex = 0;

        pq.insert(Edge(startVertex, startVertex, 0));
        key[startVertex] = 0;

        while (!pq.isEmpty()) {
            Edge minEdge = pq.removeFront();
            int u = minEdge.dest;
            visited[u] = true;

            ListNode* node = adjList[u].head;
            while (node != nullptr) {
                int v = node->data.dest;
                int weight = node->data.weight;

                if (!visited[v] && weight < key[v]) {
                    parent[v] = u;
                    key[v] = weight;
                    pq.insert(Edge(u, v, weight));
                }

                node = node->next;
            }

            int minWeight = INT_MAX;
            int minIndex = -1;
            for (int i = 0; i < numVertices; i++) {
                if (!visited[i] && key[i] < minWeight) {
                    minWeight = key[i];
                    minIndex = i;
                }
            }

            if (minIndex != -1) {
                pq.insert(Edge(parent[minIndex], minIndex, minWeight));
            }
        }

        std::vector<Edge> minimumSpanningTree;
        for (int i = 1; i < numVertices; i++) {
            minimumSpanningTree.push_back(Edge(parent[i], i, key[i]));
        }

        return minimumSpanningTree;
    }
};

std::vector<Edge> generateGraph(int numVertices) {
    std::vector<Edge> graph;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(1, 100);

    for (int i = 0; i < numVertices; i++) {
        for (int j = i + 1; j < numVertices; j++) {
            int weight = dist(gen);
            graph.emplace_back(i, j, weight);
        }
    }

    return graph;
}

int main() {
    std::vector<int> numVerticesList = { 50, 100, 500 };
    std::vector<double> densityList = { 0.25, 0.5, 0.75, 1.0 };

    std::ofstream outputFile("results.csv");
    outputFile << "NumVertices,Density,QueueTime,ImTime\n";

    for (int numVertices : numVerticesList) {
        for (double density : densityList) {
            for (int i = 0; i < 10; i++) {
                int numEdges = density * numVertices * (numVertices - 1) / 2;

                Graph g(numVertices);
                std::vector<Edge> graph = generateGraph(numVertices);

                for (const Edge& edge : graph) {
                    g.addEdge(edge.src, edge.dest, edge.weight);
                }

                auto start = std::chrono::high_resolution_clock::now();
                std::vector<Edge> mstQueue = g.primMSTQueue();
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> durationQueue = end - start;

                start = std::chrono::high_resolution_clock::now();
                std::vector<Edge> mstIm = g.primMSTim();
                end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> durationIm = end - start;

                outputFile << numVertices << "," << density << "," << durationQueue.count() << "," << durationIm.count() << "\n";
            }
        }
    }

    outputFile.close();

    std::cout << "Results saved to results.csv" << std::endl;

    return 0;
}
