#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <climits>

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

    void primMST() {
        std::vector<bool> visited(numVertices, false);
        std::vector<int> parent(numVertices, -1);
        std::vector<int> key(numVertices, INT_MAX);

        struct Vertex {
            int index;
            int weight;
            Vertex(int i, int w) : index(i), weight(w) {}
        };

        auto compare = [](const Vertex& a, const Vertex& b) {
            return a.weight > b.weight;
        };

        List pq;

        int startVertex = 0;

        pq.insert(Edge(startVertex, startVertex, 0));
        key[startVertex] = 0;

        while (!pq.isEmpty()) {
            int u = pq.removeFront().dest;
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
        }

        std::cout << "Minimal Spanning Tree:\n";
        for (int i = 1; i < numVertices; i++) {
            std::cout << parent[i] << " - " << i << std::endl;
        }
    }
    void primMSTim() {
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
            // Wybieranie krawędzi o minimalnej wadze spośród nieodwiedzonych wierzchołków
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

        std::cout << "Minimal Spanning Tree:\n";
        for (int i = 1; i < numVertices; i++) {
            std::cout << parent[i] << " - " << i << std::endl;
        }
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
    int numVertices = 10;

  

    Graph g(numVertices);

    g.addEdge(0, 1, 4);
    g.addEdge(0, 2, 8);
    g.addEdge(1, 2, 11);
    g.addEdge(1, 3, 3);
    g.addEdge(2, 4, 7);
    g.addEdge(2, 5, 1);
    g.addEdge(3, 4, 2);
    g.addEdge(3, 6, 5);
    g.addEdge(4, 5, 6);
    g.addEdge(4, 7, 12);
    g.addEdge(5, 8, 13);
    g.addEdge(6, 7, 14);
    g.addEdge(6, 9, 9);
    g.addEdge(7, 8, 15);
    g.addEdge(8, 9, 10);

    auto start = std::chrono::high_resolution_clock::now();

    g.primMST();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
    g.primMSTim();
    return 0;
}

