#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <climits>
#include <iomanip>
#include <algorithm>

using namespace std;

#define ll long long
#define INF LLONG_MAX
#define EXIT (-2)
#define PASS (-3)
#define OBSTRUCTION (-9)
#define BLANK "\n\n\n"


class Graph;

class Vertex;

class Node;


void initial_vertices();

void load_data_to_graph();

void create_connection();

// calculate

vector<vector<ll>> calculate_shortest_path_to_exit();

vector<ll> pick_best_exit(const vector<vector<ll>> &);

vector<ll> pick_closest_exit(const vector<vector<ll>> &);

vector<ll> pick_from_exit();

vector<vector<vector<Node>>> pick_from_A_star();

vector<vector<ll>> A_star_to_Dijkstra(vector<vector<vector<Node>>> &);

// print

void print_distance_to_exit(const vector<vector<ll>> &);

void print_index_of_vertices(const vector<vector<string>> &);

void print_index_of_vertices_for_file(const vector<vector<string>> &);

void print_exit_of_vertices(const vector<ll> &);

void print_exit_with_count(const vector<ll> &);

void print_exit_with_astar(const vector<vector<vector<Node>>> &);


class Node {
public:
    ll index, g, h;
    Node *parent;


    Node(ll index, ll g = 0, ll h = 0) : index(index), g(g), h(h), parent(nullptr) {};


    Node(ll index, ll g, ll h, Node &parent) : index(index), g(g), h(h), parent(&parent) {};


    ll length() const {
        return g + h;
    }


    static ll calculate_heuristic(ll x, ll y, ll target_x, ll target_y) {
        return llabs(x - target_x) + llabs(y - target_y);
    }
};


class Graph {
public:
    ll row, col, pass_count, exit_count, index_count;
    vector<ll> pass_index, exit_index;
    vector<vector<ll>> data;
    vector<vector<Vertex>> vertices;
    map<ll, Vertex *> index_to_vertex;

    Graph() : row(0), col(0), pass_count(0), exit_count(0), index_count(0), pass_index({}), exit_index({}), data({}),
              vertices({}), index_to_vertex({}) {}

    void print_with_symbol() {
        string graph;
        for (ll i = 0; i < row; i++) {
            for (ll j = 0; j < col; j++) {
                switch (data[i][j]) {
                    case PASS:
                        graph.push_back('#');
                        break;
                    case EXIT:
                        graph.push_back('O');
                        break;
                    case OBSTRUCTION:
                        graph.push_back('X');
                        break;
                }
            }
            graph.push_back('\n');
        }
        cout << graph << BLANK;
    }

    vector<vector<string>> print_with_circle() {
        // build vertex
        vector<vector<string>> vertex;
        for (ll i = 0; i < row; i++) {
            vertex.emplace_back(vector<string>());
            for (ll j = 0; j < col; j++) {
                switch (data[i][j]) {
                    case PASS:
                        vertex[i].push_back("°?");
                        break;
                    case EXIT:
                        vertex[i].push_back("°?");
                        break;
                    case OBSTRUCTION:
                        if (!vertex[i].empty() and vertex[i].back() == "°X")
                            vertex[i].back() = "°@";
                        vertex[i].push_back("°@");
                        break;
                }
                vertex[i].push_back(vertex[i].back() != "°@" ? "°X" : "°@");
            }
            if (vertex[i].back() == "°X")
                vertex[i].back() = "\n";
            else
                vertex[i].push_back("\n");
        }

        // build edge
        vector<vector<string>> edge;
        for (ll i = 0; i < row; i++) {
            edge.emplace_back(vector<string>());
            for (ll j = 0; j < vertex[i].size() - 1; j++) {
                if (i + 1 < row and (vertex[i][j] == "°?" or vertex[i][j] == "°?") and
                    (vertex[i + 1][j] == "°?" or vertex[i + 1][j] == "°?"))
                    edge[i].push_back("°U");
                else
                    edge[i].push_back("°@");
            }
            edge[i].push_back("\n");
        }

        // combine vertex and edge to graph
        vector<vector<string>> res;
        for (int i = 0; i < vertex.size(); i++) {
            res.push_back(vertex[i]);
            res.push_back(edge[i]);
        }
        res.pop_back();

        // print graph
        for (auto &x: res) for (auto &y: x) cout << y;
        cout << BLANK;

        return res;
    }


    vector<vector<string>> print_with_circle_for_file() {
        // build vertex
        vector<vector<string>> vertex;
        for (ll i = 0; i < row; i++) {
            vertex.emplace_back(vector<string>());
            for (ll j = 0; j < col; j++) {
                switch (data[i][j]) {
                    case PASS:
                        vertex[i].push_back("°?");
                        break;
                    case EXIT:
                        vertex[i].push_back("°?");
                        break;
                    case OBSTRUCTION:
                        if (!vertex[i].empty() and vertex[i].back() == "°X")
                            vertex[i].back() = " ";
                        vertex[i].push_back(" ");
                        break;
                }
                vertex[i].push_back(vertex[i].back() != " " ? "°X" : " ");
            }
            if (vertex[i].back() == "°X")
                vertex[i].back() = "\n";
            else
                vertex[i].push_back("\n");
        }

        // build edge
        vector<vector<string>> edge;
        for (ll i = 0; i < row; i++) {
            edge.emplace_back(vector<string>());
            for (ll j = 0; j < vertex[i].size() - 1; j++) {
                if (i + 1 < row and (vertex[i][j] == "°?" or vertex[i][j] == "°?") and
                    (vertex[i + 1][j] == "°?" or vertex[i + 1][j] == "°?"))
                    edge[i].push_back("|");
                else
                    edge[i].push_back(" ");
            }
            edge[i].push_back("\n");
        }

        // combine vertex and edge to graph
        vector<vector<string>> res;
        for (int i = 0; i < vertex.size(); i++) {
            res.push_back(vertex[i]);
            res.push_back(edge[i]);
        }
        res.pop_back();

        // print graph
        for (auto &x: res) for (auto &y: x) cout << y;
        cout << BLANK;

        return res;
    }
} G;


class Vertex {
public:
    ll index;
    ll type;
    ll x, y;
    map<ll, ll> adj; // index : distance

    Vertex() : index(-1), type(OBSTRUCTION), x(0), y(0) {}

    Vertex(ll index, ll type, ll x, ll y) : index(index), type(type), x(x), y(y), adj() {}

    bool connectable() const {
        return type != OBSTRUCTION;
    }

    void addEdge(ll v, ll weight) {
        adj[v] = weight;
    }

    // from wikipedia
    vector<ll> dijkstra() {
        priority_queue<pair<ll, ll>, vector<pair<ll, ll>>, greater<>> pq;

        vector<ll> distance(G.index_count, INF);
        vector<bool> visited(G.index_count, false);

        pq.push(make_pair(0, this->index));
        distance[this->index] = 0;

        while (!pq.empty()) {
            ll u = pq.top().second;
            pq.pop();

            if (visited[u]) continue;
            visited[u] = true;


            for (auto &i: G.index_to_vertex[u]->adj) {
                ll v = i.first, weight = i.second;

                // relax
                if (distance[v] > distance[u] + weight) {
                    distance[v] = distance[u] + weight;
                    pq.push(make_pair(distance[v], v));
                }
            }
        }
        return distance;
    }


    vector<Node> a_star(ll target) {
        auto cmp = [](Node a, Node b) {
            return a.length() > b.length();
        };
        ll target_x = G.index_to_vertex[target]->x, target_y = G.index_to_vertex[target]->y;
        priority_queue<Node, vector<Node>, decltype(cmp)> open(cmp);
        set<ll> closed;
        set<ll> count;
        Node *point = new Node(index, 0, llabs(target_x - x) + llabs(target_y - y));
        open.push(*point);
        while (!open.empty()) {
            Node tmp = open.top();
            point = new Node(tmp.index, tmp.g, tmp.h, *tmp.parent);
            open.pop();
            if (closed.count(point->index) or !G.index_to_vertex[point->index]->connectable())
                continue;

            ll point_x = G.index_to_vertex[point->index]->x, point_y = G.index_to_vertex[point->index]->y;
            closed.insert(point->index);

            if (point->index == target)
                break;

            if (point_x > 0 and !count.count(G.vertices[point_x - 1][point_y].index)) {
                count.insert(G.vertices[point_x - 1][point_y].index);
                open.push(Node(G.vertices[point_x - 1][point_y].index, point->g + 1,
                               Node::calculate_heuristic(point_x - 1, point_y, target_x, target_y), *point));
            }
            if (point_x < G.row - 1 and !count.count(G.vertices[point_x + 1][point_y].index)) {
                count.insert(G.vertices[point_x + 1][point_y].index);
                open.push(Node(G.vertices[point_x + 1][point_y].index, point->g + 1,
                               Node::calculate_heuristic(point_x + 1, point_y, target_x, target_y), *point));
            }
            if (point_y > 0 and !count.count(G.vertices[point_x][point_y - 1].index)) {
                count.insert(G.vertices[point_x][point_y - 1].index);
                open.push(Node(G.vertices[point_x][point_y - 1].index, point->g + 1,
                               Node::calculate_heuristic(point_x, point_y - 1, target_x, target_y), *point));
            }
            if (point_y < G.col - 1 and !count.count(G.vertices[point_x][point_y + 1].index)) {
                count.insert(G.vertices[point_x][point_y + 1].index);
                open.push(Node(G.vertices[point_x][point_y + 1].index, point->g + 1,
                               Node::calculate_heuristic(point_x, point_y + 1, target_x, target_y), *point));
            }
        }
        vector<Node> result;
        Node *p = point;
        while (p) {
            result.push_back(*p);
            p = p->parent;
        }
        reverse(result.begin(), result.end());
        return result;
    }
};


int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // init
    load_data_to_graph();
    initial_vertices();

    // print for input
    // G.print_with_symbol();
    // G.print_with_circle();
    // G.print_with_circle_for_file();

    // build
    create_connection();

    // calculate
    vector<vector<ll>> vertices_distance_to_exit = calculate_shortest_path_to_exit();
    // vector<ll> result = pick_best_exit(vertices_distance_to_exit);
    // vector<ll> result = pick_closest_exit(vertices_distance_to_exit);
    // vector<ll> result = pick_from_exit();
    vector<vector<vector<Node>>> result = pick_from_A_star();

    // print for output
    // print_distance_to_exit(vertices_distance_to_exit);
    // print_index_of_vertices_for_file(G.print_with_circle_for_file());
    // print_exit_of_vertices(result);
    // print_exit_with_count(result);
    print_distance_to_exit(A_star_to_Dijkstra(result));
    // print_exit_with_astar(result);

    return 0;
}


void load_data_to_graph() {
    ll row, col;
    cin >> row >> col;

    vector<vector<ll>> input_data(row, vector<ll>(col, 0));
    for (ll i = 0; i < row; i++)
        for (ll j = 0; j < col; j++)
            cin >> input_data[i][j];

    G.row = row, G.col = col, G.data = input_data;
}


void initial_vertices() {
    G.vertices.assign(G.row, vector<Vertex>(G.col, Vertex()));
    for (ll i = 0, count = 0, index, type; i < G.row; i++) {
        for (ll j = 0; j < G.col; j++) {
            type = G.data[i][j];
            switch (type) {
                case PASS:
                    index = count++;
                    G.pass_index.push_back(index);
                    break;
                case EXIT:
                    index = count++;
                    G.exit_index.push_back(index);
                    break;
                case OBSTRUCTION:
                    index = -1;
                    break;
                default:
                    break;
            }
            G.vertices[i][j] = Vertex(index, type, i, j);
            G.index_to_vertex[index] = &G.vertices[i][j];
        }
    }
    G.pass_count = G.pass_index.size(), G.exit_count = G.exit_index.size();
    G.index_count = G.pass_count + G.exit_count;
}


void create_connection() {
    for (int i = 0; i < G.row; i++) {
        for (int j = 0; j < G.col; j++) {
            if (G.vertices[i][j].connectable()) {
                if (i and G.vertices[i - 1][j].connectable()) {
                    G.vertices[i][j].addEdge(G.vertices[i - 1][j].index, 1);
                    G.vertices[i - 1][j].addEdge(G.vertices[i][j].index, 1);
                }
                if (i + 1 < G.row and G.vertices[i + 1][j].connectable()) {
                    G.vertices[i][j].addEdge(G.vertices[i + 1][j].index, 1);
                    G.vertices[i + 1][j].addEdge(G.vertices[i][j].index, 1);
                }
                if (j and G.vertices[i][j - 1].connectable()) {
                    G.vertices[i][j].addEdge(G.vertices[i][j - 1].index, 1);
                    G.vertices[i][j - 1].addEdge(G.vertices[i][j].index, 1);
                }
                if (j + 1 < G.col and G.vertices[i][j + 1].connectable()) {
                    G.vertices[i][j].addEdge(G.vertices[i][j + 1].index, 1);
                    G.vertices[i][j + 1].addEdge(G.vertices[i][j].index, 1);
                }
            }
        }
    }
}


// return each vertex's distance to exit.
vector<vector<ll>> calculate_shortest_path_to_exit() {
    vector<vector<ll>> res;
    for (auto &i: G.exit_index) {
        ll x = G.index_to_vertex[i]->x, y = G.index_to_vertex[i]->y;
        res.push_back(G.vertices[x][y].dijkstra());
    }
    return res;
}


vector<ll> pick_best_exit(const vector<vector<ll>> &ar) {
    vector<ll> distance(G.index_count, -1);
    for (ll index = 0; index < G.index_count; index++) {
        ll tmp = INF, exit = -1;
        for (ll i = 0; i < ar.size(); i++) {
            if (ar[i][index] < tmp) {
                tmp = ar[i][index], exit = i;
            }
        }
        distance[index] = exit;
    }
    return distance;
}

vector<ll> pick_closest_exit(const vector<vector<ll>> &ar) {
    vector<ll> res(G.index_count, -1);
    vector<vector<pair<ll, ll>>> table(G.exit_count);
    for (ll i = 0; i < G.exit_count; i++) {
        for (ll j = 0; j < G.index_count; j++)
            table[i].push_back(make_pair(j, ar[i][j]));
        sort(table[i].begin(), table[i].end(), [](pair<ll, ll> &a, pair<ll, ll> &b) {
            return a.second == b.second ? a.first > b.first : a.second > b.second;
        });
    }
    for (ll index = 0; index < G.index_count; index++) {
        for (ll exit_index = 0; exit_index < G.exit_count; exit_index++) {
            while (!table[exit_index].empty()) {
                pair<ll, ll> p = table[exit_index].back();
                table[exit_index].pop_back();
                if (res[p.first] != -1)
                    continue;
                res[p.first] = exit_index;
                break;
            }
        }
    }
    return res;
}

vector<ll> pick_from_exit() {
    vector<ll> res(G.index_count, -1);
    vector<queue<ll>> table(G.exit_count);
    for (ll i = 0; i < G.exit_index.size(); i++)
        table[i].push(G.exit_index[i]);
    for (ll i = 0; i < G.index_count; i++) {
        for (ll exit_index = 0; exit_index < G.exit_count; exit_index++) {
            while (!table[exit_index].empty()) {
                ll index = table[exit_index].front();
                table[exit_index].pop();
                if (!G.index_to_vertex[index]->connectable())
                    res[index] = G.index_to_vertex[index]->type;
                if (res[index] != -1)
                    continue;
                res[index] = exit_index;
                if (G.index_to_vertex[index]->x > 0)
                    table[exit_index].push(
                            G.vertices[G.index_to_vertex[index]->x - 1][G.index_to_vertex[index]->y].index);
                if (G.index_to_vertex[index]->x < G.row - 1)
                    table[exit_index].push(
                            G.vertices[G.index_to_vertex[index]->x + 1][G.index_to_vertex[index]->y].index);
                if (G.index_to_vertex[index]->y > 0)
                    table[exit_index].push(
                            G.vertices[G.index_to_vertex[index]->x][G.index_to_vertex[index]->y - 1].index);
                if (G.index_to_vertex[index]->y < G.col - 1)
                    table[exit_index].push(
                            G.vertices[G.index_to_vertex[index]->x][G.index_to_vertex[index]->y + 1].index);
                break;
            }
        }
    }
    return res;
}


vector<vector<vector<Node>>> pick_from_A_star() {
    vector<vector<vector<Node>>> res;
    for (auto &index: G.exit_index) {
        res.push_back(vector<vector<Node>>());
        for (auto &target: G.pass_index) {
            res.back().push_back(G.index_to_vertex[index]->a_star(target));
        }
    }
    return res;
}

vector<vector<ll>> A_star_to_Dijkstra(vector<vector<vector<Node>>>& ar) {
    vector<vector<ll>> res;
    for (auto& x : ar) {
        res.push_back(vector<ll>());
        for (auto& y : x) {
            res.back().push_back(y.size() - 1);
        }
    }
    return res;
}


void print_distance_to_exit(const vector<vector<ll>> &ar) {
    for (ll i = 0; i < ar.size(); i++) {
        for (auto &x: ar[i])
            cout << setw(2) << setfill(' ') << (x == INF ? -1 : x) << ' ';
        cout << " ---------- from " << G.exit_index[i] << " exit" << "\n";
    }
    cout << BLANK;
}


void print_index_of_vertices(const vector<vector<string>> &ar) {
    ll count = 0;
    for (const auto &i: ar) {
        for (const auto &j: i) {
            if (j == "\n" or j == "°U" or j == "°X")
                cout << j;
            else if (j == "°?" or j == "°?") {
                cout << setw(3) << setfill(' ') << count;
                count++;
            } else
                cout << "°@";
        }
    }
    cout << BLANK;
}


void print_index_of_vertices_for_file(const vector<vector<string>> &ar) {
    ll count = 0;
    for (const auto &i: ar) {
        for (const auto &j: i) {
            if (j == "\n" or j == "|" or j == "°X")
                cout << j;
            else if (j == "°?" or j == "°?") {
                cout << count++;
            } else
                cout << " ";
        }
    }
    cout << BLANK;
}


void print_exit_of_vertices(const vector<ll> &ar) {
    ll digit = 0;
    map<ll, ll> mp;
    for (const auto &i: G.vertices) {
        for (auto &j: i) {
            if (j.type == OBSTRUCTION)
                cout << OBSTRUCTION << ' ';
            else {
                if (!mp.count(ar[j.index]))
                    mp[ar[j.index]] = digit++;
                cout << mp[ar[j.index]] << ' ';
            }
        }
    }
    cout << BLANK;
}


void print_exit_with_count(const vector<ll> &ar) {
    map<ll, ll> mp;
    for (auto &i: ar)
        if (i >= 0)
            mp[i]++;

    for (auto &i: mp)
        cout << "出口 " << i.first << " 共有 " << i.second << " 個點\n";
}

void print_exit_with_astar(const vector<vector<vector<Node>>> &ar) {
    for (auto &x: ar) {
        for (auto &y: x) {
            cout << y.size() << '\n';
//            for (auto node : y) {
//                 cout << node.index << ' ';
//            }
//            cout << '\n';
        }
        cout << '\n';
    }
}
