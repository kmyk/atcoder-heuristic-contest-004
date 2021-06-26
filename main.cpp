#include <bits/stdc++.h>
#define REP(i, n) for (int i = 0; (i) < (int)(n); ++ (i))
#define REP3(i, m, n) for (int i = (m); (i) < (int)(n); ++ (i))
#define REP_R(i, n) for (int i = (int)(n) - 1; (i) >= 0; -- (i))
#define REP3R(i, m, n) for (int i = (int)(n) - 1; (i) >= (int)(m); -- (i))
#define ALL(x) std::begin(x), std::end(x)
using namespace std;

class xor_shift_128 {
public:
    typedef uint32_t result_type;
    xor_shift_128(uint32_t seed = 42) {
        set_seed(seed);
    }
    void set_seed(uint32_t seed) {
        a = seed = 1812433253u * (seed ^ (seed >> 30));
        b = seed = 1812433253u * (seed ^ (seed >> 30)) + 1;
        c = seed = 1812433253u * (seed ^ (seed >> 30)) + 2;
        d = seed = 1812433253u * (seed ^ (seed >> 30)) + 3;
    }
    uint32_t operator() () {
        uint32_t t = (a ^ (a << 11));
        a = b; b = c; c = d;
        return d = (d ^ (d >> 19)) ^ (t ^ (t >> 8));
    }
    static constexpr uint32_t max() { return numeric_limits<result_type>::max(); }
    static constexpr uint32_t min() { return numeric_limits<result_type>::min(); }
private:
    uint32_t a, b, c, d;
};

constexpr int N = 20;
constexpr int D = 8;
constexpr int LEN_MIN = 2;
constexpr int LEN_MAX = 12;

array<array<char, N>, N> get_empty_board() {
    array<array<char, N>, N> f;
    REP (y, N) {
        fill(ALL(f[y]), '.');
    }
    return f;
}

int64_t calculate_score(int m, int c, int d) {
    return c < m ? 1.0e8 * c / m : 1.0e8 * 2 * N * N / (2 * N * N - d);
}

void print_field(ostream &out, const array<array<char, N>, N>& f) {
}

struct trie_node {
    array<trie_node *, D> children;
    vector<int> indices;
};

trie_node *construct_trie(const vector<pair<string, int>> &s) {
    if (s.empty()) {
        return nullptr;
    }
    array<vector<pair<string, int>>, D> ts;
    vector<int> indices;
    for (auto &s_i : s) {
        if (s_i.first.empty()) {
            indices.push_back(s_i.second);
        } else {
            ts[s_i.first[0] - 'A'].emplace_back(s_i.first.substr(1), s_i.second);
        }
    }
    array<trie_node *, D> children;
    REP (d, D) {
        children[d] = construct_trie(ts[d]);
    }
    return new trie_node((trie_node) {children, indices});
}

trie_node *construct_trie(const vector<string> &s) {
    vector<pair<string, int>> t(s.size());
    REP (i, s.size()) {
        t[i] = make_pair(s[i], i);
    }
    return construct_trie(t);
}

const trie_node *peek_trie(char c, const trie_node *trie) {
    if (trie == nullptr) {
        return nullptr;
    }
    if (c == '.') {
        return nullptr;
    }
    return trie->children[c - 'A'];
}

inline int modadd(int a, int b, int m) {
    int c = a + b;
    return c >= m ? c - m : c;
}

inline int modsub(int a, int b, int m) {
    int c = a - b;
    return c < 0 ? c + m : c;
}

template <class RandomEngine>
array<array<char, N>, N> solve(const int m, const vector<string> &s, RandomEngine& gen, chrono::high_resolution_clock::time_point clock_end) {
    chrono::high_resolution_clock::time_point clock_begin = chrono::high_resolution_clock::now();

    int len_min = INT_MAX;
    int len_max = INT_MIN;
    for (const string &s_i : s) {
        len_min = min<int>(len_min, s_i.length());
        len_max = max<int>(len_max, s_i.length());
    }

    array<array<char, N>, N> ans = get_empty_board();
    int ans_c = 0;
    int ans_d = N * N;
    int64_t highscore = calculate_score(m, ans_c, ans_d);

    const trie_node *trie = construct_trie(s);

    array<array<char, N>, N> cur = get_empty_board();
    int cur_c = 0;
    int cur_d = N * N;
    vector<int> used(m);

    auto update = [&](int y, int x, char c) {
        assert (0 <= y and y < N);
        assert (0 <= x and x < N);
        if (cur[y][x] == c) {
            return;
        }

        REP (is_hr, 2) {
            REP3 (offset, - len_max + 1, 0 + 1) {
                const trie_node *ptr = trie;
                REP (i, len_max) {
                    ptr = peek_trie(is_hr ? cur[y][modadd(modsub(x, - offset, N), i, N)] : cur[modadd(modsub(y, - offset, N), i, N)][x], ptr);
                    if (ptr == nullptr) {
                        break;
                    }
                    if (offset + i >= 0) {
                        for (int j : ptr->indices) {
                            used[j] -= 1;
                            if (not used[j]) {
                                cur_c -= 1;
                            }
                        }
                    }
                }
            }
        }

        if (cur[y][x] == '.') {
            cur_d -= 1;
        }
        cur[y][x] = c;
        if (cur[y][x] == '.') {
            cur_d += 1;
        }

        REP (is_hr, 2) {
            REP3 (offset, - len_max + 1, 0 + 1) {
                const trie_node *ptr = trie;
                REP (i, len_max) {
                    ptr = peek_trie(is_hr ? cur[y][modadd(modsub(x, - offset, N), i, N)] : cur[modadd(modsub(y, - offset, N), i, N)][x], ptr);
                    if (ptr == nullptr) {
                        break;
                    }
                    if (offset + i >= 0) {
                        for (int j : ptr->indices) {
                            if (not used[j]) {
                                cur_c += 1;
                            }
                            used[j] += 1;
                        }
                    }
                }
            }
        }
    };

    int64_t iteration = 0;
    double temperature = 1.0;
    for (; ; ++ iteration) {
        if (iteration % 64 == 0) {
            chrono::high_resolution_clock::time_point clock_now = chrono::high_resolution_clock::now();
            temperature = static_cast<long double>((clock_end - clock_now).count()) / (clock_end - clock_begin).count();
            if (temperature <= 0.0) {
                break;
            }
        }

        int64_t prv_score = calculate_score(m, cur_c, cur_d);

        int y = uniform_int_distribution<int>(0, N - 1)(gen);
        int x = uniform_int_distribution<int>(0, N - 1)(gen);
        std::function<void ()> accept = [&]() {};
        std::function<void ()> reject = [&]() {};

        double choice = uniform_real_distribution<double>()(gen);
        if (choice < 0.99 or cur_c == m) {
            char c = uniform_int_distribution<char>('A', 'H')(gen);
            char preserved = cur[y][x];
            update(y, x, c);
            reject = [&, preserved]() {
                update(y, x, preserved);
            };

        } else {
            int i;
            while (true) {
                i = uniform_int_distribution<int>(0, m - 1)(gen);
                if (not used[i]) {
                    break;
                }
            }
            bool is_hr = uniform_int_distribution<int>(0, 4 - 1)(gen);
            string preserved;
            REP (z, s[i].length()) {
                int ny = (y + (is_hr ? z : 0)) % N;
                int nx = (x + (is_hr ? 0 : z)) % N;
                preserved += cur[ny][nx];
                update(ny, nx, s[i][z]);
            }
            reject = [&, i, is_hr, preserved]() {
                REP (z, s[i].length()) {
                    int ny = (y + (is_hr ? z : 0)) % N;
                    int nx = (x + (is_hr ? 0 : z)) % N;
                    update(ny, nx, preserved[z]);
                }
            };
        }

        int64_t nxt_score = calculate_score(m, cur_c, cur_d);

        int64_t delta = nxt_score - prv_score;
        auto probability = [&]() {
            constexpr long double boltzmann = 0.00001;
            return exp(boltzmann * delta / temperature);
        };
        if (delta >= 0 or bernoulli_distribution(probability())(gen)) {
            accept();
            if (highscore < nxt_score) {
                highscore = nxt_score;
                ans = cur;
                ans_c = cur_c;
                ans_d = cur_d;
            }
        } else {
            reject();
        }
    }

    int sum_length = 0;
    for (auto &s_i : s) {
        sum_length += s_i.length();
    }
    int average_length = round(sum_length / m);

    cerr << "m = " << m << endl;
    cerr << "average length = " << average_length << endl;
    cerr << "iteration = " << iteration << endl;
    cerr << "c = " << ans_c << endl;
    cerr << "d = " << ans_d << endl;
    cerr << "score = " << highscore / 1.0e8 << "e8" << endl;
    return ans;
}

int main() {
    constexpr auto TIME_LIMIT = chrono::milliseconds(3000);
    chrono::high_resolution_clock::time_point clock_begin = chrono::high_resolution_clock::now();
    xor_shift_128 gen(20210425);

    int n, m; cin >> n >> m;
    assert (n == N);
    vector<string> s(m);
    REP (i, m) {
        cin >> s[i];
    }
    array<array<char, N>, N> ans = solve(m, s, gen, clock_begin + chrono::duration_cast<chrono::milliseconds>(TIME_LIMIT * 0.95));
    REP (y, N) {
        REP (x, N) {
            cout << ans[y][x];
        }
        cout << endl;
    }
    return 0;
}
