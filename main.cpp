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
constexpr int LEN_MIN = 2;
constexpr int LEN_MAX = 12;

array<array<char, N>, N> get_empty_board() {
    array<array<char, N>, N> f;
    REP (y, N) {
        fill(ALL(f[y]), '.');
    }
    return f;
}

bool is_horizontal_subarray_at(const string& s, int y, int x, const array<array<char, N>, N>& f) {
    REP (i, s.length()) {
        if (s[i] != f[y][(x + i) % N]) {
            return false;
        }
    }
    return true;
}

bool is_vertical_subarray_at(const string& s, int y, int x, const array<array<char, N>, N>& f) {
    REP (i, s.length()) {
        if (s[i] != f[(y + i) % N][x]) {
            return false;
        }
    }
    return true;
}

string get_horizontal_subarray_at(int y, int x, int len, const array<array<char, N>, N>& f) {
    string s(len, '\0');
    REP (i, len) {
        s[i] = f[y][(x + i) % N];
    }
    return s;
}

string get_vertical_subarray_at(int y, int x, int len, const array<array<char, N>, N>& f) {
    string s(len, '\0');
    REP (i, len) {
        s[i] = f[(y + i) % N][x];
    }
    return s;
}

int64_t calculate_score(int m, int c, int d) {
    return c < m ? 1.0e8 * c / m : 1.0e8 * 2 * N * N / (2 * N * N - d);
}

void print_field(ostream &out, const array<array<char, N>, N>& f) {
    REP (y, N) {
        REP (x, N) {
            out << f[y][x];
        }
        out << endl;
    }
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

    unordered_map<string, vector<int>> occur;
    REP (i, m) {
        occur[s[i]].push_back(i);
    }

    vector<vector<int>> max_common_length(m, vector<int>(m));
    REP (i, m) {
        REP (j, m) {
            REP_R (k, min(s[i].length(), s[j].length()) + 1) {
                if (s[i].substr(s[i].length() - k) == s[j].substr(0, k)) {
                    max_common_length[i][j] = k;
                    break;
                }
            }
        }
    }

    array<vector<vector<int>>, LEN_MAX + 1> g_from;
    REP (len, len_max + 1) {
        g_from[len].resize(m);
    }
    REP (i, m) {
        REP (j, m) {
            int dist = s[i].length() - max_common_length[i][j];
            assert (0 <= dist and dist < len_max + 1);
            g_from[dist][i].push_back(j);
        }
    }

    array<vector<vector<int>>, LEN_MAX + 1> g_to;
    REP (len, len_max + 1) {
        g_to[len].resize(m);
    }
    REP (i, m) {
        REP (j, m) {
            int dist = s[j].length() - max_common_length[i][j];
            assert (0 <= dist and dist < len_max + 1);
            g_to[dist][i].push_back(j);
        }
    }

    array<array<char, N>, N> cur = get_empty_board();
    int cur_c = 0;
    int cur_d = N * N;
    vector<int> used(m);

    auto update = [&](int y, int x, char c) {
        assert (0 <= y and y < N);
        assert (0 <= x and x < N);

        REP (is_hr, 2) {
            REP3 (len, len_min, len_max + 1) {
                REP3 (delta, - len + 1, 0 + 1) {
                    string t = is_hr
                        ? get_horizontal_subarray_at(y, (x + delta + N) % N, len, cur)
                        : get_vertical_subarray_at((y + delta + N) % N, x, len, cur);
                    auto it = occur.find(t);
                    if (it != occur.end()) {
                        for (int j : it->second) {
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
            REP3 (len, LEN_MIN, LEN_MAX + 1) {
                REP3 (delta, - len + 1, 0 + 1) {
                    string t = is_hr
                        ? get_horizontal_subarray_at(y, (x + delta + N) % N, len, cur)
                        : get_vertical_subarray_at((y + delta + N) % N, x, len, cur);
                    auto it = occur.find(t);
                    if (it != occur.end()) {
                        for (int j : it->second) {
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
        if (choice < 0.8) {
            char c = uniform_int_distribution<char>('A', 'H')(gen);
            char preserved = cur[y][x];
            update(y, x, c);
            reject = [&, preserved]() {
                update(y, x, preserved);
            };

        } else if (choice < 1.0 and cur_c < m) {
            int i;
            while (true) {
                i = uniform_int_distribution<int>(0, m - 1)(gen);
                if (not used[i]) {
                    break;
                }
            }
            int is_hr = uniform_int_distribution<int>(0, 4 - 1)(gen);
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

        } else {
            bool is_hr = bernoulli_distribution(0.5)(gen);
            if (is_hr) {
                char c = cur[y][0];
                REP (z, N - 1) {
                    update(y, z, cur[y][z + 1]);
                }
                update(y, N - 1, c);
            } else {
                char c = cur[0][x];
                REP (z, N - 1) {
                    update(z, x, cur[z + 1][x]);
                }
                update(N - 1, x, c);
            }
            reject = [&, is_hr]() {
                if (is_hr) {
                    char c = cur[y][N - 1];
                    REP_R (z, N - 1) {
                        update(y, z + 1, cur[y][z]);
                    }
                    update(y, 0, c);
                } else {
                    char c = cur[N - 1][x];
                    REP_R (z, N - 1) {
                        update(z + 1, x, cur[z][x]);
                    }
                    update(0, x, c);
                }
            };
        }

        int64_t nxt_score = calculate_score(m, cur_c, cur_d);

        int64_t delta = nxt_score - prv_score;
        auto probability = [&]() {
            constexpr long double boltzmann = 0.00002;
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
    print_field(cout, ans);
    return 0;
}
