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

template <class RandomEngine>
array<array<char, N>, N> solve(const int m, const vector<string>& s, RandomEngine& gen, chrono::high_resolution_clock::time_point clock_end) {
    chrono::high_resolution_clock::time_point clock_begin = chrono::high_resolution_clock::now();

    array<array<char, N>, N> ans = get_empty_board();
    pair<int, int> highscore = make_pair(- m, N * N);  // (m - c, d)

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
    REP (len, LEN_MAX + 1) {
        g_from[len].resize(m);
    }
    REP (i, m) {
        REP (j, m) {
            int dist = s[i].length() - max_common_length[i][j];
            assert (0 <= dist and dist < LEN_MAX + 1);
            g_from[dist][i].push_back(j);
        }
    }

    array<vector<vector<int>>, LEN_MAX + 1> g_to;
    REP (len, LEN_MAX + 1) {
        g_to[len].resize(m);
    }
    REP (i, m) {
        REP (j, m) {
            int dist = s[j].length() - max_common_length[i][j];
            assert (0 <= dist and dist < LEN_MAX + 1);
            g_to[dist][i].push_back(j);
        }
    }

    int64_t iteration = 0;
    double temperature = 1.0;
    for (; ; ++ iteration) {
        if (iteration % 64 == 0) {
            chrono::high_resolution_clock::time_point clock_now = chrono::high_resolution_clock::now();
            temperature = static_cast<long double>((clock_end - clock_now).count()) / (clock_end - clock_begin).count();
            if (temperature <= 0.0) {
                cerr << "done  (iteration = " << iteration << ")" << endl;
                break;
            }
        }

        array<array<char, N>, N> f = get_empty_board();
        int c = 0;
        int d = N * N;
        vector<bool> used(m);
        int y = 0;
        int x = 0;

        int last = -1;
        while (c < m) {
            int i = -1;
            if (last == -1) {
                i = uniform_int_distribution<int>(0, m - 1)(gen);
            } else {
                vector<int> cands;
                REP (len, LEN_MAX + 1) {
                    vector<int> cands;
                    for (int j : g_to[len][last]) {
                        if (not used[j]) {
                            cands.push_back(j);
                        }
                    }
                    if (not cands.empty()) {
                        i = cands[uniform_int_distribution<int>(0, (int)cands.size() - 1)(gen)];
                        break;
                    }
                }
            }
            assert (i != -1);

            int offset = -1;
            REP_R (k, (int)s[i].length()) {
                if (x - k >= 0 and x - k + s[i].length() <= N and is_horizontal_subarray_at(s[i].substr(0, k), y, x - k, f)) {
                    offset = k;
                    break;
                }
            }
            if (offset == -1) {
                ++ y;
                x = 0;
                if (y >= N) {
                    break;
                }
                offset = 0;
            }

            REP3 (z, offset, s[i].length()) {
                assert (x < N);
                f[y][x] = s[i][z];

                REP3 (len, LEN_MIN, LEN_MAX + 1) {
                    if (x - len >= 0) {
                        string t = get_horizontal_subarray_at(y, x - len, len, f);
                        auto it = occur.find(t);
                        if (it != occur.end()) {
                            for (int j : it->second) {
                                if (not used[j]) {
                                    c += 1;
                                    used[j] = true;
                                }
                            }
                        }
                    }
                    if (y - len >= 0) {
                        string t = get_vertical_subarray_at(y - len, x, len, f);
                        auto it = occur.find(t);
                        if (it != occur.end()) {
                            for (int j : it->second) {
                                if (not used[j]) {
                                    c += 1;
                                    used[j] = true;
                                }
                            }
                        }
                    }
                }

                x += 1;
                d -= 1;
            }

            // assert (used[i]);
            last = i;
        }

        pair<int, int> score = make_pair(m - c, d);
        if (highscore < score) {
            highscore = score;
            ans = f;
        }
    }

    int sum_length = 0;
    for (auto &s_i : s) {
        sum_length += s_i.length();
    }
    int average_length = round(sum_length / m);

    cerr << "m = " << m << endl;
    cerr << "average length = " << average_length << endl;
    cerr << "c = " << m - highscore.first << endl;
    cerr << "d = " << highscore.second << endl;
    cerr << "score = " << calculate_score(m, m - highscore.first, highscore.second) / 1.0e8 << "e8" << endl;
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
