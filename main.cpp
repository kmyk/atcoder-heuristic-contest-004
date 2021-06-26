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
    REP (i, s.size()) {
        if (s[i] != f[y][(x + i) % N]) {
            return false;
        }
    }
    return true;
}

bool is_vertical_subarray_at(const string& s, int y, int x, const array<array<char, N>, N>& f) {
    REP (i, s.size()) {
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

template <class RandomEngine>
array<array<char, N>, N> solve(const int m, const vector<string>& s, RandomEngine& gen, chrono::high_resolution_clock::time_point clock_end) {
    chrono::high_resolution_clock::time_point clock_begin = chrono::high_resolution_clock::now();

    array<array<char, N>, N> ans = get_empty_board();
    pair<int, int> highscore = make_pair(- m, N * N);  // (m - c, d)

    unordered_map<string, vector<int>> occur;
    REP (i, m) {
        occur[s[i]].push_back(i);
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

        array<array<char, N>, N> cur = get_empty_board();
        int d = N * N;
        vector<bool> used(m);
        int y = 0;
        int x = 0;

        vector<int> order(m);
        iota(ALL(order), 0);
        shuffle(ALL(order), gen);
        for (int i : order) {
            if (used[i]) {
                continue;
            }
            int offset = -1;
            REP_R (k, (int)s[i].size()) {
                if (x - k >= 0 and x - k + s[i].size() <= N and is_horizontal_subarray_at(s[i].substr(0, k), y, x - k, cur)) {
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
            REP3 (z, offset, s[i].size()) {
                assert (x < N);
                cur[y][x] = s[i][z];

                REP3 (len, LEN_MIN, LEN_MAX + 1) {
                    if (x - len >= 0) {
                        string t = get_horizontal_subarray_at(y, x - len, len, cur);
                        auto it = occur.find(t);
                        if (it != occur.end()) {
                            for (int j : it->second) {
                                used[j] = true;
                            }
                        }
                    }
                    if (y - len >= 0) {
                        string t = get_vertical_subarray_at(y - len, x, len, cur);
                        auto it = occur.find(t);
                        if (it != occur.end()) {
                            for (int j : it->second) {
                                used[j] = true;
                            }
                        }
                    }
                }

                x += 1;
                d -= 1;
            }
        }

        int c = count(ALL(used), true);
        pair<int, int> score = make_pair(m - c, d);
        if (highscore < score) {
            highscore = score;
            ans = cur;
        }
    }

    cerr << "ans =" << endl;
    cerr << "m = " << m << endl;
    REP (y, N) {
        cerr << "    ";
        REP (x, N) {
            cerr << ans[y][x];
        }
        cerr << endl;
    }
    cerr << "score = (c, d) = (" << m - highscore.first << ", " << highscore.second << ")" << endl;
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
