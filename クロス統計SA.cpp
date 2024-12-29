#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iomanip> 
#include <unordered_map>  // ★ これを追加


// ---------------------------
//  グローバル定数 or サイズ (後でCSVから動的に判明する場合もある)
// ---------------------------
static const std::vector<int> ownership_types = {0, 1, 2, 3, 4};
static const std::vector<int> bldg_types = {0, 1, 2, 3};
static const std::vector<int> floor_types = {0, 1, 2, 3, 4, 5};

// ---------------------------
//  CSV から家族類型×所有の世帯数 M_{i,j} を読み込む
// ---------------------------
bool read_family_ownership_csv(
    const std::string& filename,
    std::vector<std::string>& family_types, // 動的に家族類型を追加
    std::vector<std::vector<double>>& M,    // M[i][j] が世帯数
    double& total_households                // 合計世帯数を格納
) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        std::cerr << "Error: Failed to open file: " << filename << std::endl;
        return false;
    }

    // ファイルを一行ずつ読む
    // CSV の列は「家族類型,所有,世帯数」を想定
    // ヘッダーがある場合は一行読み飛ばすなど調整してください
    std::string line;
    bool is_header = true;

    // まず family_types を動的に収集するために使う
    // family_type -> index マッピング
    std::vector<std::string> unique_family_types;
    std::unordered_map<std::string, int> family_type_to_index;

    total_households = 0.0;

    while (std::getline(ifs, line)) {
        // 必要に応じてヘッダーを読み飛ばす
        // ここでは最初の行がヘッダーと仮定
        if (is_header) {
            is_header = false;
            continue; 
        }

        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        std::string token;
        std::string family_type;
        int ownership = 0;
        double households = 0.0;

        // 列をパース
        // family_type
        std::getline(ss, token, ',');
        family_type = token;

        // ownership
        std::getline(ss, token, ',');
        ownership = std::stoi(token);

        // 世帯数
        std::getline(ss, token, ',');
        households = std::stod(token);

        // family_type がこれまでに無い場合は追加
        if (family_type_to_index.find(family_type) == family_type_to_index.end()) {
            int idx = static_cast<int>(unique_family_types.size());
            family_type_to_index[family_type] = idx;
            unique_family_types.push_back(family_type);
        }

        // 集計
        total_households += households;
    }

    // family_types をソートしてから元のインデックスを確定
    std::sort(unique_family_types.begin(), unique_family_types.end());
    family_types = unique_family_types;

    // family_type_to_index を再構築（ソート後に再マッピング）
    family_type_to_index.clear();
    for (int i = 0; i < (int)family_types.size(); i++) {
        family_type_to_index[family_types[i]] = i;
    }

    // M配列の確保
    int I = (int)family_types.size();
    int J = (int)ownership_types.size();
    M.assign(I, std::vector<double>(J, 0.0));

    // 再度ファイルを読み込み、今度は M に値を入れる
    ifs.clear();
    ifs.seekg(0, std::ios::beg);

    is_header = true;
    while (std::getline(ifs, line)) {
        if (is_header) {
            is_header = false;
            continue; 
        }
        if (line.empty()) {
            continue;
        }
        std::stringstream ss(line);
        std::string token;
        std::string family_type;
        int ownership = 0;
        double households = 0.0;

        std::getline(ss, token, ',');
        family_type = token;
        std::getline(ss, token, ',');
        ownership = std::stoi(token);
        std::getline(ss, token, ',');
        households = std::stod(token);

        int i_idx = family_type_to_index[family_type];
        int j_idx = std::distance(
            ownership_types.begin(),
            std::find(ownership_types.begin(), ownership_types.end(), ownership)
        );
        M[i_idx][j_idx] += households;
    }

    ifs.close();
    return true;
}

// ---------------------------
//  CSV から 所有×建て方×面積 の世帯数 T_{j,k,l} を読み込む
// ---------------------------
bool read_ownership_bldg_floor_csv(
    const std::string& filename,
    std::vector<std::vector<std::vector<double>>>& T
) {
    // T[j][k][l]
    // j: ownership, k: bldg, l: floor
    // 事前に次元を確保
    int J = (int)ownership_types.size();
    int K = (int)bldg_types.size();
    int L = (int)floor_types.size();
    T.assign(J, std::vector<std::vector<double>>(K, std::vector<double>(L, 0.0)));

    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        std::cerr << "Error: Failed to open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    bool is_header = true;
    while (std::getline(ifs, line)) {
        // ヘッダー読み飛ばし
        if (is_header) {
            is_header = false;
            continue;
        }
        if (line.empty()) {
            continue;
        }
        std::stringstream ss(line);
        std::string token;
        int ownership = 0;
        int bldg = 0;
        int floor = 0;
        double households = 0.0;

        // ownership
        std::getline(ss, token, ',');
        ownership = std::stoi(token);

        // 建て方
        std::getline(ss, token, ',');
        bldg = std::stoi(token);

        // 面積
        std::getline(ss, token, ',');
        floor = std::stoi(token);

        // 世帯数
        std::getline(ss, token, ',');
        households = std::stod(token);

        // j_idx, k_idx, l_idx
        int j_idx = std::distance(
            ownership_types.begin(),
            std::find(ownership_types.begin(), ownership_types.end(), ownership)
        );
        int k_idx = std::distance(
            bldg_types.begin(),
            std::find(bldg_types.begin(), bldg_types.end(), bldg)
        );
        int l_idx = std::distance(
            floor_types.begin(),
            std::find(floor_types.begin(), floor_types.end(), floor)
        );

        T[j_idx][k_idx][l_idx] += households;
    }

    ifs.close();
    return true;
}

// ---------------------------
//  エネルギー関数の定義
//  Python での実装と同様に、
//  E(X) = Σ_{j,k,l} |(Σ_i X_{i,j,k,l}) - T_{j,k,l}|
// ---------------------------
double calc_energy(
    const std::vector<std::vector<std::vector<std::vector<int>>>>& X,
    const std::vector<std::vector<std::vector<double>>>& T
) {
    // X[i][j][k][l]
    // T[j][k][l]
    double sum_diff = 0.0;

    int I = (int)X.size();
    int J = (int)X[0].size();
    int K = (int)X[0][0].size();
    int L = (int)X[0][0][0].size();

    for (int j = 0; j < J; j++) {
        for (int k = 0; k < K; k++) {
            for (int l = 0; l < L; l++) {
                double sum_over_i = 0.0;
                for (int i = 0; i < I; i++) {
                    sum_over_i += (double)X[i][j][k][l];
                }
                double diff = std::abs(sum_over_i - T[j][k][l]);
                sum_diff += diff;
            }
        }
    }
    return sum_diff;
}

// ---------------------------
//  近傍解の生成
// ---------------------------
std::vector<std::vector<std::vector<std::vector<int>>>> generate_neighbor(
    const std::vector<std::vector<std::vector<std::vector<int>>>>& X,
    std::mt19937& rng,
    int max_delta = 5
) {
    // X[i][j][k][l]
    // 近傍解を生成する (整数制約付き)
    auto X_new = X; // ディープコピー

    int I = (int)X.size();
    int J = (int)X[0].size();
    int K = (int)X[0][0].size();
    int L = (int)X[0][0][0].size();

    std::uniform_int_distribution<int> dist_i(0, I - 1);
    std::uniform_int_distribution<int> dist_j(0, J - 1);
    std::uniform_int_distribution<int> dist_k(0, K - 1);
    std::uniform_int_distribution<int> dist_l(0, L - 1);
    std::uniform_int_distribution<int> dist_delta(-max_delta, max_delta);

    int i_idx = dist_i(rng);
    int j_idx = dist_j(rng);
    int k1 = dist_k(rng);
    int l1 = dist_l(rng);
    int k2 = dist_k(rng);
    int l2 = dist_l(rng);

    while (k1 == k2 && l1 == l2) {
        k2 = dist_k(rng);
        l2 = dist_l(rng);
    }

    int delta = dist_delta(rng);

    // 更新
    int new_val_1 = X_new[i_idx][j_idx][k1][l1] + delta;
    int new_val_2 = X_new[i_idx][j_idx][k2][l2] - delta;

    // 非負制約
    if (new_val_1 < 0 || new_val_2 < 0) {
        // 破綻する場合は変更しないでそのまま返す
        return X;
    } else {
        X_new[i_idx][j_idx][k1][l1] = new_val_1;
        X_new[i_idx][j_idx][k2][l2] = new_val_2;
    }

    return X_new;
}

// ---------------------------
//  メイン
// ---------------------------
int main() {
    // 1. データ読み込みと初期化
    std::string file_family_ownership = "family_ownership.csv";
    std::string file_ownership_bldg_floor = "ownership_building_floor.csv";

    std::vector<std::string> family_types;              // i軸 (家族類型)
    std::vector<std::vector<double>> M;                 // M[i][j]: 家族類型×所有
    double N = 0.0;                                     // 世帯数の総和
    if (!read_family_ownership_csv(file_family_ownership, family_types, M, N)) {
        return 1;
    }
    int I = (int)family_types.size();
    int J = (int)ownership_types.size();

    std::vector<std::vector<std::vector<double>>> T;    // T[j][k][l]
    if (!read_ownership_bldg_floor_csv(file_ownership_bldg_floor, T)) {
        return 1;
    }
    int K = (int)bldg_types.size();
    int L = (int)floor_types.size();

    std::cout << "世帯数の総和: " << N << std::endl;

    // 2. 初期解 X_{i,j,k,l} の生成
    //    (「総数をランダムに保ちながら分配する」方式)
    //    Xは整数配列
    std::vector<std::vector<std::vector<std::vector<int>>>> X(
        I,
        std::vector<std::vector<std::vector<int>>>(
            J,
            std::vector<std::vector<int>>(
                K,
                std::vector<int>(L, 0)
            )
        )
    );

    std::random_device rd;
    std::mt19937 rng(rd());

    // ランダム分配
    for (int i_idx = 0; i_idx < I; i_idx++) {
        for (int j_idx = 0; j_idx < J; j_idx++) {
            int total_ij = (int)std::round(M[i_idx][j_idx]);
            if (total_ij <= 0) {
                continue;
            }
            // K*L 個の初期値をすべて 0 にして、1ずつ割り振る
            std::vector<int> rand_vec(K * L, 0);
            std::uniform_int_distribution<int> dist_index(0, K * L - 1);
            for (int p = 0; p < total_ij; p++) {
                int rand_idx = dist_index(rng);
                rand_vec[rand_idx] += 1;
            }
            // X に割り振り
            int count = 0;
            for (int k_idx = 0; k_idx < K; k_idx++) {
                for (int l_idx = 0; l_idx < L; l_idx++) {
                    X[i_idx][j_idx][k_idx][l_idx] = rand_vec[count];
                    count++;
                }
            }
        }
    }

    // 3. エネルギー関数 (既に定義済み)
    double current_energy = calc_energy(X, T);

    // ベスト解の保持
    auto best_state = X;
    double best_energy = current_energy;

    // 4. シミュレーテッド・アニーリングの設定
    double initial_temp = 0.5;
    double final_temp   = 0.1;
    // Python では max_iterations = N * 1000 としているが、N が大きいとC++で時間がかかるため調整してください
    long long max_iterations = static_cast<long long>(N * 1000.0);

    // 進捗確認用
    std::cout << "max_iterations = " << max_iterations << std::endl;

    // 時間計測開始
    auto start_time = std::chrono::steady_clock::now();

    // 5. シミュレーテッド・アニーリング
    // エネルギー履歴 (大きすぎる場合は格納しないほうが良い)
    std::vector<double> energy_history;
    energy_history.reserve(max_iterations);

    for (long long iteration = 0; iteration < max_iterations; iteration++) {
        // 温度設定 (指数冷却)
        double T_current = initial_temp * std::pow((final_temp / initial_temp), (double)iteration / (double)max_iterations);

        // 近傍解生成
        auto new_state = generate_neighbor(X, rng, 5);
        double new_energy = calc_energy(new_state, T);

        double dE = new_energy - current_energy;
        // メトロポリス基準
        if (dE < 0.0) {
            // 良化
            X = new_state;
            current_energy = new_energy;
            if (current_energy < best_energy) {
                best_energy = current_energy;
                best_state = X;
            }
        } else {
            // 悪化だが、確率で受容
            double prob = std::exp(-dE / T_current);
            std::uniform_real_distribution<double> dist_01(0.0, 1.0);
            double r = dist_01(rng);
            if (r < prob) {
                X = new_state;
                current_energy = new_energy;
                if (current_energy < best_energy) {
                    best_energy = current_energy;
                    best_state = X;
                }
            }
        }

        // 進捗状況 (一定間隔で表示) 
        if (iteration % (max_iterations / 100 + 1) == 0) {
            double ratio = (double)iteration / (double)max_iterations * 100.0;
            std::cout << "[" << std::fixed << std::setprecision(1) 
                    << ratio << "%] "
                    << "Iteration " << iteration << " / " << max_iterations
                    << "  Current E: " << current_energy
                    << "  Best E: " << best_energy << std::endl;
        }


        // エネルギー履歴
        energy_history.push_back(current_energy);
    }

    // 時間計測終了
    auto end_time = std::chrono::steady_clock::now();
    double elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();

    // 6. 結果の出力
    std::cout << "最終エネルギー: " << current_energy << std::endl;
    std::cout << "ベストエネルギー: " << best_energy << std::endl;
    std::cout << "計算にかかった時間: " << elapsed_time << "秒" << std::endl;

    // 結果をCSVに出力 (output_households.csv)
    {
        std::ofstream ofs("output_households.csv");
        if (!ofs.is_open()) {
            std::cerr << "Error: Cannot open output file." << std::endl;
            return 1;
        }
        // ヘッダー
        ofs << "家族類型,所有,建て方,延床面積,世帯数\n";

        for (int i = 0; i < I; i++) {
            for (int j = 0; j < J; j++) {
                for (int k = 0; k < K; k++) {
                    for (int l = 0; l < L; l++) {
                        ofs << family_types[i] << ","
                            << ownership_types[j] << ","
                            << bldg_types[k] << ","
                            << floor_types[l] << ","
                            << best_state[i][j][k][l]
                            << "\n";
                    }
                }
            }
        }
        ofs.close();
    }

    // 必要ならエネルギー履歴もファイルに保存
    // ここでは簡易的に出力する例のみ
    /*
    {
        std::ofstream ofs_e("energy_history.csv");
        ofs_e << "iteration,energy\n";
        for (size_t i = 0; i < energy_history.size(); i++) {
            ofs_e << i << "," << energy_history[i] << "\n";
        }
    }
    */

    return 0;
}
