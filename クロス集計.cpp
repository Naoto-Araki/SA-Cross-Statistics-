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
#include <unordered_map>

// ------------------------------------------------------------
// グローバル定数 (今回のサンプル。実際には動的に作るなど適宜対応)
// ------------------------------------------------------------
static const std::vector<int> ownership_types = {0, 1, 2, 3, 4};
static const std::vector<int> bldg_types = {0, 1, 2, 3};
static const std::vector<int> floor_types = {0, 1, 2, 3, 4, 5};

// ------------------------------------------------------------
// 1) family_ownership.csv の読み込み
//    列: 「家族類型, 所有, 世帯数, 世帯人員数」
// ------------------------------------------------------------
bool read_family_ownership_csv(
    const std::string& filename,
    std::vector<std::string>& family_types,           // [i]  家族類型の名前 (ソート済み)
    std::vector<std::vector<double>>& M_households,   // [i][j] 家族類型×所有 の世帯数
    std::vector<std::vector<double>>& M_persons,      // [i][j] 家族類型×所有 の人数
    double& total_households,                         // 全体の世帯数合計
    double& total_persons                             // 全体の人数合計
) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        std::cerr << "Error: Failed to open file: " << filename << std::endl;
        return false;
    }

    // まず family_types を動的に収集する
    std::unordered_map<std::string, int> family_type_map;
    std::vector<std::string> unique_family_types;

    total_households = 0.0;
    total_persons    = 0.0;

    bool is_header = true;
    std::string line;
    // 第1回読み込み: 家族類型のユニーク取得と合計の算出
    while (std::getline(ifs, line)) {
        if (is_header) {
            // ヘッダ行を読み飛ばす(想定)
            is_header = false;
            continue;
        }
        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        std::string family_type_str;
        int ownership_val;
        double hh, pp; // 世帯数, 人数

        // 列: 家族類型, 所有, 世帯数, 世帯人員数
        // 例: "A,0,1000,2500"
        std::getline(ss, family_type_str, ',');
        ss >> ownership_val;  // 所有
        ss.ignore(1, ',');    // カンマを飛ばす

        ss >> hh;             // 世帯数
        ss.ignore(1, ',');
        ss >> pp;             // 人数

        // 家族類型を登録
        if (family_type_map.find(family_type_str) == family_type_map.end()) {
            int idx = (int)unique_family_types.size();
            family_type_map[family_type_str] = idx;
            unique_family_types.push_back(family_type_str);
        }

        total_households += hh;
        total_persons    += pp;
    }

    // ソート
    std::sort(unique_family_types.begin(), unique_family_types.end());
    family_types = unique_family_types;

    // family_type_map を再構築 (ソート後)
    family_type_map.clear();
    for (int i = 0; i < (int)family_types.size(); i++) {
        family_type_map[family_types[i]] = i;
    }

    int I = (int)family_types.size();
    int J = (int)ownership_types.size();

    // M_households, M_persons のサイズを確保
    M_households.assign(I, std::vector<double>(J, 0.0));
    M_persons.assign(I,    std::vector<double>(J, 0.0));

    // ファイルを再度先頭へ
    ifs.clear();
    ifs.seekg(0, std::ios::beg);
    is_header = true;

    // 第2回読み込み: 具体的に M_households, M_persons を埋める
    while (std::getline(ifs, line)) {
        if (is_header) {
            is_header = false;
            continue;
        }
        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        std::string family_type_str;
        int ownership_val;
        double hh, pp;
        std::getline(ss, family_type_str, ',');
        ss >> ownership_val;
        ss.ignore(1, ',');
        ss >> hh;
        ss.ignore(1, ',');
        ss >> pp;

        int i_idx = family_type_map[family_type_str];
        // 所有のインデックス j
        int j_idx = (int)std::distance(
            ownership_types.begin(),
            std::find(ownership_types.begin(), ownership_types.end(), ownership_val)
        );

        M_households[i_idx][j_idx] += hh;
        M_persons[i_idx][j_idx]    += pp;
    }

    ifs.close();
    return true;
}

// ------------------------------------------------------------
// 2) ownership_bldg_floor.csv の読み込み
//    列: 「所有, 建て方, 面積, 世帯数, 世帯人員数」
// ------------------------------------------------------------
bool read_ownership_bldg_floor_csv(
    const std::string& filename,
    std::vector<std::vector<std::vector<double>>>& T_households,  // [j][k][l]
    std::vector<std::vector<std::vector<double>>>& T_persons      // [j][k][l]
) {
    int J = (int)ownership_types.size();
    int K = (int)bldg_types.size();
    int L = (int)floor_types.size();

    T_households.assign(J, std::vector<std::vector<double>>(K, std::vector<double>(L, 0.0)));
    T_persons.assign(J,    std::vector<std::vector<double>>(K, std::vector<double>(L, 0.0)));

    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        std::cerr << "Error: Failed to open file: " << filename << std::endl;
        return false;
    }

    bool is_header = true;
    std::string line;
    while (std::getline(ifs, line)) {
        if (is_header) {
            // ヘッダ行を読み飛ばす(想定)
            is_header = false;
            continue;
        }
        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        int own_val, bldg_val, floor_val;
        double hh, pp; // 世帯数, 人数

        // 列: 所有, 建て方, 面積, 世帯数, 世帯人員数
        // 例: "0,0,0,5000,12000"
        ss >> own_val;
        ss.ignore(1, ',');
        ss >> bldg_val;
        ss.ignore(1, ',');
        ss >> floor_val;
        ss.ignore(1, ',');
        ss >> hh;
        ss.ignore(1, ',');
        ss >> pp;

        int j_idx = (int)std::distance(
            ownership_types.begin(),
            std::find(ownership_types.begin(), ownership_types.end(), own_val)
        );
        int k_idx = (int)std::distance(
            bldg_types.begin(),
            std::find(bldg_types.begin(), bldg_types.end(), bldg_val)
        );
        int l_idx = (int)std::distance(
            floor_types.begin(),
            std::find(floor_types.begin(), floor_types.end(), floor_val)
        );

        T_households[j_idx][k_idx][l_idx] += hh;
        T_persons[j_idx][k_idx][l_idx]    += pp;
    }

    ifs.close();
    return true;
}

// ------------------------------------------------------------
// 3) エネルギー関数
//    E(X,Y) = Σ_{j,k,l}|Σ_i X[i,j,k,l] - T_households[j,k,l]|
//           + Σ_{j,k,l}|Σ_i Y[i,j,k,l] - T_persons[j,k,l]|
//    (必要に応じて重みを掛けるなど調整可能)
// ------------------------------------------------------------
double calc_energy(
    const std::vector<std::vector<std::vector<std::vector<int>>>>& X,  // 世帯数
    const std::vector<std::vector<std::vector<std::vector<int>>>>& Y,  // 人数
    const std::vector<std::vector<std::vector<double>>>& T_households,
    const std::vector<std::vector<std::vector<double>>>& T_persons,
    double alpha = 1.0,  // 世帯数誤差の重み
    double beta  = 1.0   // 人数誤差の重み
) {
    int I = (int)X.size();           // 家族類型の数
    int J = (int)X[0].size();        // 所有の数
    int K = (int)X[0][0].size();     // 建て方の数
    int L = (int)X[0][0][0].size();  // 延床面積の数

    double sum_diff_hh = 0.0;  // 世帯数
    double sum_diff_pp = 0.0;  // 人数

    for (int j = 0; j < J; j++) {
        for (int k = 0; k < K; k++) {
            for (int l = 0; l < L; l++) {
                double sumX = 0.0; // Xの合計(世帯数)
                double sumY = 0.0; // Yの合計(人数)
                for (int i = 0; i < I; i++) {
                    sumX += (double)X[i][j][k][l];
                    sumY += (double)Y[i][j][k][l];
                }
                double diff_hh = std::abs(sumX - T_households[j][k][l]);
                double diff_pp = std::abs(sumY - T_persons[j][k][l]);

                sum_diff_hh += diff_hh;
                sum_diff_pp += diff_pp;
            }
        }
    }
    return alpha * sum_diff_hh + beta * sum_diff_pp;
}

// ------------------------------------------------------------
// 4) 近傍解生成 (X, Y のどちらかを操作する例)
//    ・ 50%の確率で X を、50%の確率で Y を動かす
// ------------------------------------------------------------
std::pair<
    std::vector<std::vector<std::vector<std::vector<int>>>>,
    std::vector<std::vector<std::vector<std::vector<int>>>>
> generate_neighbor(
    const std::vector<std::vector<std::vector<std::vector<int>>>>& X,
    const std::vector<std::vector<std::vector<std::vector<int>>>>& Y,
    std::mt19937& rng,
    int max_delta = 5
) {
    auto X_new = X;
    auto Y_new = Y;

    int I = (int)X.size();
    int J = (int)X[0].size();
    int K = (int)X[0][0].size();
    int L = (int)X[0][0][0].size();

    // XかYかをランダムに
    std::uniform_int_distribution<int> dist_xy(0, 1);
    int which = dist_xy(rng);  // 0 => X, 1 => Y

    // (i,j) を固定して (k1,l1), (k2,l2) を選んで移動するサンプル
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

    if (which == 0) {
        // Xを操作
        int old_v1 = X_new[i_idx][j_idx][k1][l1];
        int old_v2 = X_new[i_idx][j_idx][k2][l2];
        int new_v1 = old_v1 + delta;
        int new_v2 = old_v2 - delta;

        // 非負チェック
        if (new_v1 < 0 || new_v2 < 0) {
            // 破綻なら元を返す
            return {X, Y};
        }
        X_new[i_idx][j_idx][k1][l1] = new_v1;
        X_new[i_idx][j_idx][k2][l2] = new_v2;
    } else {
        // Yを操作
        int old_v1 = Y_new[i_idx][j_idx][k1][l1];
        int old_v2 = Y_new[i_idx][j_idx][k2][l2];
        int new_v1 = old_v1 + delta;
        int new_v2 = old_v2 - delta;

        if (new_v1 < 0 || new_v2 < 0) {
            return {X, Y};
        }
        Y_new[i_idx][j_idx][k1][l1] = new_v1;
        Y_new[i_idx][j_idx][k2][l2] = new_v2;
    }

    return {X_new, Y_new};
}

// ------------------------------------------------------------
// メイン
// ------------------------------------------------------------
int main() {
    // 1) CSV 読み込み
    std::string file_family_ownership = "family_ownership.csv";
    std::string file_ownership_bldg_floor = "ownership_building_floor.csv";

    std::vector<std::string> family_types; 
    // M_households[i][j], M_persons[i][j]
    std::vector<std::vector<double>> M_households, M_persons;
    double total_households = 0.0;
    double total_persons    = 0.0;

    if (!read_family_ownership_csv(
        file_family_ownership,
        family_types,
        M_households,
        M_persons,
        total_households,
        total_persons
    )) {
        std::cerr << "Error: read_family_ownership_csv failed." << std::endl;
        return 1;
    }

    int I = (int)family_types.size();             // 家族類型の数
    int J = (int)ownership_types.size();          // 所有の数

    // T_households[j][k][l], T_persons[j][k][l]
    std::vector<std::vector<std::vector<double>>> T_households, T_persons;

    if (!read_ownership_bldg_floor_csv(
        file_ownership_bldg_floor,
        T_households,
        T_persons
    )) {
        std::cerr << "Error: read_ownership_bldg_floor_csv failed." << std::endl;
        return 1;
    }

    int K = (int)bldg_types.size(); // 建て方
    int L = (int)floor_types.size(); // 面積

    std::cout << "family_ownership.csv 読込完了 (I=" << I << ", J=" << J << ")\n"
              << "  総世帯数=" << total_households
              << ", 総人数="   << total_persons << std::endl;
    std::cout << "ownership_bldg_floor.csv 読込完了 (K=" << K << ", L=" << L << ")\n";

    // 2) X[i][j][k][l], Y[i][j][k][l] を初期生成
    //   - (i,j)単位の世帯数, 人数を (k,l)にランダム分配
    //   - X は世帯数, Y は人数
    std::vector<std::vector<std::vector<std::vector<int>>>> X(
        I, std::vector<std::vector<std::vector<int>>>(
            J, std::vector<std::vector<int>>(
               K, std::vector<int>(L, 0)
        ))
    );
    std::vector<std::vector<std::vector<std::vector<int>>>> Y(
        I, std::vector<std::vector<std::vector<int>>>(
            J, std::vector<std::vector<int>>(
               K, std::vector<int>(L, 0)
        ))
    );

    std::random_device rd;
    std::mt19937 rng(rd());

    // 家族類型 i, 所有 j の合計を (k,l) に乱数割り振り
    for(int i_idx=0; i_idx<I; i_idx++){
        for(int j_idx=0; j_idx<J; j_idx++){
            int hh = (int)std::round(M_households[i_idx][j_idx]); // 世帯数
            int pp = (int)std::round(M_persons[i_idx][j_idx]);    // 人数

            // 世帯数 X
            if(hh > 0){
                std::vector<int> tmp(K*L, 0);
                std::uniform_int_distribution<int> dist_kl(0, K*L - 1);
                for(int n=0; n<hh; n++){
                    int idx = dist_kl(rng);
                    tmp[idx]++;
                }
                // Xに書き込み
                int c=0;
                for(int k_idx=0; k_idx<K; k_idx++){
                    for(int l_idx=0; l_idx<L; l_idx++){
                        X[i_idx][j_idx][k_idx][l_idx] = tmp[c];
                        c++;
                    }
                }
            }

            // 人数 Y
            if(pp > 0){
                std::vector<int> tmp(K*L, 0);
                std::uniform_int_distribution<int> dist_kl(0, K*L - 1);
                for(int n=0; n<pp; n++){
                    int idx = dist_kl(rng);
                    tmp[idx]++;
                }
                // Yに書き込み
                int c=0;
                for(int k_idx=0; k_idx<K; k_idx++){
                    for(int l_idx=0; l_idx<L; l_idx++){
                        Y[i_idx][j_idx][k_idx][l_idx] = tmp[c];
                        c++;
                    }
                }
            }
        }
    }

    // 3) エネルギー関数
    double alpha = 1.0; // 世帯数誤差の重み
    double beta  = 1.0; // 人数誤差の重み

    double current_energy = calc_energy(X, Y, T_households, T_persons, alpha, beta);
    auto best_X = X;
    auto best_Y = Y;
    double best_energy = current_energy;

    // 4) SAパラメータ
    //   - とりあえず max_iterations = total_households * 500 とかに設定
    //   - 大きすぎると時間がかかるので要調整
    long long max_iterations = (long long)(total_households * 1000);

    double initial_temp = 1.0;
    double final_temp   = 0.1;

    std::cout << "初期エネルギー = " << current_energy << "\n"
              << "max_iterations = " << max_iterations << std::endl;

    auto start_time = std::chrono::steady_clock::now();

    // 5) SAループ
    for(long long iter = 0; iter < max_iterations; iter++){
        // 温度 (指数冷却)
        double ratio = (double)iter / (double)max_iterations;
        double temp_now = initial_temp * std::pow(final_temp / initial_temp, ratio);

        // 近傍解生成
        auto [X_new, Y_new] = generate_neighbor(X, Y, rng, 5);

        // 新エネルギー
        double new_energy = calc_energy(X_new, Y_new, T_households, T_persons, alpha, beta);

        double dE = new_energy - current_energy;
        if(dE < 0.0){
            // 改善
            X = X_new;
            Y = Y_new;
            current_energy = new_energy;
            if(current_energy < best_energy){
                best_energy = current_energy;
                best_X = X;
                best_Y = Y;
            }
        } else {
            // 悪化
            double prob = std::exp(-dE / temp_now);
            std::uniform_real_distribution<double> dist01(0.0, 1.0);
            double r = dist01(rng);
            if(r < prob){
                // 受容
                X = X_new;
                Y = Y_new;
                current_energy = new_energy;
                if(current_energy < best_energy){
                    best_energy = current_energy;
                    best_X = X;
                    best_Y = Y;
                }
            }
        }

        // 進捗表示 (例: 1%刻み)
        if(iter % (max_iterations/100 + 1) == 0){
            double progress = 100.0 * (double)iter / (double)max_iterations;
            std::cout << "[" << std::fixed << std::setprecision(1)
                      << progress << "%] "
                      << "Iter=" << iter << " / " << max_iterations
                      << "  E=" << current_energy
                      << "  best=" << best_energy
                      << std::endl;
        }
    }

    auto end_time = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    // ループ終了後

    std::cout << "最終エネルギー: " << current_energy << std::endl;
    std::cout << "ベストエネルギー: " << best_energy << std::endl;
    std::cout << "計算時間: " << elapsed << " 秒" << std::endl;

    // 6) 結果の出力 (例)
    {
        // 世帯数 (X)
        std::ofstream ofsX("output_households_クロス集計.csv");
        ofsX << "家族類型,所有,建て方,面積,世帯数\n";
        for(int i=0; i<I; i++){
            for(int j=0; j<J; j++){
                for(int k=0; k<K; k++){
                    for(int l=0; l<L; l++){
                        ofsX << family_types[i]       << ","
                             << ownership_types[j]    << ","
                             << bldg_types[k]         << ","
                             << floor_types[l]        << ","
                             << best_X[i][j][k][l]
                             << "\n";
                    }
                }
            }
        }
        ofsX.close();
    }
    {
        // 人数 (Y)
        std::ofstream ofsY("output_persons_クロス集計.csv");
        ofsY << "家族類型,所有,建て方,面積,人数\n";
        for(int i=0; i<I; i++){
            for(int j=0; j<J; j++){
                for(int k=0; k<K; k++){
                    for(int l=0; l<L; l++){
                        ofsY << family_types[i]       << ","
                             << ownership_types[j]    << ","
                             << bldg_types[k]         << ","
                             << floor_types[l]        << ","
                             << best_Y[i][j][k][l]
                             << "\n";
                    }
                }
            }
        }
        ofsY.close();
    }

    return 0;
}
