import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import time
from tqdm import tqdm 

# ---------------------------
# 1. データ読み込みと初期化
# ---------------------------

# --- 家族類型×所有 の合計値 M_{i,j} ---
df_family_ownership = pd.read_csv("family_ownership.csv")

# 家族類型と所有の軸
family_types = sorted(df_family_ownership["家族類型"].unique())
ownership_types = [0, 1, 2, 3, 4]

# 軸サイズ
I = len(family_types)   # 家族類型の数
J = len(ownership_types)  # 所有区分の数

# M_{i,j} の2次元配列を作成
M = np.zeros((I, J), dtype=float)
for row in df_family_ownership.itertuples():
    i_idx = family_types.index(row.家族類型)
    j_idx = ownership_types.index(row.所有)
    M[i_idx, j_idx] = row.世帯数

# 世帯数の総和を計算
N = df_family_ownership["世帯数"].sum()
print("世帯数の総和:", N)

# --- 所有×建て方×延床面積 の分布 T_{j,k,l} ---
df_ownership_bldg_floor = pd.read_csv("ownership_building_floor.csv")

# 建て方と延床面積の軸
bldg_types = [0, 1, 2, 3]
floor_types = [0, 1, 2, 3, 4, 5]

K = len(bldg_types)  # 建て方の種類
L = len(floor_types)  # 延床面積の種類

# T_{j,k,l} の3次元配列を作成
T = np.zeros((J, K, L), dtype=float)
for row in df_ownership_bldg_floor.itertuples():
    j_idx = ownership_types.index(row.所有)
    k_idx = bldg_types.index(row.建て方)
    l_idx = floor_types.index(row.面積)
    T[j_idx, k_idx, l_idx] = row.世帯数

# ---------------------------
# 2. 初期解 X_{i,j,k,l} の生成
# ---------------------------
X = np.zeros((I, J, K, L), dtype=float)

# スケールを用いた分配
# for i_idx in range(I):
#     for j_idx in range(J):
#         total_ij = int(M[i_idx, j_idx])  # 世帯数を整数にキャスト
#         if total_ij <= 0:
#             continue

#         # ランダム整数の初期分布を生成
#         rand_vec = np.random.randint(1, 101, K * L)  # 1〜100のランダム整数
#         rand_vec = (rand_vec / rand_vec.sum() * total_ij).astype(int)  # 割合を考慮して整数化

#         # 合計値を調整して誤差をゼロに
#         diff = total_ij - rand_vec.sum()  # 総和と元の値との差分を計算
#         if diff != 0:
#             rand_vec[0] += diff  # 差分を最初の要素に加算

#         # 配列に割り振る
#         count = 0
#         for k_idx in range(K):
#             for l_idx in range(L):
#                 X[i_idx, j_idx, k_idx, l_idx] = rand_vec[count]
#                 count += 1


# 総数をランダムに保ちながら分配
for i_idx in range(I):
    for j_idx in range(J):
        total_ij = int(M[i_idx, j_idx])  # 世帯数を整数にキャスト
        if total_ij <= 0:
            continue

        # K * L 個の初期値をすべてゼロにする
        rand_vec = np.zeros(K * L, dtype=int)

        # 総数をランダムに分配する
        for _ in range(total_ij):
            rand_idx = np.random.randint(0, K * L)  # ランダムにセルを選択
            rand_vec[rand_idx] += 1

        # 配列に割り振る
        count = 0
        for k_idx in range(K):
            for l_idx in range(L):
                X[i_idx, j_idx, k_idx, l_idx] = rand_vec[count]
                count += 1


# ---------------------------
# 3. エネルギー関数の定義
# ---------------------------
# def calc_energy(X, T):
#     """
#     エネルギー関数: E(X) = Σ_{j,k,l} (Σ_i X_{i,j,k,l} - T_{j,k,l})^2
#     """
#     diff = np.zeros((J, K, L), dtype=float)
#     for j in range(J):
#         for k in range(K):
#             for l in range(L):
#                 diff[j, k, l] = X[:, j, k, l].sum() - T[j, k, l]
#     return np.sum(diff**2)

def calc_energy(X, T):
    """
    エネルギー関数: E(X) = Σ_{j,k,l} (Σ_i X_{i,j,k,l} - T_{j,k,l})
    """
    diff = np.zeros((J, K, L), dtype=float)
    for j in range(J):
        for k in range(K):
            for l in range(L):
                diff[j, k, l] = np.abs(X[:, j, k, l].sum() - T[j, k, l])
    return np.sum(diff)

# ---------------------------
# 4. 近傍解の生成
# ---------------------------
def generate_neighbor(X, max_delta=5):
    """
    近傍解を生成する (整数制約付き)。
    max_delta: 入れ替え可能な最大人数。
    """
    X_new = X.copy()

    # ランダムに (i, j) を選ぶ
    i_idx = random.randint(0, I - 1)
    j_idx = random.randint(0, J - 1)

    # ランダムに (k1, l1) と (k2, l2) を選ぶ
    k1, l1 = random.randint(0, K - 1), random.randint(0, L - 1)
    k2, l2 = random.randint(0, K - 1), random.randint(0, L - 1)
    while (k1 == k2 and l1 == l2):
        k2, l2 = random.randint(0, K - 1), random.randint(0, L - 1)

    # ランダムな整数の変化量 (最小単位は1)
    delta = random.randint(-max_delta, max_delta)

    # 更新
    X_new[i_idx, j_idx, k1, l1] += delta
    X_new[i_idx, j_idx, k2, l2] -= delta

    # 非負制約をチェック (負になる場合は元の解を返す)
    if X_new[i_idx, j_idx, k1, l1] < 0 or X_new[i_idx, j_idx, k2, l2] < 0:
        return X

    return X_new

# ---------------------------
# 5. シミュレーテッド・アニーリングの実装
# ---------------------------
# 修正版 generate_neighbor を利用したシミュレーテッド・アニーリング
initial_temp = 0.5
final_temp = 0.1
max_iterations = N * 1000

current_state = X.copy()
current_energy = calc_energy(current_state, T)
best_state = current_state.copy()
best_energy = current_energy

energy_history = []
temperature_history = []

# 時間計測の開始
start_time = time.time()

from tqdm import tqdm  # 進捗バーを表示するライブラリ

# tqdmを使った進捗バー
for iteration in tqdm(range(max_iterations), desc="進捗状況"):
    T_current = initial_temp * (final_temp / initial_temp) ** (iteration / max_iterations)
    temperature_history.append(T_current)

    # 近傍解を生成 (整数制約付き)
    new_state = generate_neighbor(current_state, max_delta=5)
    new_energy = calc_energy(new_state, T)

    # メトロポリス基準
    dE = new_energy - current_energy
    if dE < 0 or random.random() < np.exp(-dE / T_current):
        current_state = new_state
        current_energy = new_energy
        if current_energy < best_energy:
            best_energy = current_energy
            best_state = current_state.copy()

    energy_history.append(current_energy)


# 結果をCSVに保存
df_output = pd.DataFrame([
    [family_types[i], ownership_types[j], bldg_types[k], floor_types[l], X[i, j, k, l]]
    for i in range(I) for j in range(J) for k in range(K) for l in range(L)
], columns=["家族類型", "所有", "建て方", "延床面積", "世帯数"])
df_output.to_csv("output_households.csv", index=False)

# 時間計測の終了
end_time = time.time()

# 総計算時間
elapsed_time = end_time - start_time

# ---------------------------
# 6. 結果の表示
# ---------------------------
print("最終エネルギー:", current_energy)
print("ベストエネルギー:", best_energy)
print(f"計算にかかった時間: {elapsed_time:.2f}秒")

# エネルギー推移をプロット
plt.figure(figsize=(10, 5))
plt.plot(energy_history, label="Energy")
plt.xlabel("Iteration")
plt.ylabel("Energy")
plt.title("Energy Convergence")
plt.legend()
plt.grid(True)
plt.show()

# 温度推移をプロット
plt.figure(figsize=(10, 5))
plt.plot(temperature_history, label="Temperature")
plt.xlabel("Iteration")
plt.ylabel("Temperature")
plt.title("Temperature Decay (Exponential Cooling)")
plt.legend()
plt.grid(True)
plt.show()