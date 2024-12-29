import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import time
from tqdm import tqdm 

# 必要な統計データをCSVから読み込む
df_family_ownership = pd.read_csv("family_ownership.csv")  # 家族類型×所有
df_ownership_building = pd.read_csv("ownership_building.csv")  # 所有×建て方
df_building_floor = pd.read_csv("building_floor.csv")  # 建て方×延床面積

# 家族類型、所有、建て方、延床面積の軸
family_types = sorted(df_family_ownership["家族類型"].unique())
ownership_types = sorted(df_family_ownership["所有"].unique())
building_types = sorted(df_ownership_building["建て方"].unique())
floor_types = sorted(df_building_floor["延床面積"].unique())

# 軸サイズ
I = len(family_types)
J = len(ownership_types)
K = len(building_types)
L = len(floor_types)

# 統計データを配列に変換
M = np.zeros((I, J), dtype=float)  # 家族類型×所有
T = np.zeros((J, K), dtype=float)  # 所有×建て方
U = np.zeros((K, L), dtype=float)  # 建て方×延床面積

# データの整形
for row in df_family_ownership.itertuples():
    i = family_types.index(row.家族類型)
    j = ownership_types.index(row.所有)
    M[i, j] = row.世帯数

for row in df_ownership_building.itertuples():
    j = ownership_types.index(row.所有)
    k = building_types.index(row.建て方)
    T[j, k] = row.世帯数

for row in df_building_floor.itertuples():
    k = building_types.index(row.建て方)
    l = floor_types.index(row.延床面積)
    U[k, l] = row.世帯数

# 初期値の設定 (M に基づきランダムに割り当て)
X = np.zeros((I, J, K, L), dtype=float)
for i in range(I):
    for j in range(J):
        if M[i, j] > 0:
            rand_vec = np.random.random(K * L)
            rand_vec /= rand_vec.sum()
            count = 0
            for k in range(K):
                for l in range(L):
                    X[i, j, k, l] = M[i, j] * rand_vec[count]
                    count += 1

def calc_energy(X, T, U, lambda_household=1.0, lambda_floor=1.0):
    """
    エネルギー関数:
    E(X) = λ_household * Σ_{j,k} (Σ_{i,l} X_{i,j,k,l} - T_{j,k})^2
         + λ_floor * Σ_{k,l} (Σ_{i,j} X_{i,j,k,l} - U_{k,l})^2
    """
    # 所有×建て方の誤差
    diff_household = np.zeros((J, K), dtype=float)
    for j in range(J):
        for k in range(K):
            diff_household[j, k] = np.sum(X[:, j, k, :]) - T[j, k]
    energy_household = np.sum(diff_household**2)

    # 建て方×延床面積の誤差
    diff_floor = np.zeros((K, L), dtype=float)
    for k in range(K):
        for l in range(L):
            diff_floor[k, l] = np.sum(X[:, :, k, l]) - U[k, l]
    energy_floor = np.sum(diff_floor**2)

    return lambda_household * energy_household + lambda_floor * energy_floor

def generate_neighbor(X, max_delta=5):
    """
    近傍解を生成する。
    """
    X_new = X.copy()

    # ランダムに (i, j, k, l) を選択
    i = random.randint(0, I - 1)
    j = random.randint(0, J - 1)
    k1, l1 = random.randint(0, K - 1), random.randint(0, L - 1)
    k2, l2 = random.randint(0, K - 1), random.randint(0, L - 1)
    while (k1, l1) == (k2, l2):
        k2, l2 = random.randint(0, K - 1), random.randint(0, L - 1)

    # ランダムな整数の変化量
    delta = random.randint(-max_delta, max_delta)

    # 再分配
    if X_new[i, j, k1, l1] + delta >= 0 and X_new[i, j, k2, l2] - delta >= 0:
        X_new[i, j, k1, l1] += delta
        X_new[i, j, k2, l2] -= delta

    return X_new

initial_temp = 0.5
final_temp = 0.1
max_iterations = 1000

current_state = X.copy()
current_energy = calc_energy(current_state, T, U)
best_state = current_state.copy()
best_energy = current_energy

for iteration in range(max_iterations):
    T_current = initial_temp * (final_temp / initial_temp) ** (iteration / max_iterations)

    # 近傍解生成
    new_state = generate_neighbor(current_state)
    new_energy = calc_energy(new_state, T, U)

    # メトロポリス基準
    dE = new_energy - current_energy
    if dE < 0 or random.random() < np.exp(-dE / T_current):
        current_state = new_state
        current_energy = new_energy
        if current_energy < best_energy:
            best_state = current_state.copy()
            best_energy = current_energy

# 結果をCSVに保存
df_output = pd.DataFrame([
    [family_types[i], ownership_types[j], building_types[k], floor_types[l], best_state[i, j, k, l]]
    for i in range(I) for j in range(J) for k in range(K) for l in range(L)
], columns=["家族類型", "所有", "建て方", "延床面積", "世帯数"])
df_output.to_csv("output_sa_3.csv", index=False)
