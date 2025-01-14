import pandas as pd
import numpy as np

# ------------------------------
# 1) データ読み込み
# ------------------------------
df_town = pd.read_csv('merged_town_population_data.csv')

df_city = pd.read_csv('formatted_city_population_data.csv')

# ------------------------------
# 2) 差を計算
# ------------------------------
df_town['Diff'] = df_town['NonResSumPopulation'] - df_town['NonResSumHouseholds']

# ------------------------------
# 3) 差=0 / 差=1 / 差>=2 に分割
# ------------------------------
df_zero = df_town.query('Diff == 0').copy()
df_one  = df_town.query('Diff == 1').copy()
df_rest = df_town.query('Diff >= 2').copy()

# ------------------------------
# 4) 差=0: 全部単独世帯（他のFamilyTypeは0に設定）
# ------------------------------
def handle_diff_zero(group):
    result = []
    for _, row in group.iterrows():
        if row['FamilyType'] == '単独世帯':
            result.append({
                'Town': row['Town'],
                'FamilyType': row['FamilyType'],
                'NonResHouseholds': row['NonResSumHouseholds'],
                'NonResPopulation': row['NonResSumPopulation']
            })
        else:
            result.append({
                'Town': row['Town'],
                'FamilyType': row['FamilyType'],
                'NonResHouseholds': 0,
                'NonResPopulation': 0
            })
    return pd.DataFrame(result)

df_zero_adjusted = df_zero.groupby('Town', group_keys=False).apply(handle_diff_zero).reset_index(drop=True)

# ------------------------------
# 5) 差=1: 「2人世帯=1 + 残り単独」
# ------------------------------
def handle_diff_one(group):
    result = []
    remaining_households = group['NonResSumHouseholds'].iloc[0]
    remaining_population = group['NonResSumPopulation'].iloc[0]

    for _, row in group.iterrows():
        if row['FamilyType'] == '夫婦のみの世帯':
            two_person_households = 1
            two_person_population = 2
            result.append({
                'Town': row['Town'],
                'FamilyType': row['FamilyType'],
                'NonResHouseholds': two_person_households,
                'NonResPopulation': two_person_population
            })
            remaining_households -= two_person_households
            remaining_population -= two_person_population
        elif row['FamilyType'] == '単独世帯':
            result.append({
                'Town': row['Town'],
                'FamilyType': row['FamilyType'],
                'NonResHouseholds': max(remaining_households, 0),
                'NonResPopulation': max(remaining_population, 0)
            })
            remaining_households = 0
            remaining_population = 0
        else:
            result.append({
                'Town': row['Town'],
                'FamilyType': row['FamilyType'],
                'NonResHouseholds': 0,
                'NonResPopulation': 0
            })
    return pd.DataFrame(result)

df_one_adjusted = df_one.groupby('Town', group_keys=False).apply(handle_diff_one).reset_index(drop=True)

# ------------------------------
# 6) 消費合計を集計し、市区町村ファイルから差し引き
# ------------------------------
# (a) 単独世帯に関して
used_single = df_zero_adjusted.query("FamilyType == '単独世帯'")['NonResHouseholds'].sum() + \
              df_one_adjusted.query("FamilyType == '単独世帯'")['NonResHouseholds'].sum()
used_single_pop = used_single

# (b) 2人世帯に関して (「夫婦のみ世帯」として消費すると仮定)
used_twop = df_one_adjusted.query("FamilyType == '夫婦のみの世帯'")['NonResHouseholds'].sum()
used_twop_pop = 2 * used_twop

print("==== (差=0,1) の町丁目 確定割当 ====")
print("単独世帯  消費 =", used_single, ", 人口 =", used_single_pop)
print("2人世帯  消費 =", used_twop, ", 人口 =", used_twop_pop)

# ------------------------------
# 7) 市区町村ファイルの合計を更新(残り枠計算)
# ------------------------------
df_city_remain = df_city.copy()

def subtract_usage(row):
    """
    row: city_familytype.csv の1行 (FamilyType, NonResHouseholds, NonResPopulation)
    単独世帯(=1人) → 'FamilyType'=='単独世帯'
    2人世帯(=夫婦のみ？) → 'FamilyType'=='夫婦のみの世帯'
    """
    ft = row['FamilyType']
    if ft == '単独世帯':
        new_hh = max(row['NonResHouseholds'] - used_single, 0)
        new_pop = max(row['NonResPopulation'] - used_single_pop, 0)
    elif ft == '夫婦のみの世帯':
        new_hh = max(row['NonResHouseholds'] - used_twop, 0)
        new_pop = max(row['NonResPopulation'] - used_twop_pop, 0)
    else:
        new_hh = row['NonResHouseholds']  # 変化なし
        new_pop = row['NonResPopulation']
    return pd.Series([new_hh, new_pop], index=['NonResHouseholds','NonResPopulation'])

df_city_remain[['NonResHouseholds','NonResPopulation']] = df_city_remain.apply(subtract_usage, axis=1)

print("\n==== 市区町村レベル 残り枠 (df_city_remain) ====")
print(df_city_remain)

# ------------------------------
# 8) 差 >= 2 の町丁目をIPFで割り当て
# ------------------------------
# 初期値を設定し IPF の準備
df_rest['TotalHouseholds_sum'] = df_rest.groupby('Town')['TotalHouseholds'].transform('sum')
df_rest['Ratio_HH'] = df_rest['TotalHouseholds'] / df_rest['TotalHouseholds_sum']

if 'NonResSumHouseholds' in df_rest.columns:
    df_rest['init_nonres_households'] = df_rest['NonResSumHouseholds'] * df_rest['Ratio_HH']
else:
    raise KeyError("Column 'NonResSumHouseholds' not found in df_rest")

df_rest['TotalPopulation_sum'] = df_rest.groupby('Town')['TotalPopulation'].transform('sum')
df_rest['Ratio_Pop'] = df_rest['TotalPopulation'] / df_rest['TotalPopulation_sum']

if 'NonResSumPopulation' in df_rest.columns:
    df_rest['init_nonres_population'] = df_rest['NonResSumPopulation'] * df_rest['Ratio_Pop']
else:
    raise KeyError("Column 'NonResSumPopulation' not found in df_rest")

# ピボットデータ生成
df_pivot_hh = df_rest.pivot_table(index='Town', columns='FamilyType', values='init_nonres_households', aggfunc='sum', fill_value=0)
df_pivot_pop = df_rest.pivot_table(index='Town', columns='FamilyType', values='init_nonres_population', aggfunc='sum', fill_value=0)

row_targets_hh = df_rest.groupby('Town')['NonResSumHouseholds'].first()
row_targets_pop = df_rest.groupby('Town')['NonResSumPopulation'].first()

col_targets_hh = df_city_remain.set_index('FamilyType')['NonResHouseholds']
col_targets_pop = df_city_remain.set_index('FamilyType')['NonResPopulation']

# IPFの実装
def ipf_2d(init_mat, row_targets, col_targets, max_iter=500, tol=1e-12):
    mat = init_mat.copy()
    for _ in range(max_iter):
        # 行方向調整
        row_sums = mat.sum(axis=1)
        mat = (mat.T * (row_targets / np.clip(row_sums, tol, None))).T

        # 列方向調整
        col_sums = mat.sum(axis=0)
        mat = mat * (col_targets / np.clip(col_sums, tol, None))

        # 負の値をクリップ
        mat = np.clip(mat, 0, None)

        # 収束条件
        if np.allclose(mat.sum(axis=1), row_targets, atol=tol) and \
           np.allclose(mat.sum(axis=0), col_targets, atol=tol):
            break

    return mat

# 丸め後の調整
def adjust_after_rounding(mat, row_targets, col_targets):
    mat_rounded = np.rint(mat).astype(int)

    # 行方向誤差
    row_diff = row_targets - mat_rounded.sum(axis=1)
    for i, diff in enumerate(row_diff):
        if diff > 0:
            idx = np.argmax(mat[i, :] - mat_rounded[i, :])
            mat_rounded[i, idx] += diff
        elif diff < 0:
            idx = np.argmin(mat[i, :] - mat_rounded[i, :])
            mat_rounded[i, idx] += diff

    # 列方向誤差
    col_diff = col_targets - mat_rounded.sum(axis=0)
    for j, diff in enumerate(col_diff):
        if diff > 0:
            idx = np.argmax(mat[:, j] - mat_rounded[:, j])
            mat_rounded[idx, j] += diff
        elif diff < 0:
            idx = np.argmin(mat[:, j] - mat_rounded[:, j])
            mat_rounded[idx, j] += diff

    return np.clip(mat_rounded, 0, None)

# IPF実行
adjusted_hh = ipf_2d(df_pivot_hh.values, row_targets_hh.values, col_targets_hh.values)
adjusted_pop = ipf_2d(df_pivot_pop.values, row_targets_pop.values, col_targets_pop.values)

# 丸め後調整
adjusted_hh = adjust_after_rounding(adjusted_hh, row_targets_hh.values, col_targets_hh.values)
adjusted_pop = adjust_after_rounding(adjusted_pop, row_targets_pop.values, col_targets_pop.values)

# データフレーム化
df_adjusted_hh = pd.DataFrame(adjusted_hh, index=df_pivot_hh.index, columns=df_pivot_hh.columns).reset_index().melt(
    id_vars='Town',
    var_name='FamilyType',
    value_name='NonResHouseholds'
)
df_adjusted_pop = pd.DataFrame(adjusted_pop, index=df_pivot_pop.index, columns=df_pivot_pop.columns).reset_index().melt(
    id_vars='Town',
    var_name='FamilyType',
    value_name='NonResPopulation'
)

# 合体
df_adjusted = pd.merge(df_adjusted_hh, df_adjusted_pop, on=['Town', 'FamilyType'], how='left')

# ------------------------------
# 9) 差=0,1 の町丁目と差>=2 の結果を結合
# ------------------------------
# 結合
df_final = pd.concat([df_zero_adjusted, df_one_adjusted, df_adjusted], ignore_index=True)

# 元の順序を維持するためにソート
ordered_towns = df_town['Town'].drop_duplicates()
df_final['Town'] = pd.Categorical(df_final['Town'], categories=ordered_towns, ordered=True)
df_final = df_final.sort_values(['Town', 'FamilyType']).reset_index(drop=True)

# 出力形式を調整
df_final = pd.merge(
    df_town[['Town', 'NonResSumHouseholds', 'NonResSumPopulation']].drop_duplicates(),
    df_final.groupby(['Town', 'FamilyType'], observed=True)[['NonResHouseholds', 'NonResPopulation']].sum().reset_index(),
    on='Town',
    how='left'
)

# 結果をCSV出力
df_final.to_csv('非住宅世帯数_人口.csv', index=False, encoding='utf-8-sig')