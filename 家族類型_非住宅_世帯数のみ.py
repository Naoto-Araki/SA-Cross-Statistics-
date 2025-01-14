import pandas as pd
import numpy as np

# === CSVファイルの読み込み ===
city_file = "formatted_city_population_data_非住宅.csv"  # 市区町村データ
town_file = "merged_town_population_data_SA.csv"  # 町丁目データ

df_city = pd.read_csv(city_file)
df_town = pd.read_csv(town_file)

# === 元の順序を保持する列を追加 ===
df_town['OriginalOrder'] = df_town.index  # 入力時の行インデックスを記録

# === STEP 1: 非住宅世帯数と非住宅人口が一致している地域を特定 ===
df_town['IsSingleOnly'] = df_town['NonResSumHouseholds'] == df_town['NonResSumPopulation']

def assign_single_only(row):
    if row['IsSingleOnly'] and row['FamilyType'] == '単独世帯':
        return round(row['NonResSumHouseholds'])
    elif row['IsSingleOnly']:
        return 0
    else:
        return np.nan

df_town['NonResHouseholds_Init'] = df_town.apply(assign_single_only, axis=1)

# === STEP 2: 初期割り当てを実施（非一致地域） ===
df_town_non_matching = df_town[df_town['NonResHouseholds_Init'].isna()].copy()

def calculate_family_ratio(row):
    total_households = df_town_non_matching[df_town_non_matching['Town'] == row['Town']]['TotalHouseholds'].sum()
    return row['TotalHouseholds'] / total_households

df_town_non_matching['FamilyRatio'] = df_town_non_matching.apply(calculate_family_ratio, axis=1)

df_town_non_matching['NonResHouseholds_Init'] = (
    df_town_non_matching['NonResSumHouseholds'] * df_town_non_matching['FamilyRatio']
).round()

df_town.update(df_town_non_matching[['Town', 'FamilyType', 'NonResHouseholds_Init']])

# === STEP 3: 町丁目単位の整合性確認と再調整 ===
def adjust_town_discrepancy(df):
    """
    町丁目単位で「非住宅世帯数の合計 = 町丁目全体の非住宅世帯数」を調整。
    """
    for town, group in df.groupby('Town'):
        total_nonres = group['NonResHouseholds_Init'].sum()
        target_nonres = group['NonResSumHouseholds'].iloc[0]

        # 誤差計算
        discrepancy = round(target_nonres - total_nonres)

        if discrepancy != 0 and total_nonres > 0:
            for idx in group.index:
                # 固定地域をスキップ
                if df.loc[idx, 'IsSingleOnly']:
                    continue

                init_value = df.loc[idx, 'NonResHouseholds_Init']
                adjust_value = discrepancy * (init_value / total_nonres)
                df.loc[idx, 'NonResHouseholds_Init'] += round(adjust_value)

        # 微調整で最終的なズレを修正
        total_nonres_corrected = df.loc[group.index, 'NonResHouseholds_Init'].sum()
        final_discrepancy = target_nonres - total_nonres_corrected
        if abs(final_discrepancy) > 0:
            max_idx = group[~group['IsSingleOnly']]['NonResHouseholds_Init'].idxmax()
            df.loc[max_idx, 'NonResHouseholds_Init'] += round(final_discrepancy)

    # 負の値をクリップ
    df['NonResHouseholds_Init'] = df['NonResHouseholds_Init'].clip(lower=0)
    return df

def adjust_city_discrepancy(df, df_city):
    """
    市区町村単位で整合性を調整する。
    """
    df_city_grouped = df.groupby('FamilyType')['NonResHouseholds_Init'].sum().reset_index()
    df_city_grouped = df_city_grouped.merge(df_city, on='FamilyType')
    df_city_grouped['Difference'] = (
        df_city_grouped['NonResHouseholds_Target'] - df_city_grouped['NonResHouseholds_Init']
    )

    for _, row in df_city_grouped.iterrows():
        family_type = row['FamilyType']
        difference = row['Difference']

        if difference != 0:
            family_df = df[(df['FamilyType'] == family_type) & (~df['IsSingleOnly'])]
            total_init = family_df['NonResHouseholds_Init'].sum()

            for idx in family_df.index:
                init_value = df.loc[idx, 'NonResHouseholds_Init']
                adjust_value = difference * (init_value / total_init) if total_init > 0 else 0
                df.loc[idx, 'NonResHouseholds_Init'] += round(adjust_value)

    return df

def finalize_discrepancy(df, df_city):
    """
    最終的に町丁目と市区町村の整合性を完全一致させる。
    """
    # 市区町村単位のズレを確認
    city_discrepancy = df.groupby('FamilyType')['NonResHouseholds_Init'].sum().reset_index()
    city_discrepancy = city_discrepancy.merge(df_city, on='FamilyType')
    city_discrepancy['Error'] = city_discrepancy['NonResHouseholds_Target'] - city_discrepancy['NonResHouseholds_Init']

    # 差分を調整
    for _, row in city_discrepancy.iterrows():
        family_type = row['FamilyType']
        error = int(row['Error'])  # 整数化

        while error != 0:
            # 差分を調整する対象の行を選択
            if error > 0:
                # 誤差がプラスの場合: 最大値を持つ町丁目に1を追加
                target_df = df[df['FamilyType'] == family_type]
                if target_df.empty:
                    print(f"Error: No data found for family_type={family_type}")
                    break
                target_idx = target_df['NonResHouseholds_Init'].idxmax()
                df.loc[target_idx, 'NonResHouseholds_Init'] += 1
                error -= 1
            elif error < 0:
                # 誤差がマイナスの場合: 最小値を持つ町丁目から1を減少
                target_df = df[(df['FamilyType'] == family_type) & (df['NonResHouseholds_Init'] > 0)]
                if target_df.empty:
                    print(f"Error: No data found for family_type={family_type} with NonResHouseholds_Init > 0")
                    break
                target_idx = target_df['NonResHouseholds_Init'].idxmin()
                df.loc[target_idx, 'NonResHouseholds_Init'] -= 1
                error += 1

    return df


def adjust_town_and_city_discrepancy_loop(df, df_city, max_iterations=100):
    """
    町丁目単位と市区町村単位の整合性を交互に調整。
    """
    for iteration in range(max_iterations):
        # 町丁目単位の調整
        df = adjust_town_discrepancy(df)

        # 市区町村単位の調整
        df = adjust_city_discrepancy(df, df_city)

        # 整合性確認
        town_discrepancy = df.groupby('Town', group_keys=False).apply(
            lambda x: round(x['NonResHouseholds_Init'].sum() - x['NonResSumHouseholds'].iloc[0])
        ).sum()
        city_discrepancy = df.groupby('FamilyType')['NonResHouseholds_Init'].sum().reset_index()
        city_discrepancy = city_discrepancy.merge(df_city, on='FamilyType')
        total_city_discrepancy = (city_discrepancy['NonResHouseholds_Target'] - 
                                  city_discrepancy['NonResHouseholds_Init']).sum()

        if town_discrepancy == 0 and total_city_discrepancy == 0:
            print(f"整合性が取れました: {iteration+1} 回のループで達成")
            break
    else:
        print("最大ループ回数に達しましたが、完全な整合性は取れていない可能性があります。")

    return df

def finalize_town_discrepancy(df):
    """
    町丁目単位の最終整合性を保証する。
    各町丁目で「非住宅世帯数の合計 = 町丁目全体の非住宅世帯数」を確認し、ズレがあれば修正。
    """
    for town, group in df.groupby('Town'):
        total_nonres = group['NonResHouseholds_Init'].sum()
        target_nonres = group['NonResSumHouseholds'].iloc[0]

        # ズレを計算
        discrepancy = round(target_nonres - total_nonres)

        if discrepancy != 0:
            # 誤差がプラスの場合（不足を補う）
            if discrepancy > 0:
                max_idx = group['NonResHouseholds_Init'].idxmax()
                df.loc[max_idx, 'NonResHouseholds_Init'] += discrepancy
            # 誤差がマイナスの場合（超過を減らす）
            elif discrepancy < 0:
                max_idx = group[group['NonResHouseholds_Init'] > 0]['NonResHouseholds_Init'].idxmax()
                df.loc[max_idx, 'NonResHouseholds_Init'] += discrepancy  # 減算なので加算符号そのまま

    # 負の値をクリップしておく
    df['NonResHouseholds_Init'] = df['NonResHouseholds_Init'].clip(lower=0)
    return df

df_town = adjust_town_and_city_discrepancy_loop(df_town, df_city)
df_town = finalize_discrepancy(df_town, df_city)  # 市区町村単位の微調整
df_town = finalize_town_discrepancy(df_town)  # 町丁目単位の最終調整

# === 元の順序で並べ替え ===
df_town = df_town.sort_values('OriginalOrder').drop(columns=['OriginalOrder'])

# === 結果を保存 ===
df_town['NonResHouseholds_Init'] = df_town['NonResHouseholds_Init'].astype(int)
df_town.to_csv("town_adjusted_final.csv", index=False)
print("処理完了: 結果は 'town_adjusted_final.csv' に保存されました。")
