import pandas as pd
import numpy as np
import random

# データの読み込み
# 町丁目単位データ
town_family = pd.read_csv("town_adjusted_final.csv")  # 家族類型 × 世帯数
town_ownership = pd.read_csv("町丁目所有_博多区.csv")  # 所有 × 世帯数
town_building = pd.read_csv("町丁目建て方_博多区.csv")  # 建て方 × 世帯数
town_area = pd.read_csv("町丁目面積_博多区.csv")  # 面積 × 世帯数

# 市区町村単位データ
city_family_ownership = pd.read_csv("表25-2_博多区_編集_2.csv")  # 家族類型 × 所有
city_building_area = pd.read_csv("表19-2_博多区_編集.csv")  # 所有 × 建て方 × 面積

# 分配型丸め処理
def round_with_total_constraint(values, target_sum):
    """
    values: ndarray or pandas.Series
        丸め処理対象の値のリスト（小数点を含む）。
    target_sum: int
        丸め後に一致させるべき合計値。
    """
    # numpy 配列に変換して処理を統一
    values_array = np.asarray(values)
    rounded_values = np.floor(values_array).astype(int)

    # 差分計算
    difference = target_sum - rounded_values.sum()
    fractional_parts = values_array - np.floor(values_array)

    # 調整するインデックスを取得
    indices_to_adjust = np.argsort(-fractional_parts)[:difference]

    # インデックスを適用して調整
    rounded_values[indices_to_adjust] += 1

    # 元の形式が pandas.Series の場合は再変換
    if isinstance(values, pd.Series):
        rounded_values = pd.Series(rounded_values, index=values.index)

    return rounded_values

# 間借り世帯の建て方推計
def estimate_rent_building(town, town_ownership, city_building_area):
    # 町丁目の間借り世帯数を取得
    rent_households = town_ownership.loc[
        (town_ownership['Town'] == town) & (town_ownership['OwnershipType'] == '間借り'),
        'TotalHouseholds'
    ].values[0]

    # 市区町村レベルの「間借り × 建て方」の分布を計算
    rent_building_distribution = city_building_area[city_building_area['OwnershipType'] == '間借り']
    rent_building_ratios = rent_building_distribution['TotalHouseholds'] / rent_building_distribution['TotalHouseholds'].sum()

    # 町丁目の間借り世帯数を建て方別に割り振り
    estimated_rent_building = rent_building_ratios * rent_households
    estimated_rent_building = round_with_total_constraint(estimated_rent_building, rent_households)

    # 推計結果をDataFrameに格納
    rent_building_df = rent_building_distribution.copy()
    rent_building_df['TotalHouseholds'] = estimated_rent_building
    rent_building_df['Town'] = town

    return rent_building_df[['Town', 'BuildingType', 'TotalHouseholds']]

# 初期割り当て
def initialize_solution(town, town_family, town_ownership, town_building, town_area, city_family_ownership, city_building_area, estimated_rent_building):
    solution = []

    # 主世帯の建て方別周辺表を取得
    building_totals = town_building.loc[town_building['Town'] == town].set_index('BuildingType')['TotalHouseholds']

    # 間借り世帯の建て方別推計結果を取得
    estimated_rent_building_totals = estimated_rent_building.groupby('BuildingType')['TotalHouseholds'].sum()

    # 主世帯 + 間借り世帯の建て方別合算
    combined_building_totals = building_totals.add(estimated_rent_building_totals, fill_value=0)

    for family_type in town_family['HouseholdType'].unique():
        town_family_total = town_family.loc[(town_family['Town'] == town) & (town_family['HouseholdType'] == family_type), 'TotalHouseholds'].values[0]

        family_ownership_ratios = city_family_ownership[city_family_ownership['HouseholdType'] == family_type]
        ownership_counts = []

        for _, row in family_ownership_ratios.iterrows():
            ratio = row['TotalHouseholds'] / family_ownership_ratios['TotalHouseholds'].sum()
            ownership_counts.append(town_family_total * ratio)

        ownership_counts = round_with_total_constraint(np.array(ownership_counts), town_family_total)

        for i, (_, row) in enumerate(family_ownership_ratios.iterrows()):
            ownership = row['OwnershipType']
            building_area_ratios = city_building_area[city_building_area['OwnershipType'] == ownership]
            building_counts = []

            for _, ba_row in building_area_ratios.iterrows():
                ba_ratio = ba_row['TotalHouseholds'] / building_area_ratios['TotalHouseholds'].sum()
                building_counts.append(ownership_counts[i] * ba_ratio)

            building_counts = round_with_total_constraint(np.array(building_counts), ownership_counts[i])

            for j, ba_row in enumerate(building_area_ratios.itertuples()):
                solution.append({
                    'Town': town,
                    'FamilyType': family_type,
                    'OwnershipType': ownership,
                    'BuildingType': ba_row.BuildingType,
                    'AreaRange': ba_row.AreaRange,
                    'Households': building_counts[j]
                })

    # # 合算した建て方別周辺表を solution に反映
    # for building_type, total in combined_building_totals.items():
    #     solution.append({
    #         'Town': town,
    #         'FamilyType': None,  # 必要に応じて家族類型を設定
    #         'OwnershipType': None,  # 必要に応じて所有を設定
    #         'BuildingType': building_type,
    #         'AreaRange': None,  # 面積は後で割り当てる
    #         'Households': total
    #     })

    return pd.DataFrame(solution)

# シミュレーテッドアニーリング
def simulated_annealing(town, solution, town_family, town_ownership, town_building, town_area):
    def calculate_error(solution, town_family, town_ownership, town_building, town_area):
        family_totals = town_family.loc[town_family['Town'] == town].groupby('HouseholdType')['TotalHouseholds'].sum()
        family_sums = solution.groupby('FamilyType')['Households'].sum()
        family_error = np.abs(family_totals - family_sums).sum()

        ownership_totals = town_ownership.loc[town_ownership['Town'] == town].groupby('OwnershipType')['TotalHouseholds'].sum()
        ownership_sums = solution.groupby('OwnershipType')['Households'].sum()
        ownership_error = np.abs(ownership_totals - ownership_sums).sum()

        area_totals = town_area.loc[town_area['Town'] == town].groupby('AreaRange')['TotalHouseholds'].sum()
        area_sums = solution.groupby('AreaRange')['Households'].sum()
        area_error = np.abs(area_totals - area_sums).sum()

        building_totals = town_building.loc[town_building['Town'] == town].groupby('BuildingType')['TotalHouseholds'].sum()
        building_sums = solution.groupby('BuildingType')['Households'].sum()
        building_error = np.abs(building_totals - building_sums).sum()

        total_error = ownership_error + area_error + building_error  # family_errorは常にゼロにする
        return total_error, {
            'family_error': family_error,
            'ownership_error': ownership_error,
            'area_error': area_error,
            'building_error': building_error
        }

    # 初期設定
    initial_temp = 0.5
    final_temp = 0.1
    max_iterations = int(town_family.loc[town_family['Town'] == town, 'TotalHouseholds'].sum() * 1000)
    current_error, error_details = calculate_error(solution, town_family, town_ownership, town_building, town_area)

    for iteration in range(max_iterations):
        T_current = initial_temp * (final_temp / initial_temp) ** (iteration / max_iterations)

        # ランダムに2つのセルを選択して入れ替え操作を試行
        idx1, idx2 = random.sample(range(len(solution)), 2)

        # 入れ替え操作の制限条件
        if solution.loc[idx1, 'FamilyType'] != solution.loc[idx2, 'FamilyType']:
            continue  # 家族類型が異なる場合はスキップ

        if solution.loc[idx1, 'Households'] > 0:
            solution.loc[idx1, 'Households'] -= 1
            solution.loc[idx2, 'Households'] += 1

            new_error, _ = calculate_error(solution, town_family, town_ownership, town_building, town_area)

            # メトロポリス基準による受け入れ判定
            if new_error < current_error or random.random() < np.exp((current_error - new_error) / T_current):
                current_error = new_error
            else:
                # 変更を元に戻す
                solution.loc[idx1, 'Households'] += 1
                solution.loc[idx2, 'Households'] -= 1

    # 最終的な誤差を再計算
    final_error, final_error_details = calculate_error(solution, town_family, town_ownership, town_building, town_area)

    return solution, final_error, final_error_details



# 最適化プロセス
towns = town_family['Town'].unique()
optimized_solutions = []
error_summary = []

for town in towns:
    print(f"最適化中: {town}")
    try:
        estimated_rent_building = estimate_rent_building(town, town_ownership, city_building_area)
        initial_solution = initialize_solution(
            town, town_family, town_ownership, town_building, town_area,
            city_family_ownership, city_building_area, estimated_rent_building
        )

        # 最適化プロセスを実行
        optimized_solution, final_error, error_details = simulated_annealing(
            town, initial_solution, town_family, town_ownership, town_building, town_area
        )

        optimized_solution['Town'] = town
        optimized_solutions.append(optimized_solution)

        # 誤差を記録
        error_summary.append({
            'Town': town,
            'TotalError': final_error,
            **error_details
        })

    except Exception as e:
        print(f"{town} の処理中にエラーが発生しました: {e}")

# 結果を保存
final_solution = pd.concat(optimized_solutions)
final_solution.to_csv("4次元表_町丁目_3.csv", index=False)

# 誤差情報を保存
error_df = pd.DataFrame(error_summary)
error_df.to_csv("4次元表_町丁目_誤差_3.csv", index=False)
