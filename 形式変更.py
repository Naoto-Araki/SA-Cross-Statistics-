import pandas as pd

# CSVファイルを読み込む
# file_path = 'updated_town_population_data.csv'  # CSVファイルのパスを指定
# file_path = '町丁目所有_博多区.csv'  # CSVファイルのパスを指定
# file_path = '町丁目面積_博多区.csv'  # CSVファイルのパスを指定
file_path = '町丁目建て方_博多区.csv'  # CSVファイルのパスを指定
df = pd.read_csv(file_path)

# .1 が付いている列と付いていない列を分ける
population_columns = [col for col in df.columns if '.1' in col]
household_columns = [col for col in df.columns if col not in population_columns and col != '町丁目']

# .1 を削除して列名を統一
population_columns_cleaned = [col.replace('.1', '').strip() for col in population_columns]

# TotalHouseholdsのデータフレーム作成
df_households = df[['町丁目'] + household_columns].melt(
    id_vars=['町丁目'], var_name='FamilyType', value_name='TotalHouseholds'
)

# TotalPopulationのデータフレーム作成
df_population = df[['町丁目'] + population_columns]
df_population.columns = ['町丁目'] + population_columns_cleaned
df_population = df_population.melt(
    id_vars=['町丁目'], var_name='FamilyType', value_name='TotalPopulation'
)

# TotalHouseholdsとTotalPopulationを結合
result = pd.merge(df_households, df_population, on=['町丁目', 'FamilyType'])

# 町丁目の順序を元データに基づいて保持
result['町丁目'] = pd.Categorical(result['町丁目'], categories=df['町丁目'].unique(), ordered=True)

# FamilyTypeの順序を入力データの列順に基づいて保持
family_type_order = household_columns  # household_columns の順序を保持
result['FamilyType'] = pd.Categorical(result['FamilyType'], categories=family_type_order, ordered=True)

# 並べ替え
result = result.sort_values(['町丁目', 'FamilyType']).reset_index(drop=True)

# 結果を表示
print(result)

# 必要に応じてCSVファイルとして保存
# result.to_csv('formatted_town_population_data.csv', index=False)
# result.to_csv('町丁目所有_博多区.csv', index=False)
# result.to_csv('町丁目面積_博多区.csv', index=False)
result.to_csv('町丁目建て方_博多区.csv', index=False)
