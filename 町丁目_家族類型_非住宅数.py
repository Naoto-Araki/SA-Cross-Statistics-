import pandas as pd

# CSVファイルの読み込み（適宜ファイル名を変更してください）
file1_path = 'formatted_town_population_data.csv'  # 最初の表のファイル
file2_path = '非住宅数_町丁目_修正.csv'  # 2つ目の表のファイル

df1 = pd.read_csv(file1_path)  # 最初の表
df2 = pd.read_csv(file2_path)  # 2つ目の表

# 町丁目を基準に結合
result = pd.merge(df1, df2, on='町丁目', how='left')

# 結果を表示
print(result)

# 必要に応じてCSVファイルに保存
result.to_csv('merged_town_population_data.csv', index=False)
