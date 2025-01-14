import pandas as pd

# CSVファイルの読み込み
# file_path = 'town_population_data.csv'  # CSVファイルのパス
# file_path = '非住宅数_町丁目_修正.csv'  # CSVファイルのパス
# file_path = '町丁目所有_博多区.csv'  # CSVファイルのパス
# file_path = '町丁目面積_博多区.csv'  # CSVファイルのパス
file_path = '町丁目建て方_博多区.csv'  # CSVファイルのパス
df = pd.read_csv(file_path)

# 2列を結合して '町丁目' 列を作成
df['町丁目'] = df.apply(
    lambda row: row['大字・町名'] + row['字・丁目名'] if pd.notna(row['字・丁目名']) else row['大字・町名'],
    axis=1
)

# 元の2列を削除
df.drop(['大字・町名', '字・丁目名'], axis=1, inplace=True)

# '町丁目' を一番左に移動
columns = ['町丁目'] + [col for col in df.columns if col != '町丁目']
df = df[columns]

# 結果を表示
print(df)

# 必要に応じてCSVファイルに保存
# df.to_csv('updated_town_population_data.csv', index=False)
# df.to_csv('非住宅数_町丁目_修正.csv', index=False)
# df.to_csv('町丁目所有_博多区.csv', index=False)
# df.to_csv('町丁目面積_博多区.csv', index=False)
df.to_csv('町丁目建て方_博多区.csv', index=False)
