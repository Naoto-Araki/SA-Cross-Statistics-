import pandas as pd

# === CSVファイルを読み込む ===
# input_file = "町丁目所有_博多区.csv"  # 入力ファイル名
# output_file = "町丁目所有_博多区.csv"  # 出力ファイル名

# input_file = "町丁目面積_博多区.csv"  # 入力ファイル名
# output_file = "町丁目面積_博多区.csv"  # 出力ファイル名

input_file = "町丁目建て方_博多区.csv"  # 入力ファイル名
output_file = "町丁目建て方_博多区.csv"  # 出力ファイル名

# DataFrameを読み込み
df = pd.read_csv(input_file)

# === "-" を "0" に置き換える ===
df.replace("-", 0, inplace=True)

# === 必要に応じて数値型に変換 ===
# 文字列だった列を数値型に変換
df = df.apply(pd.to_numeric, errors='ignore')  # 変換可能な列だけ数値型に

# === 結果を保存 ===
df.to_csv(output_file, index=False)
print(f"変換が完了しました。結果を '{output_file}' に保存しました。")
