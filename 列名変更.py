import pandas as pd

# 合成人口データの読み込み
data1 = pd.read_csv("表19-2_博多区_編集.csv")
data2 = pd.read_csv("表25-2_博多区_編集.csv")

# 条件に基づき新しい列を変更
data1["所有"] = data1["所有"].apply(
    lambda x: 0 if x == "持ち家" 
              else 1 if x == "公営の借家" 
              else 2 if x == "民営の借家" 
              else 3 if x == "給与住宅"
              else 4
)

data1["建て方"] = data1["建て方"].apply(
    lambda x: 0 if x == "一戸建" 
              else 1 if x == "長屋建" 
              else 2 if x == "共同住宅" 
              else 3 
)

data1["面積"] = data1["面積"].apply(
    lambda x: 0 if x == "0～29m2" 
              else 1 if x == "30～49m2" 
              else 2 if x == "50～69m2" 
              else 3 if x == "70～99m2"
              else 4 if x == "100～149m2"
              else 5
)

data2["家族類型"] = data2["家族類型"].apply(
    lambda x: 0 if x == "夫婦のみの世帯" 
              else 1 if x == "夫婦と子供から成る世帯" 
              else 2 if x == "男親と子供から成る世帯" 
              else 3 if x == "女親と子供から成る世帯"
              else 4 if x == "夫婦と両親から成る世帯"
              else 5 if x == "夫婦とひとり親から成る世帯"
              else 6 if x == "夫婦，子供と両親から成る世帯"
              else 7 if x == "夫婦，子供とひとり親から成る世帯"
              else 8 if x == "単独世帯"
              else 9 if x == "夫婦と他の親族（親，子供を含まない）から成る世帯"
              else 10 if x == "夫婦，子供と他の親族（親を含まない）から成る世帯"
              else 11 if x == "夫婦，親と他の親族（子供を含まない）から成る世帯"
              else 12 if x == "夫婦，子供，親と他の親族から成る世帯"
              else 13 if x == "兄弟姉妹のみから成る世帯"
              else 14 if x == "他に分類されない世帯"
              else 15
)

data2["所有"] = data2["所有"].apply(
    lambda x: 0 if x == "持ち家" 
              else 1 if x == "公営の借家" 
              else 2 if x == "民営の借家" 
              else 3 if x == "給与住宅"
              else 4
)

# 結果を保存
data1.to_csv("19-2_博多区_SA用.csv", index=False, encoding="utf-8")
data2.to_csv("25-2_博多区_SA用.csv", index=False, encoding="utf-8")