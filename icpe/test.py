import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Dữ liệu
data = {
    'ok': [15, 10],
    'nok': [85, 90]
}

# Chuyển dữ liệu thành DataFrame
df = pd.DataFrame(data, index=['Row 1', 'Row 2'])

# Vẽ heatmap
plt.figure(figsize=(8, 4))  # Điều chỉnh kích thước bảng
ax = sns.heatmap(df, annot=True, cmap='YlOrRd', fmt='d', linewidths=0.5, linecolor='black')

# Tinh chỉnh tiêu đề và nhãn
plt.title('Heatmap of Values')
plt.xlabel('Columns')
plt.ylabel('Rows')

plt.show()
