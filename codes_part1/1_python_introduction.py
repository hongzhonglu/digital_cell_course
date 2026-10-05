# 基本数据类型 (Data Types)
# 整数 (int)
num_int = 100
print(f"整数: {num_int}, 类型: {type(num_int)}")

# 浮点数 (float)
num_float = 10.5
print(f"浮点数: {num_float}, 类型: {type(num_float)}")

# 布尔值 (bool)
is_active = True
print(f"布尔值: {is_active}, 类型: {type(is_active)}")

# 字符串 (str) - 可以使用单引号或双引号
greeting = 'Hello Python!'
print(f"字符串: {greeting}, 类型: {type(greeting)}")

# 字符串拼接 (连接)
full_message = greeting + " Welcome."
print(f"拼接后的字符串: {full_message}")





# 核心数据结构 (Data Structures)
# list
# 创建一个列表 (使用方括号 [])
my_list = [10, "apple", 3.14, True]
print(f"\n初始列表: {my_list}")

# 访问元素 (索引从 0 开始)
first_element = my_list[0]
print(f"第一个元素 (索引 0): {first_element}")

# 负数索引: -1 表示最后一个元素
last_element = my_list[-1]
print(f"最后一个元素 (索引 -1): {last_element}")

# 切片 (Slicing): [start:end] - 不包含 end 索引的元素
slice_data = my_list[1:3]  # 从索引 1 到索引 3 (不包括 3)
print(f"切片 [1:3]: {slice_data}")

# 修改元素 (列表是可变的)
my_list[1] = "banana"
print(f"修改后的列表: {my_list}")

# 添加元素到末尾
my_list.append(50)
print(f"添加元素后的列表: {my_list}")

# 删除元素 (根据值删除第一个匹配项)
my_list.remove(10)
print(f"删除 10 后的列表: {my_list}")

b = {'one': 1, 'two': 2, 'three': 3} # 字典

tup = (4,5,6,4) # 元组，固定长度 不可变

ss = {1,1,2,3,4,4} # 集合，无序且唯一



# 函数 (Functions)
# 使用 def 关键字定义一个函数，括号内是参数
def greet(name, time="morning"):
    """
    这是一个打招呼的函数。
    name 是必需参数，time 是可选参数 (有默认值)。
    """
    # 函数体需要缩进
    message = f"Good {time}, {name}! Welcome to Python world."
    return message # 使用 return 返回结果

# 调用函数 (传入必需参数)
greeting_message = greet("Professor")
print(f"\n函数调用结果 1: {greeting_message}")

# 调用函数 (传入可选参数)
greeting_message_evening = greet("Students", time="evening")
print(f"函数调用结果 2: {greeting_message_evening}")








# some popular package
import numpy as np
# 1. 创建一个 NumPy 数组 (通常是数据分析的起点)
a = np.array([10, 20, 30, 40, 50])
print(f"NumPy 数组 a: {a}")

# 2. 向量化操作 (快速且简洁)
# 数组中的每个元素都加上 5
b = a + 5
print(f"向量化加法 b: {b}")

# 3. 创建一个 2x3 的随机数矩阵
matrix = np.random.rand(2, 3)
print("\n随机 2x3 矩阵:")
print(matrix)

# 4. 矩阵的基本属性
print(f"矩阵的形状 (行, 列): {matrix.shape}")
print(f"矩阵的维度: {matrix.ndim}")



import pandas as pd
# 1. 使用字典创建 DataFrame (表格)
data = {
    'Name': ['Alice', 'Bob', 'Charlie'],
    'Age': [25, 30, 22],
    'City': ['NY', 'LA', 'SF']
}
df = pd.DataFrame(data)
print("初始 DataFrame (数据表):")
print(df)

# 2. 选择一列 (得到一个 Series)
ages = df['Age']
print(f"\n年龄列 (Series):\n{ages}")

# 3. 简单统计分析
average_age = df['Age'].mean()
print(f"\n平均年龄: {average_age:.2f}")

# 4. 数据筛选 (只看年龄大于 23 岁的行)
filtered_df = df[df['Age'] > 23]
print("\n筛选后的数据 (年龄 > 23):")
print(filtered_df)



### 绘图 ###
import matplotlib.pyplot as plt
import numpy as np # 经常与 NumPy 结合使用

plt.plot([1,2,3], [3, -4, 20])
plt.title('Awesome plot!')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()


# 1. 准备数据
# 生成 x 坐标 (0 到 10，共 50 个点)
x = np.linspace(0, 10, 50)
# 生成 y 坐标 (正弦函数)
y = np.sin(x)

# 2. 绘制折线图
plt.figure(figsize=(8, 4)) # 创建画布并设置大小
plt.plot(x, y, label='sin(x)', color='blue', linestyle='--')

# 3. 添加图表元素
plt.title('Simple Sine Wave Plot')
plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.legend() # 显示图例
plt.grid(True) # 显示网格线

# 4. 显示图表
# 注意：在 jupyter/colab 环境中，plt.show() 通常是可选的。
# 在脚本中运行，需要调用 plt.show()
# plt.show()
print("\n已生成 Matplotlib 图表，请在运行环境中查看或使用 plt.show() 显示。")