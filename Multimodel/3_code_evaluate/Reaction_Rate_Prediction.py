import os
import json
import statistics
import re, pandas as pd

# 定义一个解析范围的函数
def parse_range(range_str):
    # 将字符串中的'-'替换为空格
    range_str = range_str.replace("-", " ")
    # 使用正则表达式找到字符串中的数字
    numbers = re.findall(r'[+-]?[0-9]*\.?[0-9]+', range_str)
    # 如果找到的数字小于2，返回None
    if len(numbers) < 2:
        return None
    # 将找到的数字转换为浮点数
    low, high = map(float, numbers[:2])
    # 返回最低值和最高值
    return low, high


def calculate_overlap(range_dict):
    '''
    计算range_dict中'answer'和'gold'之间的重合值

    参数：
    range_dict (dict): 一个包含'answer'和'gold'的range的字典

    返回：
    float: 重合值
    '''

    answer_range = parse_range(range_dict['answer'])

    gold_range = parse_range(range_dict['gold'])
    # answer_range = parse_range(range_dict['answer'])
    # gold_range = parse_range(range_dict['target'])
    if answer_range == None:
        return 0

    # 计算交集
    intersection_low = max(answer_range[0], gold_range[0])
    intersection_high = min(answer_range[1], gold_range[1])
    intersection = max(0, intersection_high - intersection_low)

    # 计算并集
    union_low = min(answer_range[0], gold_range[0])
    union_high = max(answer_range[1], gold_range[1])
    union = union_high - union_low

    # 计算重合值
    if union == 0:
        return 0
    overlap_value = intersection / union

    return overlap_value

# 定义一个函数，将jsonl文件转换为json对象列表
def convert_jsonl_to_json(jsonl_file_path):
    # 创建一个空列表，用于存储json对象
    data_list = []
    # 以读取模式打开jsonl文件
    with open(jsonl_file_path, 'r', encoding='utf-8') as file:
        # 遍历文件的每一行
        for line in file:
            # 将每一行转换为json对象
            json_object = json.loads(line.strip())
            # 将json对象添加到列表中
            data_list.append(json_object)
    # 返回列表
    return data_list

# 定义一个函数，列出指定文件夹中的json文件
def list_json_files(folder_path):
    # 创建一个空列表，用于存储json文件名
    file_list = []
    # 遍历文件夹中的所有文件
    for filename in os.listdir(folder_path):
        # 如果文件名以.json结尾，将其添加到列表中
        if filename.endswith(".json"):
            file_list.append(str(filename))
    # 返回列表
    return file_list


def reaction_rate(root_path, excel_path):
    # 定义文件夹路径
    # foler_path = r"C:\Users\xiangli67\Downloads\3880数据集\o3-mini\范围"
    # 获取文件夹中的json文件列表
    foler_path = os.path.join(root_path, "Reaction_Rate_Prediction")
    file_list = list_json_files(foler_path)
    result_data = []

    # 遍历文件列表
    for name in file_list:
        # 定义文件路径
        path = foler_path + "\\" + name
        # 将json文件转换为字典
        data = convert_jsonl_to_json(path)
        # 定义空列表，用于存储每个样本的overlap值
        results = []
        # 遍历字典中的每个样本
        for i in range(len(data)):
            # 计算每个样本的overlap值
            overlap = calculate_overlap(data[i])
            # 将overlap值添加到results列表中
            results.append(overlap)

        # 计算所有样本的overlap平均值
        result = statistics.mean(results)
        result_data.append([name, result])
        # 打印文件名和overlap平均值
        print(name, f"Overlap Score: {result}")

    # 创建DataFrame
    df = pd.DataFrame(result_data, columns=["File Name", "Overlap Score"])
    sheet_name = "Reaction_Rate_Prediction"
    # 读取现有Excel文件（如果存在）
    if os.path.exists(excel_path):
        with pd.ExcelWriter(excel_path, mode="a", if_sheet_exists="replace", engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        with pd.ExcelWriter(excel_path, mode="w", engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
