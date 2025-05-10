import os
import nltk
from nltk.translate.meteor_score import meteor_score
import sacrebleu
import json
import numpy as np
from rdkit import Chem
import selfies as sf
import math
import re

# 定义一个函数，用于计算准确率
def calculate_acc(pre_list):
    # 初始化一个空列表，用于存储准确率
    ans_list = []
    acc_sum = 0
    for pre in pre_list:
        acc_sum += (pre-1)/4

    ans_num = acc_sum/len(pre_list)
    # 返回准确率
    return ans_num

# 定义一个函数，将字符串转换为浮点数
def safe_float_convert(value):
    # 如果值为空，返回None
    if value is None:
        return None
    # 尝试将值转换为浮点数，如果转换失败，返回None
    try:
        return float(value)
    except ValueError:
        return None

# 定义一个函数，将jsonl文件转换为json对象列表
def convert_jsonl_to_json(jsonl_file_path):
    # 初始化一个空列表，用于存储json对象
    data_list = []
    # 打开jsonl文件，读取每一行
    with open(jsonl_file_path, 'r', encoding='utf-8') as file:
        for line in file:
            # 将每一行转换为json对象，并添加到列表中
            json_object = json.loads(line.strip())
            data_list.append(json_object)
    # 返回列表
    return data_list

def extract_score(text):
    # 使用正则表达式匹配最终得分，并返回整数值
    match = re.search(r"<最终得分>\D*?([1-5])", text)
    # 如果匹配成功，返回整数值
    if match:
        return int(match.group(1))
    # 否则，使用正则表达式匹配1-5之间的数字，并返回浮点值
    else:
        match = re.search(r'\b[1-5]\b(?=\D*$)', text)
        if match:
            return float(match.group())
        # 否则，返回None
        else:
            return None

# 定义一个函数，用于从文本中提取填空题目最终得分
def extract_decimals(text):
    # 匹配 <最终得分> 标签后面的0到1之间的小数
    match = re.search(r"<最终得分>\D*?([0-1](?:\.\d+)?)", text)
    # 如果匹配成功，返回浮动值
    if match:
        return float(match.group(1))
    # 否则，使用正则表达式匹配0到1之间的其他小数
    else:
        match = re.search(r'\b([0-1](?:\.\d+)?)\b(?=\D*$)', text)
        if match:
            return float(match.group(1))
        # 否则，返回None
        else:
            return None

# 定义一个函数，用于列出给定文件夹中的json文件
def list_json_files(folder_path):
    # 初始化文件列表
    file_list = []
    # 遍历文件夹中的所有文件
    for filename in os.listdir(folder_path):
        # 如果文件以.json结尾，将其添加到文件列表中
        if filename.endswith(".json") or filename.endswith("jsonl"):
            file_list.append(str(filename))
    # 返回文件列表
    return file_list

if __name__ == "__main__":
    # 定义文件夹路径
    ### You need to modify
    foler_path = r"path/to/input/dir"
    excel_path = r"path/to/output/dir"+"other_task_eval.xlsx"

    results = []
    
    sheet_name = "任务"  # 你想写入的 Sheet 名称
    import pandas as pd
    file_list = list_json_files(foler_path) 

    # 遍历文件夹中的所有json文件
    for name in file_list:
        # 定义文件路径
        path = foler_path + "\\" + name
        # 定义空列表，用于存放json文件中的数据
        data = []
        # 打开文件，读取每一行数据
        with open(path, "r", encoding='utf-8') as file:
            for line in file:
                # 将每一行数据转换为json对象
                json_object = json.loads(line.strip())
                # 将json对象添加到data列表中
                data.append(json_object)

        # 定义空列表，用于存放预测值
        golds = []
        # 定义空列表，用于存放真实值
        pres = []
        # 遍历data列表中的每一个json对象
        for i in range(len(data)):
            # 获取预测值
            pre = data[i].get("answer", "")
            if "填空" in name or "FillBlank" in name:
                pre_num = extract_decimals(str(pre))
                # 如果预测值为空，则设置为0
                if pre_num == None:
                    pre_num = 0
                # 将预测值添加到pres列表中
                pres.append(pre_num)
            else:
                pre_num = extract_score(str(pre))
                # 如果预测值为空，则设置为1
                if pre_num == None:
                    pre_num = 1
                # 将预测值添加到pres列表中
                pres.append(pre_num)
        # 计算准确率
        if "填空" in name or "FillBlank" in name:
            ans = sum(pres) / len(pres)
        else:
            ans = calculate_acc(pres)

        # 打印文件名和准确率
        print(name, "Score:",ans)
        results.append([name, ans])
        keywords = ["选择", "填空", "判断", "简答", "计算"]
        filtered_results = [r for r in results if not any(k in r[0] for k in keywords)]
        df = pd.DataFrame(filtered_results, columns=["Model Name", "ACC Score"])
        df.to_excel(excel_path, sheet_name=sheet_name, index=False, engine="openpyxl")
















# import json
# import os
# import openpyxl

# def list_json_files(folder_path):
#     """获取文件夹中所有 JSON 文件的文件名"""
#     return [f for f in os.listdir(folder_path) if f.endswith(".json")]

# def extract_decimals(value):
#     """提取小数，如果提取不到，返回 None"""
#     try:
#         return float(value)
#     except ValueError:
#         return None

# def extract_score(value):
#     """提取分数（假设是整数），如果提取不到，返回 None"""
#     try:
#         return int(value)
#     except ValueError:
#         return None

# def calculate_acc(scores):
#     """计算准确率（示例计算方式，可按需修改）"""
#     return sum(scores) / len(scores) if scores else 0

# if __name__ == "__main__":
#     # 定义 JSON 文件夹路径
#     folder_path = r"C:\Users\xiangli67\Downloads\20250312LLMtest\主观题\task_files"
#     # 指定 Excel 文件路径
#     excel_path = r"C:\Users\xiangli67\Downloads\20250312LLMtest\excel结果\result.xlsx"
#     # 指定 Sheet 名称
#     sheet_name = "简答任务"

#     # 获取 JSON 文件列表
#     file_list = list_json_files(folder_path)

#     # 打开或创建 Excel 文件
#     if os.path.exists(excel_path):
#         wb = openpyxl.load_workbook(excel_path)  # 载入已有文件
#         if sheet_name in wb.sheetnames:
#             ws = wb[sheet_name]  # 选择已有 sheet
#         else:
#             ws = wb.create_sheet(sheet_name)  # 创建新 sheet
#     else:
#         wb = openpyxl.Workbook()  # 创建新 Excel 文件
#         ws = wb.active
#         ws.title = sheet_name  # 设置 sheet 名称

#     # 写入表头
#     if ws.max_row == 1:
#         ws.append(["File Name", "Score"])  # 只写一次表头

#     # 遍历 JSON 文件并计算准确率
#     for name in file_list:

#         if "简答" in name:

#             path = os.path.join(folder_path, name)
#             data = []

#             # 读取 JSON 数据
#             with open(path, "r", encoding='utf-8') as file:
#                 for line in file:
#                     json_object = json.loads(line.strip())
#                     data.append(json_object)

#             # 存放预测值
#             pres = []

#             for item in data:
#                 pre = item.get("answer", "")

#                 if "填空" in name or "FillBlank" in name:
#                     pre_num = extract_decimals(str(pre)) or 0
#                 else:
#                     pre_num = extract_score(str(pre)) or 1
                
#                 pres.append(pre_num)

#             # 计算准确率
#             if "填空" in name or "FillBlank" in name:
#                 ans = sum(pres) / len(pres) if pres else 0
#             else:
#                 ans = calculate_acc(pres)

#             # 写入 Excel
#             ws.append([name, ans])

#     # 保存 Excel 文件
#     wb.save(excel_path)
#     print(f"结果已写入 Excel: {excel_path}")