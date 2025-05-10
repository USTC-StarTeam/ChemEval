import nltk
from nltk.translate.meteor_score import meteor_score
import sacrebleu
import json
import numpy as np
from rdkit import Chem
import selfies as sf
import math
import os, pandas as pd

def calculate_accuracy(gold_label_list, pred_label_list):
    # 断言保证预测和金标列表长度相同
    assert len(pred_label_list) == len(
        gold_label_list), "The length of prediction and gold label lists must be the same."
    # 初始化指标字典
    num_all = len(pred_label_list)
    metrics = {}
    metrics['num_all'] = num_all
    # 初始化正确和无效的数量
    num_correct = 0
    num_invalid = 0

    # 遍历预测和金标的标签
    for pred_label, gold_label in zip(pred_label_list, gold_label_list):
        # 如果预测标签为空，则无效
        if pred_label is None:
            num_invalid += 1
            continue
        # 尝试将预测和金标的标签转换为字符串
        try:
            pred_label = str(pred_label)
            gold_label = str(gold_label)
        except ValueError:
            num_invalid += 1
            continue
        # 要求一致 答案是Yes or No
        if pred_label == gold_label:
            num_correct += 1

    # 计算准确率
    accuracy = num_correct / num_all if num_all > 0 else 0

    # 添加指标
    metrics['num_correct'] = num_correct
    metrics['num_invalid'] = num_invalid
    metrics['Accuracy'] = accuracy

    return accuracy, metrics


def calculate_accuracy2(gold_label_list, pred_label_list):
    # 断言预测和金标列表的长度相同
    assert len(pred_label_list) == len(
        gold_label_list), "The length of prediction and gold label lists must be the same."
    # 总的样本数
    num_all = len(pred_label_list)
    # 创建一个字典存储指标
    metrics = {}
    metrics['num_all'] = num_all
    # 正确的样本数
    num_correct = 0
    # 无效的样本数
    num_invalid = 0

    # 遍历预测和金标的每个样本
    for pred_label, gold_label in zip(pred_label_list, gold_label_list):
        # 如果预测的标签为空，则无效的样本数加1
        if pred_label is None:
            num_invalid += 1
            continue
        # 将预测和金标的标签都转换为字符串
        try:
            pred_label = str(pred_label)
            gold_label = str(gold_label)
        except ValueError:
            num_invalid += 1
            continue
        # 如果金标标签是预测标签的子集，则正确的样本数加1
        if gold_label in pred_label:
            num_correct += 1

    # 准确率
    accuracy = num_correct / num_all if num_all > 0 else 0

    metrics['num_correct'] = num_correct
    metrics['num_invalid'] = num_invalid
    metrics['Accuracy'] = accuracy

    return accuracy, metrics

# 定义一个函数，将jsonl文件转换为json对象列表
def convert_jsonl_to_json(jsonl_file_path):
    # 创建一个空列表，用于存储json对象
    data_list = []
    # 以读取模式打开jsonl文件，并设置编码为utf-8
    with open(jsonl_file_path, 'r', encoding='utf-8') as file:
        # 遍历文件的每一行
        for line in file:
            # 将每一行转换为json对象
            json_object = json.loads(line.strip())
            # 将json对象添加到列表中
            data_list.append(json_object)
    # 返回列表
    return data_list

# 定义一个函数，从字符串中提取答案
def answer_extract(string):
    # 如果字符串为空，则返回空
    if string is None:
        return None
    # 如果字符串中包含"Yes"和"No"，则返回空
    if "Yes" in string and "No" in string:
        return None
    # 如果字符串中包含"Yes"，则返回"Yes"
    elif "Yes" in string:
        return "Yes"
    # 如果字符串中包含"No"，则返回"No"
    elif "No" in string:
        return "No"
    # 其他情况返回空
    else:
        return None

def list_json_files(folder_path):
    file_list = []
    # 遍历文件夹中的所有文件
    for filename in os.listdir(folder_path):
        # 如果文件以.jsonl结尾
        if filename.endswith(".jsonl"):
            # 将文件名添加到文件列表中
            file_list.append(str(filename))
    # 返回文件列表
    return file_list

def classify(root_path,excel_path):
    # foler_path = r"C:\Users\xiangli67\Downloads\3880数据集\o3-mini\分类"
    foler_path = os.path.join(root_path, "classification")
    file_list = list_json_files(foler_path)
    result_data = []
    # 遍历文件列表
    for name in file_list:
        # 获取文件路径
        path = foler_path + "\\" + name
        # 将文件转换为json格式
        data = convert_jsonl_to_json(path)
        # 初始化金标准列表和预测列表
        golds = []
        pres = []
        # 遍历数据中的每一项
        for i in range(len(data)):
            # 获取金标准
            gold = data[i].get("gold", "").replace("正确", "Yes").replace("错误", "No")
            # 获取预测结果
            pre = data[i].get("answer","")
            if type(pre) == str:
                pre = pre.replace("正确", "Yes").replace("错误", "No")
            # 将预测结果添加到预测列表中
            pres.append(pre)
            # 将金标准添加到金标准列表中
            golds.append(gold)
        # 如果文件名包含"化学文献主题分类"
        if "化学文献主题分类" or "Chemical_Literature_Topic_Classification" in name:
            # 得分计算方式不一样
            result = calculate_accuracy2(golds, pres)
            result_data.append([name, result[1].get("Accuracy")])
        else:
            # 否则，使用默认的得分计算方式
            result = calculate_accuracy(golds, pres)
            result_data.append([name, result[1].get("Accuracy")])


        # 打印文件名和得分
        print(name, f"ACC: ",result[1]['Accuracy'])
        print(name, f"invalid: ",result[1]['num_invalid'])

    df = pd.DataFrame(result_data, columns=["File Name", "ACC Score"])
    sheet_name = "classification"
    # 读取现有Excel文件（如果存在）
    if os.path.exists(excel_path):
        with pd.ExcelWriter(excel_path, mode="a", if_sheet_exists="replace", engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        with pd.ExcelWriter(excel_path, mode="w", engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)