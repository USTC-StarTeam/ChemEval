import os
import json
import ast

def convert_jsonl_to_json(jsonl_file_path):
    '''将jsonl文件转换为json对象列表'''
    data_list = []
    with open(jsonl_file_path, 'r', encoding='utf-8') as file:
        '''打开jsonl文件'''
        for line in file:
            '''遍历文件每一行'''
            json_object = json.loads(line.strip())
            '''将每一行转换为json对象'''
            data_list.append(json_object)
        '''将所有json对象添加到列表中'''
    return data_list


def calculate_accuracy(jsonl_file):
    '''
    计算给定jsonl文件中问题的正确率
    :param jsonl_file: jsonl文件路径
    :return: 问题正确率
    '''
    total_correct = 0
    total_count = 0

    with open(jsonl_file, 'r', encoding='utf-8') as file:
        for index,line in enumerate(file):
            # print(index+1)
            data = json.loads(line)
            if data["answer"] is None:
                gold_labels = ast.literal_eval(data["gold"])
                answer_labels = ['X'] * len(gold_labels)
            else:
                answer_labels = data["answer"]
                if type(data["gold"]) == str:
                    gold_labels = ast.literal_eval(data["gold"])
                elif type(data["gold"]) == list:
                    gold_labels = " ".join(data["gold"])
                # elif type(data["gold"]) == list:
                #     gold_labels =ast.literal_eval(data["gold"][0])

            correct_predictions = sum(1 for x, y in zip(answer_labels, gold_labels) if x == y)
            total_predictions = len(gold_labels)

            total_correct += correct_predictions
            total_count += total_predictions
    # print(total_correct)
    overall_accuracy = total_correct / total_count if total_count > 0 else 0
    return overall_accuracy

# 定义一个函数，用于获取指定文件夹中的所有jsonl文件
def list_json_files(folder_path):
    # 初始化一个空列表，用于存储jsonl文件名
    file_list = []
    # 遍历指定文件夹中的所有文件
    for filename in os.listdir(folder_path):
        # 如果文件名以.jsonl结尾
        if filename.endswith(".jsonl"):
            # 将文件名添加到列表中
            file_list.append(str(filename))
    # 返回文件列表
    return file_list
import pandas as pd
# 如果该模块是主模块
def recognize(fold_path,excel_path):
    # 指定文件夹路径
    fold_path = os.path.join(fold_path,"entity_recognition")
    # 调用函数，获取指定文件夹中的所有jsonl文件
    file_list = list_json_files(fold_path)
    # 遍历文件列表
    result_data  =[]
    for name in file_list:
        # 指定文件路径
        path = fold_path + "\\" + name
        # 调用函数，计算准确率
        accuracy = calculate_accuracy(path)
        # 打印文件名和准确率
        result_data.append([name,accuracy])
        print(name, accuracy)
    # 创建DataFrame
    df = pd.DataFrame(result_data, columns=["File Name", "ACC Score"])
    sheet_name = "Entity_recognition"
    
    # 读取现有Excel文件（如果存在）
    if os.path.exists(excel_path):
        with pd.ExcelWriter(excel_path, mode="a", if_sheet_exists="replace", engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        with pd.ExcelWriter(excel_path, mode="w", engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
