import json
import re
import os

def calculate_f1_score(true_entities, predicted_entities):
    # 定义正则表达式，用于匹配括号内的内容
    pattern = r'\(.*?\)'
    # 将true_entities中的空格替换为空字符串，然后使用正则表达式匹配括号内的内容
    true_entities = true_entities.replace(" ", "")
    true_entities = re.findall(pattern, true_entities)
    # 将predicted_entities中的空格替换为空字符串，然后使用正则表达式匹配括号内的内容
    predicted_entities = predicted_entities.replace(" ", "")
    predicted_entities = re.findall(pattern, predicted_entities)

    # 将true_entities转换为集合
    true_entities = set(true_entities)
    # 将predicted_entities转换为集合
    predicted_entities = set(predicted_entities)
    # 计算真阳性数量
    true_positive = len(true_entities & predicted_entities)
    # 计算精确度
    precision = true_positive / len(predicted_entities) if len(predicted_entities) > 0 else 0
    # 计算召回率
    recall = true_positive / len(true_entities) if len(true_entities) > 0 else 0

    # 计算f1分数
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    return f1_score


def calculate_accuracy(golds, pres):
    '''
    计算F1分数
    :param golds: 真实输出
    :param pres: 预测输出
    :return: F1分数
    '''
    correct_count = 0
    total_count = len(golds)
    for i in range(len(golds)):
        true_output = golds[i].lower()
        if type(pres[i]) == list:
            pres[i] = ' '.join(str(pres[i]))
        my_output = pres[i].lower()
        # print(true_output, 11111, my_output)
        f1_score = calculate_f1_score(true_output, my_output)
        # print(f1_score)
        correct_count += f1_score

    return correct_count / total_count


def list_json_files(folder_path):
    file_list = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".jsonl"):
            file_list.append(str(filename))
    return file_list
import pandas as pd
def relation_extract(fold_path, excel_path):
    # 定义文件路径
    # fold_path = r"C:\Users\xiangli67\Downloads\3880数据集\o3-mini\关系抽取"
    fold_path = os.path.join(fold_path, "relation_extraction")
    # 获取文件列表
    file_list = list_json_files(fold_path)
    # 遍历文件列表
    result_data = []
    for name in file_list:
        # 定义文件路径
        path = fold_path + "\\" + name

        # 定义空列表，用于存放数据
        data = []
        # 读取文件内容
        with open(path, "r", encoding='utf-8') as file:
            for line in file:
                # 将每行内容转换为json对象
                json_object = json.loads(line.strip())
                # 将json对象添加到数据列表中
                data.append(json_object)

        # 定义空列表，用于存放金标
        golds = []
        # 定义空列表，用于存放预测值
        pres = []
        # 遍历数据列表
        for i in range(len(data)):
            # 获取金标值
            gold = data[i].get("gold", "")
            # 如果预测值是字符串类型，则将其添加到预测值列表中
            if type(data[i]) == str:
                pre = data[i]
            else:
                # 获取预测值
                pre = data[i].get("answer", "")
            # 如果预测值是字符串类型，且不为空，则进行格式化处理
            if type(pre) == str and '(' not in pre and len(pre) != 0:
                # 定义空列表，用于存放格式化后的预测值
                formatted_pairs = []
                # 将预测值按逗号分隔开，然后两两一组，用括号包围，然后用逗号连接
                pre_list = pre.replace(";", ",").split(",")
                for i in range(0, len(pre_list), 2):
                    if i + 1 < len(pre_list):
                        formatted_pairs.append(f"({pre_list[i]}, {pre_list[i + 1]})")
                # 将格式化后的预测值用逗号连接
                pre = ','.join(formatted_pairs)
            # 将预测值添加到预测值列表中
            pres.append(pre)
            # 将金标值添加到金标列表中
            golds.append(gold)
        # 计算准确率
        accuracy = calculate_accuracy(golds, pres)
        # 打印文件名和准确率
        result_data.append([name,round(accuracy,4)])
        print(name, f"F1: {accuracy:.4f}")

    # 创建DataFrame
    df = pd.DataFrame(result_data, columns=["File Name", "ACC Score"])
    sheet_name = "Relation_extraction"
    # 读取现有Excel文件（如果存在）
    if os.path.exists(excel_path):
        with pd.ExcelWriter(excel_path, mode="a", if_sheet_exists="replace", engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        with pd.ExcelWriter(excel_path, mode="w", engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
