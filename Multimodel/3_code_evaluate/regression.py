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
import logging


def calculate_rmse(gold_text_list, pred_text_list):
    # 断言pred_text_list的长度等于gold_text_list的长度
    assert len(pred_text_list) == len(gold_text_list)

    # 初始化metrics字典
    num_all = len(pred_text_list)
    metrics = {}
    metrics['num_all'] = num_all
    # 初始化num_no_answer和num_invalid
    num_no_answer = 0
    num_invalid = 0
    # 初始化new_pred_text_list和new_gold_text_list
    new_pred_text_list, new_gold_text_list = [], []
    # 遍历pred_text_list和gold_text_list
    for (pred_item, gold_item) in zip(pred_text_list, gold_text_list):
        # 如果pred_item为None，则num_no_answer加1
        if pred_item is None:
            num_no_answer += 1
            continue
        # 如果pred_item为空字符串，则num_no_answer加1
        if pred_item == '':
            num_no_answer += 1
            continue
        # 将pred_item转换为浮点数
        try:
            pred_item = float(pred_item)
        except (SyntaxError, ValueError):
            # 如果转换失败，则num_invalid加1
            num_invalid += 1
            continue
        # 将gold_item转换为浮点数
        gold_item = float(gold_item)
        # 将pred_item和gold_item添加到new_pred_text_list和new_gold_text_list中
        new_pred_text_list.append(pred_item)
        new_gold_text_list.append(gold_item)

    # 将new_pred_text_list和new_gold_text_list转换为numpy数组
    new_pred_text_list = np.array(new_pred_text_list)
    new_gold_text_list = np.array(new_gold_text_list)
    # 计算方差
    score = np.sqrt(((new_pred_text_list - new_gold_text_list) ** 2).mean())

    # 更新metrics字典
    metrics['num_no_answer'] = num_no_answer
    metrics['num_invalid'] = num_invalid
    metrics['RMSE'] = score

    return score, num_no_answer
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
    # 打开jsonl文件，读取每一行，将每一行转换为json对象，添加到data_list中
    with open(jsonl_file_path, 'r', encoding='utf-8') as file:
        for line in file:
            json_object = json.loads(line.strip())
            data_list.append(json_object)
    # 返回data_list
    return data_list

# 定义一个函数，从输入字符串中提取数字
def extract_number(input_string):
    # 使用正则表达式在输入字符串中查找数字，并返回第一个匹配的数字
    match = re.search(r'[-+]?\d*\.?\d+(?=\D*$)', input_string)
    if match:
        return float(match.group())
    else:
        # 如果找不到数字，返回None
        return None

def list_json_files(folder_path):
    '''获取指定文件夹中所有json文件的文件名'''
    file_list = []
    # 遍历文件夹中的所有文件
    for filename in os.listdir(folder_path):
        # 如果文件名以.json结尾
        if filename.endswith(".jsonl"):
            # 将文件名添加到文件列表中
            file_list.append(str(filename))
    # 返回文件列表
    return file_list


def compute_rae(gold, pre):
    # 转换为numpy数组，以便进行计算
    gold = np.array(gold)
    pre = np.array(pre)

    # 初始化计数器
    num_valid = 0

    # 初始化新的金标准和预测列表（忽略无效项）
    new_gold = []
    new_pre = []

    # 遍历gold和pre
    for g, p in zip(gold, pre):
        # 如果p是None或空字符串，跳过这个预测
        if p is None or p == '':
            continue
        try:
            p = float(p)  # 尝试将预测值转换为浮动数
        except (ValueError, TypeError):
            continue  # 如果转换失败，则忽略这个预测
        try:
            g = float(g)  # 将实际值转换为浮动数
        except (ValueError, TypeError):
            continue  # 如果实际值转换失败，也跳过

        # 只有在预处理有效时，才加入有效的列表
        new_gold.append(g)
        new_pre.append(p)
        num_valid += 1

    if num_valid == 0:
        return None, num_valid  # 如果没有有效的值，返回None

    # 计算基准模型（均值模型）预测的误差
    mean_gold = np.mean(new_gold)
    baseline_error = np.sum(np.abs(np.array(new_gold) - mean_gold))

    # 计算模型的绝对误差
    model_error = np.sum(np.abs(np.array(new_gold) - np.array(new_pre)))

    # 计算RAE
    rae = model_error / baseline_error

    return rae, num_valid

def compute_r2(gold, pre):
    # 转换为numpy数组，以便进行计算
    gold = np.array(gold)
    pre = np.array(pre)

    # 初始化有效预测值的数量
    num_valid = 0

    # 初始化新的金标准和预测列表（忽略无效项）
    new_gold = []
    new_pre = []

    # 遍历gold和pre
    for g, p in zip(gold, pre):
        # 如果p是None或空字符串，跳过这个预测
        if p is None or p == '':
            continue
        try:
            p = float(p)  # 尝试将预测值转换为浮动数
        except (ValueError, TypeError):
            continue  # 如果转换失败，则忽略这个预测
        try:
            g = float(g)  # 将实际值转换为浮动数
        except (ValueError, TypeError):
            continue  # 如果实际值转换失败，也跳过

        # 只有在预处理有效时，才加入有效的列表
        new_gold.append(g)
        new_pre.append(p)
        num_valid += 1

    if num_valid == 0:
        return None, num_valid  # 如果没有有效的值，返回None

    # 计算总平方和 (TSS)
    mean_gold = np.mean(new_gold)
    total_sum_of_squares = np.sum((new_gold - mean_gold) ** 2)

    new_pre = np.array(new_pre)
    # 计算残差平方和 (RSS)
    residual_sum_of_squares = np.sum((new_gold - new_pre) ** 2)

    # 计算R²
    r2 = 1 - (residual_sum_of_squares / total_sum_of_squares)

    return r2, num_valid

def compute_nrmse(gold, pre):
    # 转换为numpy数组，以便进行计算
    gold = np.array(gold)
    pre = np.array(pre)

    # 初始化有效预测值的数量
    num_valid = 0

    # 初始化新的金标准和预测列表（忽略无效项）
    new_gold = []
    new_pre = []

    # 遍历gold和pre
    for g, p in zip(gold, pre):
        # 如果p是None或空字符串，跳过这个预测
        if p is None or p == '':
            continue
        try:
            p = float(p)  # 尝试将预测值转换为浮动数
        except (ValueError, TypeError):
            continue  # 如果转换失败，则忽略这个预测
        try:
            g = float(g)  # 将实际值转换为浮动数
        except (ValueError, TypeError):
            continue  # 如果实际值转换失败，也跳过

        # 只有在预处理有效时，才加入有效的列表
        new_gold.append(g)
        new_pre.append(p)
        num_valid += 1

    if num_valid == 0:
        return None, num_valid  # 如果没有有效的值，返回None

    # 计算RMSE
    rmse = np.sqrt(np.mean((np.array(new_gold) - np.array(new_pre)) ** 2))

    # 计算金标准的最大值和最小值
    gold_max = np.max(new_gold)
    gold_min = np.min(new_gold)

    # 计算NRMSE
    nrmse = rmse / np.mean(np.abs(gold))  # 或 std(gold)

    return nrmse, num_valid

def add_json_files(file_path,rmse_list,rae_list,r2_list,nrmse_list):
    jsonl_data=[]
    with open(file_path, "r", encoding="utf-8") as infile:
        idx=0
        for line in infile:
            obj = json.loads(line)
            obj['RMSE']=rmse_list[idx]
            obj['RAE']=rae_list[idx]
            obj['R2']=r2_list[idx]
            obj['NRMSE']=nrmse_list[idx]
            jsonl_data.append(obj)
            idx=idx+1

    with open(file_path, "w", encoding="utf-8") as outfile:
        for data in jsonl_data:
            json.dump(data, outfile, ensure_ascii=False)
            outfile.write("\n")

    json_data=[]
    with open(file_path.replace('.jsonl','.json'), "r", encoding="utf-8") as infile:
        idx=0
        for line in infile:
            obj = json.loads(line)
            obj['RMSE']=rmse_list[idx]
            obj['RAE']=rae_list[idx]
            obj['R2']=r2_list[idx]
            obj['NRMSE']=nrmse_list[idx]
            json_data.append(obj)
            idx=idx+1

    with open(file_path.replace('.jsonl','.json'), "w", encoding="utf-8") as outfile:
        for data in json_data:
            json.dump(data, outfile, ensure_ascii=False)
            outfile.write("\n")

def reg(root_path):
    # 定义文件路径
    foler_path = os.path.join(root_path, "regression")
    logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(foler_path+r'\output.txt','w',encoding='utf-8'),  # 写入文件
                logging.StreamHandler()  # 同时输出到终端
            ],
            force=True
        )
    # 获取文件列表
    file_list = list_json_files(foler_path)
    # 遍历文件列表
    for name in file_list:
        # 定义文件路径
        path = foler_path + "\\" + name
        # 定义空列表
        data = []
        # 打开文件
        with open(path, "r", encoding='utf-8') as file:
            # 遍历文件内容
            for line in file:
                # 将每一行内容加载为json对象
                json_object = json.loads(line.strip())
                # 将json对象添加到列表中
                data.append(json_object)

        # 定义空列表
        golds = []
        pres = []
        # 遍历数据列表
        for i in range(len(data)):
            # 打印当前索引
            # logging.info(i)
            # 获取标准答案
            gold = data[i].get("gold", "")
            # 获取预测答案
            if type(data[i]) == str:
                pre = data[i]
            else:
                pre = data[i].get("answer", "")
            # 清洗不佳时使用
            pre_num = extract_number(str(pre))
            gold_num = extract_number(str(gold))
            # 如果预测答案不是数字，则将其设置为None
            if safe_float_convert(pre_num) is None or safe_float_convert(pre_num) > 1000 or safe_float_convert(pre_num) < -1000:
                num = None
            else:
                # 将预测答案转换为数字
                num = safe_float_convert(pre_num)
            # 将预测答案添加到列表中
            pres.append(num)
            # 将标准答案添加到列表中
            golds.append(gold_num)
        # 打印列表长度和无答案数量
        # logging.info(len(pres), pres)
        # 计算RMSE分数
        rmse_result, no_answer = calculate_rmse(golds, pres)
        rae_reult,rae_no_answer = compute_rae(golds,pres)
        r2_result,r2_no_answer = compute_r2(golds,pres)
        nrmse_result,nrmse_result_no_answer = compute_nrmse(golds,pres)
        # 打印结果
        logging.info(f'{name}, RMSE Score: {float(rmse_result), no_answer}')
        logging.info(f'{name}, RAE Score: {safe_float_convert(rae_reult), no_answer}')
        logging.info(f'{name}, R2 Score: {safe_float_convert(r2_result),no_answer}')
        logging.info(f'{name}, NRMSE Score: {safe_float_convert(nrmse_result), no_answer}')
        logging.info('-'*100)
reg("D:\project\\nips-chemeval\\pipline\\data\\机评结果\\按t拆分\\gemini2.5pro-4-2t\\output2_me\\")