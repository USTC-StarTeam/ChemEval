import sacrebleu
from nltk.translate.bleu_score import corpus_bleu
import numpy as np
import json
from rdkit import Chem
import selfies as sf
import re
import os
import ast
from Levenshtein import distance as lev
from smiles_canonicalization import canonicalize_molecule_smiles, get_molecule_id
from rdkit import Chem, RDLogger
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.Chem import AllChem
from tqdm.auto import tqdm
import logging


# 计算F1分数,注意分子式区分大小写
def calculate_f1_score(true_entities, predicted_entities):
    # 将真实实体和预测实体转换为集合
    true_entities = set(true_entities)
    predicted_entities = set(predicted_entities)
    # 计算真阳性数量
    true_positive = len(true_entities & predicted_entities)
    # 计算精确度
    precision = true_positive / len(predicted_entities) if len(predicted_entities) > 0 else 0
    # 计算召回率
    recall = true_positive / len(true_entities) if len(true_entities) > 0 else 0

    # 计算F1分数
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    return f1_score

# 获取指定文件夹中的json文件列表
def list_json_files(folder_path):
    # 初始化文件列表
    file_list = []
    # 遍历文件夹中的文件
    for filename in os.listdir(folder_path):
        # 如果文件以.jsonl结尾，将其添加到文件列表中
        if filename.endswith(".jsonl"):
            file_list.append(str(filename))
    # 返回文件列表
    return file_list

def add_json_files(file_path,f1_scores):
    jsonl_data=[]
    with open(file_path, "r", encoding="utf-8") as infile:
        idx=0
        for line in infile:
            obj = json.loads(line)
            obj['f1']=f1_scores[idx]
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
            obj['f1']=f1_scores[idx]
            json_data.append(obj)
            idx=idx+1

    with open(file_path.replace('.jsonl','.json'), "w", encoding="utf-8") as outfile:
        for data in json_data:
            json.dump(data, outfile, ensure_ascii=False)
            outfile.write("\n")

def select(fold_path):
    fold_path = os.path.join(fold_path, "reagent_selection")
    logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(fold_path+r'\output.txt','w',encoding='utf-8'),  # 写入文件
                logging.StreamHandler()  # 同时输出到终端
            ],
            force=True
        )
    # 加载json文件
    file_list = list_json_files(fold_path)
    for name in file_list:
        path = fold_path + "\\" + name
        data = []
        # 逐行读取json文件
        with open(path, "r", encoding='utf-8') as file:
            for line in file:
                json_object = json.loads(line.strip())
                data.append(json_object)
        f1_list = []
        for i in range(len(data)):
            correct, flag = 0, 0
            try:
                # ast.literal_eval 用于评估符串中的 Python 字面量表达式
                gold_data = ast.literal_eval(data[i].get("gold", ""))
            except (SyntaxError, ValueError) as e:
                gold_data = data[i].get("gold", "")
            # 如果是列表，则连成字符串
            if type(gold_data) == list:
                gold = '.'.join(gold_data).replace(" ", "").replace("and", ".").replace("\n", "")
            else:
                gold = gold_data.replace(" ", "").replace("and", ".").replace("\n", "")
            try:
                gold1 = sf.decoder(gold)
                # 如果是smiles会转为"" ,不报错
                if gold1 != "":
                    gold = gold1
                    flag = 1
            except:
                gold = gold
            # 统一格式 此时gold必须是smiles selfies必须进行转换
            try:
                mol = Chem.MolFromSmiles(gold)
                gold = Chem.MolToSmiles(mol)
            except:
                #logging.info(f'gold 不能转化为正确SMILES {i+1}')
                #f1_list.append(0)
                break
            # 将字符串分割成列表
            golds = gold.split(".")
            # 金标列表 并且统一了格式
            gold_list = list(set(golds))
            sample_num = len(gold_list)

            # 处理模型预测的答案
            pre_data = data[i].get("answer", "")
            if type(pre_data) == list:
                try:
                    pre = '.'.join(pre_data).replace(" ", "").replace("and", ".")
                except:
                    pre_lists=[]
                    for pre_list in pre_data:
                        for pre in pre_list:
                            pre_lists.append(pre)
                    pre = '.'.join(pre_lists).replace(" ", "").replace("and", ".")
            elif type(pre_data) == str:
                pre = pre_data.replace(" ", "").replace("and", ".")
            else:
                pre = pre_data
            #存在部分结尾为n的输出
            if type(pre) == str and pre.endswith('n'):
                pre = pre[:-1]
            if type(pre) == str and pre.endswith('\n'):
                pre = pre[:-1]
            #存在"[H2O]"的情况 去掉"[" 和"]"
            if type(pre) == str and pre.count('[') == 1 and pre.count(']') == 1 and pre.startswith('[') and pre.endswith(']'):
                pre = pre[1:-1]
            pre_list = []
            if pre is not None:
                for mol in (str(pre).split(".")):
                    try:
                        # 是否需要转为smiles
                        if flag == 1:
                            mol = sf.decoder(mol)
                        pre_s = Chem.MolFromSmiles(mol)
                        pre = Chem.MolToSmiles(pre_s)
                        pre_list.append(pre)
                    except:
                        continue
            # 一行计算一个F1分数
            f1_score = calculate_f1_score(gold_list, pre_list)
            
            f1_list.append(f1_score)
        #print(len(data))
        #print(len(f1_list))
        add_json_files(path,f1_list)
         
        logging.info(f'{name}, F1 score: {np.mean(f1_list)}')

