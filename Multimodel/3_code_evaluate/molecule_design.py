import sacrebleu
from nltk.translate.bleu_score import corpus_bleu
import numpy as np
import json
from rdkit import Chem
import selfies as sf
import re
from Levenshtein import distance as lev
import os
import statistics
import pandas as pd

def re_formers_Lla(string):
    # print(string)
    # 正则表达式，匹配<SMILES>标签内的字符串，并返回
    regex = r"<SMILES>\s*([A-Za-z0-9]+)\s*(?:</SMILES>)?"
    # 使用re.search()函数在string中查找匹配项
    match = re.search(regex, string)
    # 如果找到匹配项，返回第一个捕获组，即<SMILES>标签内的字符串
    if match:
        return match.group(1)
    # 如果没有找到匹配项，返回None
    else:
        return None


def evaluate1(outputs, verbose=False):

    references = []
    hypotheses = []

    for i, (out, gt) in enumerate(outputs):

        if i % 100 == 0:
            if verbose:
                print(i, 'processed.')
        gt_tokens = [c for c in gt]

        out_tokens = [c for c in out]

        references.append([gt_tokens])
        hypotheses.append(out_tokens)

        # mscore = meteor_score([gt], out)
        # meteor_scores.append(mscore)

    # 计算BLEU，调用包
    bleu_score = corpus_bleu(references, hypotheses)
    if verbose: print('BLEU score:', bleu_score)

    # Meteor score
    # _meteor_score = np.mean(meteor_scores)
    # print('Average Meteor score:', _meteor_score)

    rouge_scores = []

    references = []
    hypotheses = []

    levs = []

    num_exact = 0

    bad_mols = 0

    for i, (out, gt) in enumerate(outputs):

        hypotheses.append(out)
        references.append(gt)

        try:
            m_out = Chem.MolFromSmiles(out)
            m_gt = Chem.MolFromSmiles(gt)

            if Chem.MolToInchi(m_out) == Chem.MolToInchi(m_gt): num_exact += 1
            # if gt == out: num_exact += 1 #old version that didn't standardize strings
        except:
            bad_mols += 1

        levs.append(lev(out, gt))

    exact_match_score = num_exact / (i + 1)
    if verbose:
        print('Exact Match:')
        print(exact_match_score)

    levenshtein_score = np.mean(levs)
    if verbose:
        print('Levenshtein:')
        print(levenshtein_score)

    validity_score = 1 - bad_mols / len(outputs)
    if verbose:
        print('validity:', validity_score)

    return bleu_score, exact_match_score, levenshtein_score, validity_score

def evaluate_one_molecule_design(details_df):
    '''
    评估单个分子
    :param details_df: 包含预测值和gold的DataFrame
    :return: 返回评估结果，包括bleu分数、exact_match分数和levenshtein分数
    '''
    bleu_scores = []
    exact_match_scores = []
    levenshtein_scores = []

    # 获取预测值和gold
    tem = list(zip(details_df['answer'], details_df['gold']))

    bleu_score, exact_match_score, levenshtein_score, _ = evaluate1(tem)

    bleu_scores.append(bleu_score)
    exact_match_scores.append(exact_match_score)
    levenshtein_scores.append(levenshtein_score)

    bleu_mean = np.mean(bleu_scores)
    exact_match_mean = np.mean(exact_match_scores)
    levenshtein_mean = np.mean(levenshtein_scores)

    # cal std
    # bleu_variance = statistics.stdev(bleu_scores)
    # exact_match_variance = statistics.stdev(exact_match_scores)
    # levenshtein_variance = statistics.stdev(levenshtein_scores)

    # stds = [bleu_variance, exact_match_variance, levenshtein_variance]

    return [bleu_mean, exact_match_mean, levenshtein_mean]#, stds

# 定义一个函数，用于获取指定文件夹中的所有jsonl文件
def list_json_files(folder_path):
    # 初始化一个空列表，用于存储jsonl文件名
    file_list = []
    # 遍历指定文件夹中的所有文件
    for filename in os.listdir(folder_path):
        # 如果文件名以.jsonl结尾，将其添加到file_list列表中
        if filename.endswith(".jsonl"):
            file_list.append(str(filename))
    # 返回file_list列表
    return file_list

# 如果该模块是作为主模块运行
def design(root_path,excel_path):
    # 指定文件夹路径
    fold_path = os.path.join(root_path, "molecule_design")
    # 调用list_json_files函数，获取指定文件夹中的所有jsonl文件
    file_list = list_json_files(fold_path)
    # 初始化一个空列表，用于存储结果
    res = []
    # 初始化一个空列表，用于存储文件名
    names = []

    result_data = []
    # 遍历file_list列表
    for name in file_list:
        # 拼接文件路径
        path = fold_path + "\\" + name
        # 打印文件名
        print(name)
        # 打开文件
        with open(path, 'r',encoding="utf-8") as file:
            # 将文件内容读取为data列表
            data = [json.loads(line) for line in file]
        # 将data列表转换为dataframe
        details_df = pd.DataFrame(data)
        # 调用evaluate_one_molecule_design函数，计算bleu值
        [bleu_mean, exact_match_mean, levenshtein_mean] = evaluate_one_molecule_design(details_df)
        # 打印bleu值

        result_data.append([name,bleu_mean])
        print(name, "bleu", bleu_mean)
        # 将文件名添加到names列表中
        names.append(name)
        # 将bleu值添加到res列表中
        res.append(bleu_mean)

        
    # 打印结果
    print(names)
    print(res)

    # 创建DataFrame
    df = pd.DataFrame(result_data, columns=["File Name", "bleu"])
    sheet_name = "Molecule_design"
    # 读取现有Excel文件（如果存在）
    if os.path.exists(excel_path):
        with pd.ExcelWriter(excel_path, mode="a", if_sheet_exists="replace", engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        with pd.ExcelWriter(excel_path, mode="w", engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)


