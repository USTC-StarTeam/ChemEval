import nltk
# nltk.download('punkt')
# nltk.download('punkt_tab')
from nltk.translate.bleu_score import sentence_bleu
from nltk.tokenize import word_tokenize
from nltk.translate.meteor_score import meteor_score
import json
from nltk.translate.bleu_score import corpus_bleu
from rdkit import Chem
import selfies as sf
import re
import os
from molecule_design_S import smiles_to_fps, calculate_tanimoto_similarity
from py2opsin import py2opsin
import logging
import Levenshtein


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
            '''将json对象添加到列表中'''
    return data_list


def calculate_exact_match(references, predictions):
    # 检查参考答案和预测答案的数量是否相等
    assert len(references) == len(predictions), "参考答案和预测答案的数量必须相等"
    # 计算总的样本数
    total = len(references)
    # 初始化精确匹配的数量
    exact_match_count = 0

    exact_matches = []
    # 遍历参考答案和预测答案
    ### modified
    # for ref_list, pred_list in zip(references, predictions):
    for i in range(len(references)):

        # 将参考答案和预测答案都转换为小写 IUPAC不区分大小写
        # ref_str = ref_list[0].lower()
        # pred_str = pred_list[0].lower()
        ref_str = references[i].lower()
        pred_str = predictions[i].lower()
        # 统计分数
        if ref_str == pred_str:
            print()
            # 如果参考答案和预测答案相同，则将精确匹配的数量加1
            exact_match_count += 1
            exact_matches.append(1)
        else:
            exact_matches.append(0)

    # 计算精确匹配的分数
    exact_match_score = (exact_match_count / total)
    # 返回精确匹配的分数
    return exact_match_score, exact_matches


def calculate_edit_distance(references, predictions):
    edit_distances = []
    assert len(references) == len(predictions), "参考答案和预测答案的数量必须相等"
    for i in range(len(references)):
        # 分词
        ref_tokens = word_tokenize(references[i].lower())
        hyp_tokens = word_tokenize(predictions[i].lower())

        # 编辑距离
        edit_distance = Levenshtein.distance(references[i], predictions[i])
        edit_distances.append(edit_distance)

    return edit_distances


def calculate_bleu(references, predictions):
    bleu_scores = []
    assert len(references) == len(predictions), "参考答案和预测答案的数量必须相等"
    for i in range(len(references)):
        # 分词
        ref_tokens = word_tokenize(references[i].lower())
        hyp_tokens = word_tokenize(predictions[i].lower())

        # BLEU分数（使用1-gram到4-gram）
        bleu_score = sentence_bleu([ref_tokens], hyp_tokens, weights=(0.25, 0.25, 0.25, 0.25))
        bleu_scores.append(bleu_score)

    return bleu_scores


def iupac_to_smiles(preds_iupac_list):
    preds_smi_list = []
    for idx, iupac in enumerate(preds_iupac_list):
        iupac = iupac.replace(';', ' ')
        try:
            smiles_string = py2opsin(
                chemical_name=iupac,
                output_format="SMILES",
                jar_fpath="./opsin-cli-2.8.0-jar-with-dependencies.jar",
            )
            if smiles_string:
                preds_smi_list.append(smiles_string)
            else:
                preds_smi_list.append(None)
                logging.info(f"[Warning] No compound found for pre IUPAC name: {iupac}, id: {idx + 1}")
        except Exception as e:
            preds_smi_list.append(None)
            logging.info(f"[Error] Failed to fetch SMILES for pre {iupac}: {e}, id: {idx + 1}")

    return preds_smi_list


# 定义一个函数，用于从字符串中提取LlaSMOl IUPAC字符串
def re_formers_Lla(string):
    # logging.info(string)
    # 定义一个正则表达式，用于匹配IUPAC字符串
    regex = r"<IUPAC>\s*([^<]+?)\s*</IUPAC>"

    # 使用re.search()函数查找匹配项
    match = re.search(regex, string)
    # 如果找到匹配项，返回匹配的字符串
    if match:
        return match.group(1)
    # 如果没有找到匹配项，返回None
    else:
        return None


# 定义一个函数，用于从字符串中提取ChemDFM字符串
def re_formers_DFM(string):
    # 定义一个正则表达式，用于匹配DFM字符串
    regex = r'"answer":\s*"([A-Za-z0-9]+)"'
    # 使用re.search()函数查找匹配项
    match = re.search(regex, string)
    # 如果找到匹配项，返回匹配的字符串
    if match:
        return match.group(1)
    # 如果没有找到匹配项，返回原始字符串
    else:
        return string


# 定义一个函数，用于列出指定文件夹中的json文件
def list_json_files(folder_path):
    # 初始化文件列表
    file_list = []
    # 遍历指定文件夹中的所有文件
    for filename in os.listdir(folder_path):
        # 如果文件以.jsonl结尾，将其添加到文件列表中
        ### modified
        if filename.endswith(".jsonl"):
            file_list.append(str(filename))
    # 返回文件列表
    return file_list


def add_json_files(file_path, sims, exact_matches, bleus, edit_distances):
    jsonl_data = []
    with open(file_path, "r", encoding="utf-8") as infile:
        idx = 0
        for line in infile:
            obj = json.loads(line)
            obj['tanimoto'] = sims[idx]
            obj['exact_match'] = exact_matches[idx]
            obj['bleu'] = bleus[idx]
            obj['edit_distance'] = edit_distances[idx]
            jsonl_data.append(obj)
            idx = idx + 1

    with open(file_path, "w", encoding="utf-8") as outfile:
        for data in jsonl_data:
            json.dump(data, outfile, ensure_ascii=False)
            outfile.write("\n")

    json_data = []
    with open(file_path.replace('.jsonl', '.json'), "r", encoding="utf-8") as infile:
        idx = 0
        for line in infile:
            obj = json.loads(line)
            obj['tanimoto'] = sims[idx]
            obj['exact_match'] = exact_matches[idx]
            obj['bleu'] = bleus[idx]
            obj['edit_distance'] = edit_distances[idx]
            json_data.append(obj)
            idx = idx + 1

    with open(file_path.replace('.jsonl', '.json'), "w", encoding="utf-8") as outfile:
        for data in json_data:
            json.dump(data, outfile, ensure_ascii=False)
            outfile.write("\n")


def design_I(fold_path):
    # 指定文件夹路径
    fold_path = os.path.join(fold_path, "molecule_design_I")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(fold_path + r'\output.txt', 'w', encoding='utf-8'),  # 写入文件
            logging.StreamHandler()  # 同时输出到终端
        ],
        force=True
    )
    # 获取文件夹中的json文件列表
    file_list = list_json_files(fold_path)
    # 遍历文件列表
    for name in file_list:
        # 获取文件路径
        path = fold_path + "\\" + name
        # 将json文件转换为字典
        data = convert_jsonl_to_json(path)

        # 初始化gold和pre列表
        golds = []
        pres = []
        gold_smis = []

        # 遍历字典中的数据
        for i in range(len(data)):
            # 获取gold
            gold = data[i].get("gold", "")
            # 获取pre
            pre = data[i].get("answer", "")
            # 获取gold_smi
            gold_smi = data[i].get("gold_smi", "")
            # 对pre进行处理
            if pre:
                pre_f = re_formers_DFM(pre)
            else:
                pre_f = ""
            # 将处理后的pre添加到pres列表中
            pres.append(pre_f)
            # 将gold添加到golds列表中
            golds.append(gold)
            # 将gold_smi添加到gold_smis列表中
            gold_smis.append(gold_smi)

        # 初始化metrics列表
        metrics = ["tanimoto", "exact match", 'bleu', 'edit distance']

        # 计算分数
        metric_results = {}
        # metric_results = calculate_exact_match(golds, pres)
        preds_smi_list = iupac_to_smiles(pres)
        for metric in metrics:
            if 'tanimoto' == metric:
                print("pre",preds_smi_list)
                print("gold",gold_smis)
                sims = calculate_tanimoto_similarity(preds_smi_list, gold_smis)
                metric_results[metric] = sum(sims) / len(sims)
            elif 'exact match' == metric:
                metric_results[metric], exact_matches = calculate_exact_match(golds, pres)
            elif 'bleu' == metric:
                bleus = calculate_bleu(golds, pres)
                metric_results[metric] = sum(bleus) / len(bleus)
            elif 'edit distance' == metric:
                edit_distances = calculate_edit_distance(golds, pres)
                metric_results[metric] = sum(edit_distances) / len(edit_distances)

            # 打印结果
            logging.info(f'{name}, {metric}: {metric_results[metric]}')

        add_json_files(path, sims, exact_matches, bleus, edit_distances)

        logging.info("-" * 100)

    # em_score = calculate_exact_match(references, candidates)
    # logging.info(f"ACC: {em_score:.2f}%")

design_I(r"D:\project\nips-chemeval\2chemeval_code-2\output_gpt-2t\output2_me")