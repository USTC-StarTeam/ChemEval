import numpy as np
from tqdm.auto import tqdm
from smiles_canonicalization import canonicalize_molecule_smiles, get_molecule_id
import json
import re
import os
import selfies as sf
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import numpy as np
import pandas as pd
import logging

def is_too_long_or_single_smis(smiles: str, min_length: int = 4908) -> bool:
    return len(smiles) >= min_length or len(set(smiles)) == 1

# 精确匹配计算函数
def judge_exact_match(pred_can_smiles_list, gold_can_smiles_list):
    # 断言pred_can_smiles_list和gold_can_smiles_list的长度相等
    assert len(pred_can_smiles_list) == len(gold_can_smiles_list)
    # 初始化精确匹配标签列表
    exact_match_labels = []
    # 遍历pred_can_smiles_list和gold_can_smiles_list
    for pred_smiles, gold_smiles_list in zip(pred_can_smiles_list, gold_can_smiles_list):
        # 如果pred_smiles为None，则将精确匹配标签列表添加为False
        if pred_smiles is None:
            exact_match_labels.append(False)
            continue
        # 将pred_smiles转换为inchi
        pred_smiles_inchi = get_molecule_id(pred_smiles)
        # 初始化样本精确匹配为False
        sample_exact_match = False
        # 遍历gold_smiles_list
        for gold_smiles in gold_smiles_list:
            # 断言gold_smiles不为None
            assert gold_smiles is not None
            # 将gold_smiles转换为inchi
            gold_smiles_inchi = get_molecule_id(gold_smiles)
            # 判断inchi是否相等
            if pred_smiles_inchi == gold_smiles_inchi:
                # 如果相等，则将样本精确匹配改为True
                sample_exact_match = True
                break
        # 将样本精确匹配添加到精确匹配标签列表
        exact_match_labels.append(sample_exact_match)
    # 返回精确匹配标签列表
    return np.array(exact_match_labels)

# 将 SMILES 转换为 Morgan 指纹
def smiles_to_fps(preds_smiles_list, golds_smiles_list, radius=2, nBits=2048):
    preds_fps_list = []
    preds_mols_list = []
    golds_fps_list = []
    golds_mols_list = []
    
    for idx,smi in enumerate(preds_smiles_list):
        #logging.info(smi)
        if smi is None:
            logging.info(f"pre SMILES无效: {idx+1}")
            preds_fps_list.append(None)
            preds_mols_list.append(None)
            continue
        
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            logging.info(f"pre 无法解析 SMILES: {smi}")
            preds_fps_list.append(None)
            preds_mols_list.append(None)
        else:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
            preds_fps_list.append(fp)
            preds_mols_list.append(mol)

    for smi in golds_smiles_list:
        #logging.info(smi)
        if smi is None:
            logging.info(f"gold SMILES空: {idx+1}")
            golds_fps_list.append(None)
            golds_mols_list.append(None)
            continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            logging.info(f"gold 无法解析 SMILES: {smi}")
            golds_fps_list.append(None)
            golds_mols_list.append(None)
        else:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
            golds_fps_list.append(fp)
            golds_mols_list.append(mol)
    return preds_fps_list, preds_mols_list, golds_fps_list, golds_mols_list

# 计算指纹相似度（Tanimoto）
def calculate_tanimoto_similarity(preds_smiles_list,golds_smiles_list):
    preds_fps_list,preds_mols_list,golds_fps_list,golds_mols_list=smiles_to_fps(preds_smiles_list,golds_smiles_list)
    sims=[]
    for i in range(len(preds_fps_list)):
        if preds_fps_list[i] is not None and golds_fps_list[i] is not None:
            sim = DataStructs.TanimotoSimilarity(preds_fps_list[i], golds_fps_list[i])
        else:
            sim = 0
        sims.append(sim)
    return sims
    mol1 = Chem.MolFromSmiles(pred_smiles) if pred_smiles else None
    mol2 = Chem.MolFromSmiles(gold_smiles) if gold_smiles else None
    if mol1 is None or mol2 is None:
        return None
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def calculate_smiles_metrics(
        # 列表
        preds_smiles_list,
        golds_smiles_list,
        metrics=('exact_match')
):
    # 计算总的翻译个数
    num_all = len(preds_smiles_list)
    # 断言保证翻译个数不为0
    assert num_all > 0
    # 断言保证翻译个数与翻译结果个数相同
    assert num_all == len(golds_smiles_list)
    # 翻译的个数，可能为列表
    k = len(preds_smiles_list[0])

    # 创建字典，用于存储每个翻译句子的指标
    dk_pred_smiles_list_dict = {}
    dk_pred_no_answer_labels_dict = {}
    dk_pred_invalid_labels_dict = {}
    for dk in range(k):
        dk_pred_smiles_list_dict[dk] = []
        dk_pred_no_answer_labels_dict[dk] = []
        dk_pred_invalid_labels_dict[dk] = []
    # 标记预测的情况
    for idx,pred_smiles_list in enumerate(tqdm(preds_smiles_list)):
        ### modified
        #if pred_smiles_list is None:
        
        if pred_smiles_list =='':
            for dk in range(k):
                dk_pred_no_answer_labels_dict[dk].append(True)
                dk_pred_invalid_labels_dict[dk].append(False)
                dk_pred_smiles_list_dict[dk].append(None)
            continue
        elif pred_smiles_list is None:
            for dk in range(k):
                dk_pred_no_answer_labels_dict[dk].append(False)
                dk_pred_invalid_labels_dict[dk].append(True)
                dk_pred_smiles_list_dict[dk].append(None)
            continue
        assert len(pred_smiles_list) == k
        for dk, item in enumerate(pred_smiles_list):
            # item = item.strip()
            ### modified
            #if item == '' or item is None:
            if item == '':
                item = None
                dk_pred_no_answer_labels_dict[dk].append(True)
                dk_pred_invalid_labels_dict[dk].append(False)
            elif item==None:
                dk_pred_no_answer_labels_dict[dk].append(False)
                dk_pred_invalid_labels_dict[dk].append(True)
            elif is_too_long_or_single_smis(item):
                item = None
                dk_pred_no_answer_labels_dict[dk].append(False)
                dk_pred_invalid_labels_dict[dk].append(True)
            else:
                # 是否存在答案
                dk_pred_no_answer_labels_dict[dk].append(False)
                item = canonicalize_molecule_smiles(item)
                # print(item)
                if item is None:
                    dk_pred_invalid_labels_dict[dk].append(True)
                else:
                    dk_pred_invalid_labels_dict[dk].append(False)
                #print(dk_pred_invalid_labels_dict[dk][-1])
            dk_pred_smiles_list_dict[dk].append(item)
            

    # 处理golds_smiles_list
    new_list = []
    for gold_smiles_list in tqdm(golds_smiles_list):
        sample_gold_smiles_list = []
        for gold in gold_smiles_list:
            item = gold.strip()
            new_item = canonicalize_molecule_smiles(item, return_none_for_error=False)
            # if new_item is None:
            #     new_item = item  #TODO
            # assert new_item is not None, item
            #sample_gold_smiles_list.append(new_item)
            ### modified
        #new_list.append(sample_gold_smiles_list)
        new_list.append(new_item)
    golds_smiles_list = new_list

    # 初始化metric_results字典，用于存储指标结果
    metric_results = {'num_all': num_all}

    # 初始化tk_pred_no_answer_labels数组，用于存储预测为无答案的标签
    tk_pred_no_answer_labels = np.array([True] * num_all)
    # 初始化tk_pred_invalid_labels数组，用于存储预测为无效的标签
    tk_pred_invalid_labels = np.array([True] * num_all)
    # 遍历每个翻译器的预测结果
    for dk in range(k):
        # 获取当前翻译器的无答案标签
        dk_no_answer_labels = dk_pred_no_answer_labels_dict[dk]
        # 获取当前翻译器的无效标签
        dk_invalid_labels = dk_pred_invalid_labels_dict[dk]
        # 与当前翻译器的无答案标签进行与运算
        tk_pred_no_answer_labels = tk_pred_no_answer_labels & dk_no_answer_labels
        # 与当前翻译器的无效标签进行与运算
        tk_pred_invalid_labels = tk_pred_invalid_labels & dk_invalid_labels
        # 统计当前翻译器的无答案标签数量
        metric_results['num_t%d_no_answer' % (dk + 1)] = tk_pred_no_answer_labels.sum().item()
        # 统计当前翻译器的无效标签数量
        metric_results['num_t%d_invalid' % (dk + 1)] = tk_pred_invalid_labels.sum().item()

    # d1_no_answer_labels = dk_pred_no_answer_labels_dict[0]
    # # logging.info(np.array(d1_no_answer_labels).sum().item())
    # for label, item in zip(d1_no_answer_labels, preds_smiles_list):
    #     if label:
    #         logging.info(item)
    # 计算精确匹配分数
    for metric in metrics:
        if metric == 'exact_match':
            # 初始化tk_exact_match_labels，全为False
            tk_exact_match_labels = np.array([False] * num_all)
            for dk in range(k):
                dk_pred_smiles_list = dk_pred_smiles_list_dict[dk]
                # 最底层计算的函数
                dk_exact_match_labels = judge_exact_match(dk_pred_smiles_list, golds_smiles_list)
                # tk_exact_match_labels 全为0 '|'操作实际上取的是dk_exact_match_labels的值
                tk_exact_match_labels = tk_exact_match_labels | dk_exact_match_labels
                # 记录每一层的exact_match数量
                metric_results['num_t%d_exact_match' % (dk + 1)] = tk_exact_match_labels.sum().item()
        
        elif metric == 'tanimoto':
            # 初始化tk_exact_match_labels，全为False
            tk_exact_match_labels = np.array([False] * num_all)
            for dk in range(k):
                dk_pred_smiles_list = dk_pred_smiles_list_dict[dk]
                # 最底层计算的函数
                dk_sims_labels = calculate_tanimoto_similarity(dk_pred_smiles_list, golds_smiles_list)
                # tk_exact_match_labels 全为0 '|'操作实际上取的是dk_exact_match_labels的值
                #tk_exact_match_labels = tk_exact_match_labels | dk_exact_match_labels
                # 记录每一层的exact_match数量
                metric_results['num_t%d_tanimoto' % (dk + 1)] = sum(dk_sims_labels)/len(dk_sims_labels)

        else:
            raise ValueError(metric)

    return metric_results,dk_sims_labels

# 定义一个函数，将jsonl文件转换为json对象
def convert_jsonl_to_json(jsonl_file_path):
    data_list = []
    # 打开jsonl文件，并逐行读取
    with open(jsonl_file_path, 'r', encoding='utf-8') as file:
        for line in file:
            # 将每一行转换为json对象
            json_object = json.loads(line.strip())
            # 将json对象添加到列表中
            data_list.append(json_object)
    # 返回列表
    return data_list

# 定义一个函数，用于提取SMILES
def re_formers_Lla(string):
    # 打印输入字符串
    # logging.info(string)
    # 定义正则表达式，用于匹配SMILES
    if string==None:
        return None
    regex = r"<SMILES>\s*([A-Za-z0-9]+)\s*(?:</SMILES>)?"
    # 查找匹配项
    match = re.search(regex, string)
    # 如果找到匹配项，返回第一个捕获组，否则返回None
    if match:
        return match.group(1)
    else:
        return None
# 定义一个函数，用于提取DFM的答案
def re_formers_DFM(string):
    # 定义正则表达式，用于匹配答案
    regex = r'"answer":\s*"([A-Za-z0-9]+)"'
    # 查找匹配项
    match = re.search(regex, string)
    # 如果找到匹配项，返回第一个捕获组，否则返回输入字符串
    if match:
        return match.group(1)
    else:
        return string
# 定义一个函数，用于列出指定文件夹中的jsonl文件
def list_json_files(folder_path):
    # 定义一个空列表，用于存储文件名
    file_list = []
    # 遍历文件夹中的所有文件
    for filename in os.listdir(folder_path):
        # 如果文件名以.jsonl结尾，将其添加到列表中
        if filename.endswith(".jsonl"):
            file_list.append(str(filename))
    # 返回列表
    return file_list

def add_json_files(file_path,sims):
    jsonl_data=[]
    with open(file_path, "r", encoding="utf-8") as infile:
        idx=0
        for line in infile:
            obj = json.loads(line)
            obj['tanimoto']=sims[idx]
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
            obj['tanimoto']=sims[idx]
            json_data.append(obj)
            idx=idx+1

    with open(file_path.replace('.jsonl','.json'), "w", encoding="utf-8") as outfile:
        for data in json_data:
            json.dump(data, outfile, ensure_ascii=False)
            outfile.write("\n")

def design_S(fold_path):
    # 定义文件路径
    fold_path = [os.path.join(fold_path,"molecule_design_S"),os.path.join(fold_path,"molecule_design")]
    for folder_path in fold_path:
        # 获取文件列表
        logging.basicConfig(
                level=logging.INFO,
                format='%(asctime)s - %(levelname)s - %(message)s',
                handlers=[
                    logging.FileHandler(folder_path+r'\output.txt','w',encoding='utf-8'),  # 写入文件
                    logging.StreamHandler()  # 同时输出到终端
                ],
                force=True
            )
        file_list = list_json_files(folder_path)
        # 遍历文件列表
        for name in file_list:
            #logging.info(name)
            # 定义文件路径
            path = folder_path + "\\" + name
            # 读取文件内容
            data = convert_jsonl_to_json(path)

            # 初始化列表
            golds = []
            pres = []
            # 遍历文件内容
            for i in range(len(data)):
                # 获取标准答案
                # gold = data[i].get("gold", "").replace(" ", "")
                if type(data[i]['gold'])==list:
                    gold = data[i].get("gold", "")[0].replace(" ", "").replace(";", ".")
                else:
                    gold = data[i].get("gold", "").replace(" ", "").replace(";", ".")

                # 获取预测答案
                if 'llasmol' in path:
                    pre = re_formers_Lla(data[i].get("answer", ""))
                else:
                    pre = data[i].get("answer", "")
                # 标记
                flag = 0
                # 存疑 判断分子类型
                ### modified
                #if "Branch" in gold or "Ring" in gold:
                if i<25 and 'SELFIES' in name:
                    try:
                        # 将gold以smiles格式解码
                        gold1 = sf.decoder(gold)
                        # 如果是smiles会转为"" ,不报错
                        if gold1 != "":
                            gold = gold1
                            flag = 1
                        else:
                            logging.info('gold selfies wrong'+f'{i+1}')
                    except:
                        gold = gold
                # if flag == 1:
                    try:
                        # 将pre以smiles格式解码
                        pre = sf.decoder(pre)
                        if pre == "":
                            pre=None
                    except:
                        #logging.info(i)
                        pre = None

                # 将预测答案和标准答案添加到列表中
                pres.append([pre])
                golds.append([gold])

            # 定义评估指标
            #metrics = ["exact_match"]
            metrics = ["tanimoto"]
            # 计算评估指标
            metric_results1,sims = calculate_smiles_metrics(pres, golds, metrics)

            add_json_files(path,sims)

            # 打印评估指标
            #logging.info(name, (metric_results1['num_t1_exact_match'])/metric_results1['num_all'])
            logging.info(f'{name}, tanimoto: {(metric_results1['num_t1_tanimoto'])}')
            #logging.info(f'{name}, no_answer: {(metric_results1['num_t1_no_answer'])}')
            logging.info(f'{name}, invalid: {(metric_results1['num_t1_invalid'])}')
            logging.info("-"*100)
