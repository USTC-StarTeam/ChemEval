import numpy as np
import json
import re
import os
from collections import defaultdict
from sklearn.metrics.pairwise import cosine_similarity
import logging

def judge_string_exact_match(pred_string_list, golds_string_list):
    '''判断pred_string_list中的每个字符串是否与golds_string_list中的字符串完全匹配'''
    # 创建一个空列表，用于存储每个pred_string是否与golds_string_list中的字符串完全匹配
    exact_match_labels = []
    # 遍历pred_string_list中的每个字符串
    for pred_string, gold_string_list in zip(pred_string_list, golds_string_list):
        # 初始化exact_match为False
        exact_match = False
        # 遍历golds_string_list中的每个字符串
        for gold_string in gold_string_list:
            # 如果pred_string与gold_string相同，则exact_match为True，并跳出循环
            if pred_string == gold_string:
                exact_match = True
                break
        # 将exact_match添加到exact_match_labels列表中
        exact_match_labels.append(exact_match)
    # 返回exact_match_labels列表，表示每个pred_string是否与golds_string_list中的字符串完全匹配
    return np.array(exact_match_labels)

def calculate_formula_metrics(
        preds_formula_list,
        golds_formula_list,
        metrics=('exact_match')
):
    """
    计算分子式的度量。这里我们使用元素匹配（等于在论文中使用的精确匹配）作为默认值，比较原子数量并忽略顺序。
    例如，C5H8 == H8C5。
    """
    num_all = len(preds_formula_list)
    assert len(preds_formula_list) == len(golds_formula_list)
    try:
        k = len(preds_formula_list[0])
    except IndexError:
        print(preds_formula_list)
        raise
    dk_pred_formula_list_dict = dict()
    for dk in range(k):
        dk_pred_formula_list_dict[dk] = []
    for sample_formula_list in preds_formula_list:
        if sample_formula_list is None:
            for dk in range(k):
                dk_pred_formula_list_dict[dk].append('')
            continue
        assert len(sample_formula_list) == k
        for dk in range(k):
            item = sample_formula_list[dk]
            dk_pred_formula_list_dict[dk].append(item)
    golds_formula_list = [[small_item.strip() for small_item in item] for item in golds_formula_list]
    new_golds_formula_list = []
    for item in golds_formula_list:
        new_item = []
        for small_item in item:
            small_item = small_item.strip()
            assert small_item != ''
            new_item.append(small_item)
        new_golds_formula_list.append(new_item)
    golds_formula_list = new_golds_formula_list

    metric_results = {'num_all': num_all}

    tk_no_answer_labels = np.array([True] * num_all)
    for dk in range(k):
        dk_pred_formula_list = dk_pred_formula_list_dict[dk]
        dk_no_answer_labels = []
        for item in dk_pred_formula_list:
            if item == '' or item is None:
                dk_no_answer_labels.append(True)
            else:
                dk_no_answer_labels.append(False)
        dk_no_answer_labels = np.array(dk_no_answer_labels)
        tk_no_answer_labels = tk_no_answer_labels & dk_no_answer_labels
        metric_results['num_t%d_no_answer' % (dk + 1)] = tk_no_answer_labels.sum().item()

    for metric in metrics:
        if metric == 'exact_match':
            tk_exact_match_labels = np.array([False] * num_all)
            for dk in range(k):
                dk_pred_formula_list = dk_pred_formula_list_dict[dk]
                dk_exact_match_labels = judge_string_exact_match(dk_pred_formula_list, golds_formula_list)
                tk_exact_match_labels = tk_exact_match_labels | dk_exact_match_labels
                metric_results['num_t%d_exact_match' % (dk + 1)] = tk_exact_match_labels.sum().item()

        else:
            raise ValueError(metric)

    return metric_results

def parse_molecular_formula(formula):
        # 匹配元素符号（大写字母开头，后跟可选的小写字母）和可选的数字
        pattern = r'([A-Z][a-z]*)(\d*)'
        preds_atom_counts = defaultdict(int)

        for (element, count) in re.findall(pattern, formula):
            count = int(count) if count else 1
            preds_atom_counts[element] += count
    
        return dict(preds_atom_counts)

def get_union_atoms(vec1, vec2):
    # 获取两个向量中出现过的所有原子种类（合并后的维度）
    return sorted(set(vec1.keys()).union(set(vec2.keys())))

def dict_to_vector(atom_dict, atom_list):
    # 将原子频率字典转换为数值向量，按给定原子顺序排列
    return np.array([atom_dict.get(atom, 0) for atom in atom_list])

def atom_cosine_similarity(preds_vec_list, golds_vec_list):
    assert len(preds_vec_list) == len(golds_vec_list), "参考答案和预测答案的数量必须相等"
    sims=[]
    for i in range(len(preds_vec_list)):
        if preds_vec_list[i]:
            atom_list = get_union_atoms(preds_vec_list[i], golds_vec_list[i])
            v1 = dict_to_vector(preds_vec_list[i], atom_list).reshape(1, -1)
            v2 = dict_to_vector(golds_vec_list[i], atom_list).reshape(1, -1)
            sim=cosine_similarity(v1,v2)[0][0]
            sims.append(sim)
        else:
            sims.append(0)
    return sims

def l1_similarity(v1, v2):
    dist = np.sum(np.abs(v1 - v2))
    return 1 / (1 + dist)
    

def l2_similarity(v1, v2):
    dist = np.linalg.norm(v1 - v2)
    return 1 / (1 + dist)  # 同上
    

def atom_L1_similarity(preds_vec_list, golds_vec_list):
    assert len(preds_vec_list) == len(golds_vec_list), "参考答案和预测答案的数量必须相等"
    l1_similarities=[]
    for i in range(len(preds_vec_list)):
        if preds_vec_list[i]:
            atom_list = get_union_atoms(preds_vec_list[i], golds_vec_list[i])
            v1 = dict_to_vector(preds_vec_list[i], atom_list)
            v2 = dict_to_vector(golds_vec_list[i], atom_list)
            L1=l1_similarity(v1, v2)
            l1_similarities.append(L1)
        else:
            l1_similarities.append(0)
    return l1_similarities

def atom_L2_similarity(preds_vec_list, golds_vec_list):
    assert len(preds_vec_list) == len(golds_vec_list), "参考答案和预测答案的数量必须相等"
    l2_similarities=[]
    for i in range(len(preds_vec_list)):
        if preds_vec_list[i]:
            atom_list = get_union_atoms(preds_vec_list[i], golds_vec_list[i])
            v1 = dict_to_vector(preds_vec_list[i], atom_list)
            v2 = dict_to_vector(golds_vec_list[i], atom_list)
            L2=l2_similarity(v1, v2)
            l2_similarities.append(L2)
        else:
            l2_similarities.append(0)
    return l2_similarities
    

# 定义一个函数，将jsonl文件转换为json对象列表
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

# 定义一个函数，用于提取LlaSMol MOLFORMULA标签的内容
def re_formers_Lla(string):
    # print(string)
    # 定义正则表达式，用于匹配MOLFORMULA标签及其内容
    regex = r"<MOLFORMULA>\s*([A-Za-z0-9]+)\s*(?:</MOLFORMULA>)?"
    # 查找匹配项
    match = re.search(regex, string)
    # 如果找到匹配项，返回匹配的内容，否则返回None
    if match:
        return match.group(1)
    else:
        return None
# 定义一个函数，用于提取CHemDFM answer标签的内容
def re_formers_DFM(string):
    # 定义正则表达式，用于匹配answer标签及其内容
    regex = r'"answer":\s*"([A-Za-z0-9]+)"'
    # 查找匹配项
    match = re.search(regex, string)
    # 如果找到匹配项，返回匹配的内容，否则返回原始字符串
    if match:
        return match.group(1)
    else:
        return string

def list_json_files(folder_path):
    '''获取文件夹中所有jsonl文件名'''
    file_list = []
    # 遍历文件夹中的所有文件
    for filename in os.listdir(folder_path):
        # 如果文件名以.jsonl结尾
        if filename.endswith(".jsonl"):
            # 将文件名添加到文件列表中
            file_list.append(str(filename))
    # 返回文件列表
    return file_list

def add_json_files(file_path,sims,l1_similarities,l2_similarities):
    jsonl_data=[]
    with open(file_path, "r", encoding="utf-8") as infile:
        idx=0
        for line in infile:
            obj = json.loads(line)
            obj['consine']=sims[idx]
            obj['L1 similarity']=l1_similarities[idx]
            obj['L2 similarity']=l2_similarities[idx]
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
            obj['consine']=sims[idx]
            obj['L1 similarity']=l1_similarities[idx]
            obj['L2 similarity']=l2_similarities[idx]
            json_data.append(obj)
            idx=idx+1

    with open(file_path.replace('.jsonl','.json'), "w", encoding="utf-8") as outfile:
        for data in json_data:
            json.dump(data, outfile, ensure_ascii=False)
            outfile.write("\n")

def design_F(root_path):
    # 定义文件路径
    fold_path = os.path.join(root_path, "molecule_design_F")

    logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(fold_path+r'\output.txt','w',encoding='utf-8'),  # 写入文件
                logging.StreamHandler()  # 同时输出到终端
            ],
            force=True
        )

    # 获取文件夹中的所有json文件
    file_list = list_json_files(fold_path)
    # 遍历文件列表
    for name in file_list:
        # 定义文件路径
        path = fold_path + "\\" + name
        
        # 将json文件转换为字典
        data = convert_jsonl_to_json(path)

        # 初始化列表
        golds = []
        pres = []
        # 遍历字典中的数据
        for i in range(len(data)):
            # 获取gold值
            gold = data[i].get("gold", "")
            # 获取pre值
            pre = data[i].get("answer", "")
            if pre and 'llasmol' in path:
                pre = re_formers_Lla(data[i].get("answer", ""))
            
            # 对pre值进行处理
            if pre:
                pre_f = re_formers_DFM(pre)
            else:
                pre_f=pre
            # 将处理后的pre值添加到pres列表中
            pres.append([pre_f])
            # 将gold值添加到golds列表中
            golds.append([gold])

        # 定义评估指标
        metrics = ['cosine','L1','L2']
        # 计算评估指标
        metric_results={}

        preds_vec_list=[]
        golds_vec_list=[]

        for i in range(len(pres)):
            if pres[i][0]:
                preds_vec = parse_molecular_formula(pres[i][0])
                preds_vec_list.append(preds_vec)
            else:
                preds_vec_list.append(None)
            golds_vec = parse_molecular_formula(golds[i][0])
            golds_vec_list.append(golds_vec)

        for metric in metrics:
            if 'cosine'==metric:
                sims=atom_cosine_similarity(preds_vec_list, golds_vec_list)
                metric_results[metric]=sum(sims)/len(sims)
                # 打印结果
                logging.info(f'{name}, {metric}: {metric_results[metric]}')
            elif 'L1'==metric:
                l1_similarities=atom_L1_similarity(preds_vec_list, golds_vec_list)
               
                try:
                    metric_results[metric]=sum(l1_similarities)/len(l1_similarities)
                except:
                    L1_sum=0
                    L1_count=0
                    for L1 in l1_similarities:
                        if L1!=None:
                            L1_sum+=L1
                            L1_count+=1
                    try:
                        metric_results[metric]=L1_sum/L1_count
                    except:
                        metric_results[metric]='infin'
                # 打印结果
                logging.info(f'{name}, {metric}: {metric_results[metric]}')
            elif 'L2'==metric:
                l2_similarities=atom_L2_similarity(preds_vec_list, golds_vec_list)
                
                try:
                    metric_results[metric]=sum(l2_similarities)/len(l2_similarities)
                except:
                    L2_sum=0
                    L2_count=0
                    for L2 in l2_similarities:
                        if L2!=None:
                            L2_sum+=L2
                            L2_count+=1
                    try:
                        metric_results[metric]=L2_sum/L2_count
                    except:
                        metric_results[metric]='infin'
                # 打印结果
                logging.info(f'{name}, {metric}: {metric_results[metric]}')
        if len(l2_similarities)!=0 and len(l2_similarities)!=0:
            add_json_files(path,sims,l2_similarities,l2_similarities)
        logging.info("-"*100)

        #metric_results = calculate_formula_metrics(pres, golds, metrics)
        # 打印评估指标
        #print(name, (metric_results['num_t1_exact_match'])/metric_results['num_all'])