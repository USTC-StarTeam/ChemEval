import json
import re


def read_file(filename,pattern=1):
    if pattern == 1:
        data_list = []
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                if line.strip():
                    data_list.append(json.loads(line))
        return data_list
    elif pattern == 2:
        data_list = []
        with open(filename, 'r', encoding='utf-8') as f:
            data_list = json.load(f)
        return data_list

def write_file(filename, datas,pattern=1):
    if pattern == 1:
        print(filename+" : ",len(datas))
        with open(filename, 'w', encoding='utf-8') as f:
            for data in datas:
                f.write(json.dumps(data, ensure_ascii=False))
                f.write('\n')
        f.close()
    elif pattern == 2:
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(datas, f, ensure_ascii=False,indent=4)

def extract_gold(text):
    match = re.search(r'其对应的SMILES为：\s*(.+?)\n', text,re.DOTALL)
    if match:
        smiles = match.group(1)
        return smiles
    else:
        print("SMILES not found.")
        print(text)
        raise Exception

def get_goldlist():
    gold_list = []
    root_path = r'D:\project\nips-chemeval\data\t1-1\化学多模态数据集\3.分子图文理解\3.分子名称翻译\4.SMILES转IUPAC\SMILES转IUPAC_23.jsonl'
    data_list = read_file(root_path)
    for data in data_list:
        gold_list.append(extract_gold(data['query']))
    return gold_list

gold_list = get_goldlist()

import os
import glob

# 固定路径的上级目录
base_path = r"D:\project\nips-chemeval\pipline\data\机评结果\按t拆分"

# 遍历所有子文件夹（即 model_name）
model_dirs = [os.path.join(base_path, d) for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]

# 构建目标文件路径的通配符
target_rel_path = r"output2_me\molecule_design_I\3.分子名称翻译.4.SMILES转IUPAC.SMILES转IUPAC_23.jsonll.jsonl"

# 收集所有匹配的文件路径
matched_files = []
for model_dir in model_dirs:
    full_path = os.path.join(model_dir, target_rel_path)
    if os.path.exists(full_path):
        matched_files.append(full_path)
matched_files = [r"D:\project\nips-chemeval\2chemeval_code-2\output_gpt-1t\output2_me\molecule_design_I\3.分子名称翻译.4.SMILES转IUPAC.SMILES转IUPAC_23.jsonll.jsonl",
                 r"D:\project\nips-chemeval\2chemeval_code-2\output_gpt-2t\output2_me\molecule_design_I\3.分子名称翻译.4.SMILES转IUPAC.SMILES转IUPAC_23.jsonll.jsonl",
                 r"D:\project\nips-chemeval\2chemeval_code-2\output_gpt-3t\output2_me\molecule_design_I\3.分子名称翻译.4.SMILES转IUPAC.SMILES转IUPAC_23.jsonll.jsonl"]
for path in matched_files:
    ori_data = read_file(path)
    assert len(ori_data) == len(gold_list)
    for data,gold_smi in zip(ori_data,gold_list):
        data['gold_smi'] = gold_smi
    write_file(path,ori_data)
