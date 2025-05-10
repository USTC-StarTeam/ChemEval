import pandas as pd
import json
import os
import logging

labels = [
    "Hepatobiliary disorders", "Metabolism and nutrition disorders",  "Eye disorders",
    "Musculoskeletal and connective tissue disorders", "Gastrointestinal disorders",
    "Immune system disorders", "Reproductive system and breast disorders",
    "Neoplasms benign, malignant and unspecified (incl cysts and polyps)", "Endocrine disorders",
    "Vascular disorders", "Blood and lymphatic system disorders", "Skin and subcutaneous tissue disorders",
    "Congenital, familial and genetic disorders", "Respiratory, thoracic and mediastinal disorders",
    "Psychiatric disorders", "Renal and urinary disorders", "Pregnancy, puerperium and perinatal conditions",
    "Ear and labyrinth disorders", "Cardiac disorders", "Nervous system disorders"
    ]
# 定义一个函数，将jsonl文件转换为json对象列表
def convert_jsonl_to_json(jsonl_file_path):
    # 创建一个空列表，用于存储json对象
    data_list = []
    # 以读取模式打开jsonl文件
    with open(jsonl_file_path, 'r', encoding='utf-8') as file:
        # 遍历文件每一行
        for line in file:
            # 将每一行转换为json对象
            json_object = json.loads(line.strip())
            # 将json对象添加到列表中
            data_list.append(json_object)
    # 返回json对象列表
    return data_list

# 定义一个函数，列出指定文件夹中的jsonl文件
def list_json_files(folder_path):
    # 创建一个空列表，用于存储jsonl文件名
    file_list = []
    # 遍历文件夹中的文件
    for filename in os.listdir(folder_path):
        # 如果文件名以.jsonl结尾，将其添加到列表中
        if filename.endswith(".jsonl"):
            file_list.append(str(filename))
    # 返回jsonl文件名列表
    return file_list

def add_json_files(file_path,acc_list):
    jsonl_data=[]
    with open(file_path, "r", encoding="utf-8") as infile:
        idx=0
        for line in infile:
            obj = json.loads(line)
            obj['ACC']=acc_list[idx]
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
            obj['ACC']=acc_list[idx]
            json_data.append(obj)
            idx=idx+1

    with open(file_path.replace('.jsonl','.json'), "w", encoding="utf-8") as outfile:
        for data in json_data:
            json.dump(data, outfile, ensure_ascii=False)
            outfile.write("\n")

def sider(fold_path):
    # 定义文件路径  
    fold_path = os.path.join(fold_path, "sider")
    logging.basicConfig(
                level=logging.INFO,
                format='%(asctime)s - %(levelname)s - %(message)s',
                handlers=[
                    logging.FileHandler(fold_path+r'\output.txt','w',encoding='utf-8'),  # 写入文件
                    logging.StreamHandler()  # 同时输出到终端
                ],
                force=True
            )
    # 获取文件列表
    file_list = list_json_files(fold_path)
    # 遍历文件列表
    for name in file_list:
        # 定义文件路径
        path = fold_path + "\\" + name
        # 读取文件内容
        data = convert_jsonl_to_json(path)
        # 初始化金标准列表
        golds = []
        # 初始化预测结果列表
        pres = []
        acc_list=[]
        # 遍历数据列表
        for i in range(len(data)):
            # 获取金标准结果
            gold = data[i].get("gold", "")
            # 如果预测结果为空或者为字符串，则将其添加到预测结果列表
            if type(data[i]) == str:
                pre = data[i]
            else:
                # 获取预测结果
                pre = data[i].get("answer", "")
            # 将预测结果和金标准结果添加到对应列表
            pres.append(pre)
            golds.append(gold)

        # 初始化预测结果列表
        pres1 = []
        # 遍历预测结果列表
        for pre in pres:
            # 如果预测结果为空或者为字符串，则将其添加到预测结果列表
            if not pre or type(pre) == str:
                # 初始化列表
                values = []
            elif list(pre.keys()) == labels:
                # 如果预测结果的键等于labels，则将其值添加到列表
                values = list(pre.values())
            else:
                # 如果预测结果的键不等于labels，则按照labels排序
                values = []
                for label in labels:
                    if label in pre.keys():
                        values.append(pre[label])
                    else:
                        values.append("None")

            # 将预测结果添加到列表
            pres1.append(values)

        # 初始化金标准列表
        golds1 = []
        # 遍历金标准结果列表
        for gold in golds:
            # 将金标准结果添加到列表
            values = list(gold.values())
            golds1.append(values)
        # 计算总的数量
        sum = len(golds1) * len(golds1[0])
        # 初始化正确数量
        correct = 0
        # 遍历金标准结果列表
        for i in range(len(golds1)):
            single_correct=0
            single_sum=len(golds1[i])
            # 遍历金标准结果列表里的每一行
            for j in range(len(golds1[0])):
                # 如果预测结果的长度大于j，则比较预测结果和金标准结果
                if len(pres1[i]) > j and pres1[i][j] == golds1[i][j]:
                    # 如果预测结果和金标准结果相等，则正确数量加1
                    correct += 1
                    single_correct+=1
            single_acc=single_correct/single_sum
            acc_list.append(single_acc)
        
        print(len(golds1))
        add_json_files(path,acc_list)
        # 输出结果
        logging.info(f'{name}, ACC Score: {correct/sum}')