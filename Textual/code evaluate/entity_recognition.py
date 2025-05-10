import os
import json
import ast
import logging

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
    acc_list=[]
    with open(jsonl_file, 'r', encoding='utf-8') as file:
        for index,line in enumerate(file):
            # logging.info(index+1)
            data = json.loads(line)
            if data["answer"] is None:
                if type(data["gold"])!=list:
                    gold_labels = ast.literal_eval(data["gold"])
                else:
                    gold_labels=data["gold"]
                answer_labels = ['X'] * len(gold_labels)
            else:
                #if type(data["answer"])==list:
                #    data["answer"]=data["answer"][0]
                if type(data["gold"]) == str:
                    gold_labels = ast.literal_eval(data["gold"])
                elif type(data["gold"]) == list:
                    gold_labels = data["gold"]

                if type(data["answer"])==list:
                    try:
                        answer_labels = [s.replace(" ", "") for s in data["answer"]]
                    except:
                        answer_labels = ['X'] * len(gold_labels)
                elif type(data["answer"])==str:
                    data["answer"]=data["answer"].replace('<end>','').replace('```','')
                    try:
                        answer_labels = ast.literal_eval(data["answer"])
                    except:
                        if '["' not in data["answer"] or "['" not in data["answer"]:
                            try:
                                answer_labels = ast.literal_eval(data["answer"].replace('[','["').replace(']','"]').replace(', ','", "'))
                                #print(answer_labels)
                            except:
                                #print(222)
                                answer_labels = ['X'] * len(gold_labels)
                        else:
                            answer_labels = ['X'] * len(gold_labels)
                else:
                    answer_labels = ['X'] * len(gold_labels)
                
                # elif type(data["gold"]) == list:
                #     gold_labels =ast.literal_eval(data["gold"][0])

            correct_predictions = sum(1 for x, y in zip(answer_labels, gold_labels) if x == y)

            total_predictions = len(gold_labels)

            single_acc=correct_predictions/total_predictions
            acc_list.append(single_acc)

            total_correct += correct_predictions
            total_count += total_predictions
    # logging.info(total_correct)
    overall_accuracy = total_correct / total_count if total_count > 0 else 0
    return overall_accuracy,acc_list

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

# 如果该模块是主模块
def recognize(fold_path):
    # 指定文件夹路径
    fold_path = os.path.join(fold_path,"entity_recognition")
    logging.basicConfig(
                level=logging.INFO,
                format='%(asctime)s - %(levelname)s - %(message)s',
                handlers=[
                    logging.FileHandler(fold_path+r'\output.txt','w',encoding='utf-8'),  # 写入文件
                    logging.StreamHandler()  # 同时输出到终端
                ],
                force=True
            )
    # 调用函数，获取指定文件夹中的所有jsonl文件
    file_list = list_json_files(fold_path)
    # 遍历文件列表
    for name in file_list:
        # 指定文件路径
        path = fold_path + "\\" + name
        # 调用函数，计算准确率
        accuracy,acc_list = calculate_accuracy(path)
        add_json_files(path,acc_list)
        # 打印文件名和准确率
        logging.info(f'{name}, ACC: {accuracy}')
        logging.info('-'*100)