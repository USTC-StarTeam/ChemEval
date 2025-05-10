import os, pandas as pd
import nltk
from nltk.translate.meteor_score import meteor_score
import sacrebleu
import json
import numpy as np
from rdkit import Chem
import selfies as sf
import math
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


base_dir = r"D:\project\nips-chemeval\pipline\data\机评结果\按t拆分\phi-vision-134-3t\output2_me\Before_rating\\"
files = list_json_files(base_dir)
files = [file for file in files if file.endswith(".jsonl")]
for filen in files:
    # print(filen)
    if "填空" in filen:
        res = 0
        data_list = read_file(base_dir+filen)
        whole = len(data_list)
        for data in data_list:
            res+=(data['answer'])
        print(filen,"acc:",res/whole)
    elif "识别" in filen:
        data_list = read_file(base_dir + filen)
        whole = len(data_list)
        right = 0
        for data in data_list:
            if data['answer'] == "正确":
                right+=1
        print(filen,"ACC:",right/whole)
    else:
        res = 0
        data_list = read_file(base_dir+filen)
        whole = len(data_list)*5
        for data in data_list:
            if type(data['answer']) != int:
                if data['answer'] is None:
                    data['answer'] = 0
                else:
                    data['answer'] = int(data['answer'][0])
            try:
                res+=data['answer']
            except:
                print(data['answer'])
                raise Exception
        print(filen,"score:",res/whole*100)