import json
import os
import logging
import re

def calculate_f1_score(true_entities, predicted_entities):
    '''
    计算F1分数
    :param true_entities: 真实实体集合
    :param predicted_entities: 预测实体集合
    :return: F1分数
    '''
    true_entities = true_entities.replace(" ", "")
    true_entities = set(true_entities.split(','))
    predicted_entities = predicted_entities.replace(" ", "")
    predicted_entities = set(predicted_entities.split(','))
    true_positive = len(true_entities & predicted_entities)
    precision = true_positive / len(predicted_entities) if len(predicted_entities) > 0 else 0
    recall = true_positive / len(true_entities) if len(true_entities) > 0 else 0

    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    return f1_score


# 定义一个函数，用于计算准确率
def calculate_accuracy(golds, pres):
    # 初始化正确计数和总计数
    correct_count = 0
    total_count = len(golds)
    f1_scores=[]
    # 遍历golds中的每个元素
    for i in range(len(golds)):
        # 将true_output转换为小写
        true_output = golds[i].lower()
        # my_output = pres[i][:-4].lower()
        # 如果pres[i]是一个字符串
        if isinstance(pres[i], list):
            # 遍历pres[i]中的每个元素，将它们转换为小写
            my_outputs = [item[:].lower() for item in pres[i] if isinstance(item, str)]  # 处理每个字符串元素
        else:
            # 将my_output转换为小写
            my_output = pres[i][:].lower()

        # 计算f1_score
        f1_score = calculate_f1_score(true_output, my_output)
        f1_scores.append(f1_score)
        # 将f1_score加到正确计数中
        correct_count += f1_score

    # 返回正确计数占总计数的比例
    return correct_count / total_count,f1_scores

# 定义一个函数，用于列出文件夹中的json文件
def list_json_files(folder_path):
    # 初始化文件列表
    file_list = []
    # 遍历文件夹中的每个文件
    for filename in os.listdir(folder_path):
        # 如果文件以.jsonl结尾
        if filename.endswith(".jsonl"):
            # 将文件名添加到文件列表中
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

def entity_extract(fold_path):

    # 定义一个字典，用于存储实体及其对应缩写
    chem_entity_dict = {
        "H ( 2 ) O ( 2 ), SiO ( 2 )": "H2O2, SiO2",
        "4 - hydroxy - 1 - ( 3 - pyridyl ) - 1 - butanone (HPB)": "HPB",
        "side chains of compound ( - ) - 1a (SCH 900229)": "SCH 900229",
        "nitric oxide (( . ) NO)": "( . ) NO",
        "4 - hydroxy - 1 - ( 3 - pyridyl ) - 1 - butanone":	"HPB",
        "side chains of compound ( - ) - 1a": "SCH 900229",
        "nitric oxide": "( . ) NO",
        # add
        "HAND - SANT - SLIDE": "HSS",
        "glucose - dependent insulinotropic polypeptide": "GIP",
        "monoacylglycerol lipase": "MAGL",
        "short - chain fatty acid": "SCFA",
        "total polyphenols": "TP",
        "structure activity relationships": "SAR",
        " life and octanol - water distribution coefficient": "log D",
        "glycyrrhetinic acid": "GA",
        "Lysyl oxidase": "LO",
        "Manganese": "Mn"

    }
    relation_dict = {
        "tranexamic acid (tAMCA)": "tAMCA",
        "thrombotic microangiopathy (TMA)": "TMA",
        "trimethaphan (TMP)": "TMP",
        "prostaglandin E1 (PGE1)": "PGE1",
        "calcium chloride (CaCl(2))": "CaCl(2)",
        "thoracic aortic aneurysm (TAA)": "TAA",
        "tranexamic acid": "tAMCA",
        "thrombotic microangiopathy": "TMA",
        "trimethaphan": "TMP",
        "prostaglandin E1": "PGE1",
        "calcium chloride": "CaCl(2)",
        "thoracic aortic aneurysm": "TAA",
        # add
        "counter cyanide": "CN",
        "acetylcholine": "ACh",
        "organophosphorus": "OP",
        "diisopropylfluorophosphate": "DFP",
        "acetylcholinesterases": "AChEs",
        "butyrylcholinesterases": "BChEs",
        "N-methyl-D-aspartate": "NMDA",
        "venous thromboembolism": "VTE",
        "antiepileptic drug": "AED",
        "Carbamazepine": "CBZ",
        "gamma-Vinyl GABA": "GVG",
        "Parkinson's disease": "PD",
        "metalloproteinase": "ADAM",
        "metalloproteinases": "MMPs"

    }
    additive_dict = {
        "pivalic acid (PivOH)": "PivOH",
        "2,2,6,6-tetramethylpiperidinooxy (TEMPO)": "TEMPO",
        "pivalic acid": "PivOH",
        "2,2,6,6-tetramethylpiperidinooxy": "TEMPO"
    }
    solvent_dict = {
        "N,N-dimethylformamide (DMF)": "DMF",
        "1,1,1,3,3,3-hexafluoroisopropanol (HFIP)": "HFIP",
        "N,N-dimethylformamide": "DMF",
        "1,1,1,3,3,3-hexafluoroisopropanol": "HFIP"
    }

    
    # fold_path = r"C:\Users\xiangli67\Downloads\3880数据集\o3-mini\实体抽取"
    fold_path = os.path.join(fold_path,"entrty_extraction")
    logging.basicConfig(
                level=logging.INFO,
                format='%(asctime)s - %(levelname)s - %(message)s',
                handlers=[
                    logging.FileHandler(fold_path+r'\output.txt','w',encoding='utf-8'),  # 写入文件
                    logging.StreamHandler()  # 同时输出到终端
                ],
                force=True
            )
    file_list = list_json_files(fold_path)
    for name in file_list:
        path = fold_path + "\\" + name
        # path = "split_jsonl_files\\chemical_entity_recognition.json"

        data = []
        with open(path, "r", encoding='utf-8') as file:
            for line in file:
                json_object = json.loads(line.strip())
                data.append(json_object)

        golds = []
        pres = []
        for i in range(len(data)):
            gold = data[i].get("gold", "")
            if type(data[i]) == str:
                pre = str(data[i])
            else:
                pre = str(data[i].get("answer", ""))
            if "化学命名实体识别" in name or "Chemical_Named_Entity_Recognition" in name or "Chemical_Named_Entity_Recognition" in name:
                for key, value in chem_entity_dict.items():
                    gold = gold.replace(key, value)
                    pre = pre.replace(key, value)
            if "化学实体关系分类" in name or "Chemical_Entity_Relationship_Classification" in name:
                for key, value in relation_dict.items():
                    gold = gold.replace(key, value)
                    pre = pre.replace(key, value)
            if "合成反应添加剂抽取" in name or "thetic_Reaction_Additive_Extraction" in name:
                for key, value in additive_dict.items():
                    gold = gold.replace(key, value)
                    pre = pre.replace(key, value)
            if "合成反应溶剂抽取" in name or "Synthetic_Reaction_Solvent_Extraction" in name:
                for key, value in solvent_dict.items():
                    gold = gold.replace(key, value)
                    pre = pre.replace(key, value)
            if "催化类型抽取" in name or "Catalysis_Type_Extraction" in name  or "化学反应类型识别归纳" in name or "Reaction_Type_Recognition_and_Induction" in name:
                pre = pre.replace("reactions", "").replace("reaction", "").replace(",", ".")
                gold = gold.replace("reactions", "").replace("reaction", "").replace(",", ".")
            if "合成反应温度抽取" in name or "Reaction_Temperature_Extraction" in name:
                pre = pre.replace("℃", "").replace("°C", "").replace("degrees Celsius", "").replace("degrees Celsius.", "")
                gold = gold.replace("℃", "").replace("°C", "")
            if "合成反应时间抽取" in name:
                pre = pre.replace("hours", "").replace("hour", "").replace("h", "").replace('.0','')  
                if ' min' in pre:
                    try:
                        pre=str(float(pre.replace('min','').strip())/60).replace('.0','')  
                    except:
                        match=re.search(r'(\d+) min',pre)
                        if match:
                            pre=match.group(1).replace('.0','')
                gold = gold.replace("h", "")
            if "产率性能抽取" in name:
                pre=pre.replace('answer:','').replace('percent','%').replace('% to ','-').replace('%-','-')
                # match_p = re.search(r'[-+]?\d*\.?\d+(?=\D*$)', pre)
                # if match_p:
                #     pre = match_p.group()
                # else:
                #     pre = None
                # match_g = re.search(r'[-+]?\d*\.?\d+(?=\D*$)', gold)
                # if match_g:
                #     gold = match_g.group()
                # else:
                #     gold = None

            pres.append(pre)
            golds.append(gold)
        accuracy,f1_scores = calculate_accuracy(golds, pres)
        add_json_files(path,f1_scores)
        # if "化学命名实体识别" in name or "化学实体关系分类" in name or "合成反应添加剂抽取" in name or "合成反应溶剂抽取" in name:
        logging.info(f"{name}, F1 score avg: {accuracy:.4f}")