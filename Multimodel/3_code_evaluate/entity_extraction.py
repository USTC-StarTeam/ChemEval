import json
import os

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
        # 将f1_score加到正确计数中
        correct_count += f1_score

    # 返回正确计数占总计数的比例
    return correct_count / total_count

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

def entity_extract(fold_path,excel_path):

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
    file_list = list_json_files(fold_path)
    result_data = []
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
                pre = pre.replace("℃", "").replace("°C", "")
                gold = gold.replace("℃", "").replace("°C", "")
            pres.append(pre)
            golds.append(gold)
        accuracy = calculate_accuracy(golds, pres)
        result_data.append([name,round(accuracy,4)])
        # if "化学命名实体识别" in name or "化学实体关系分类" in name or "合成反应添加剂抽取" in name or "合成反应溶剂抽取" in name:
        print(name, f"ACC: {accuracy:.4f}")
    import pandas as pd
    # 创建DataFrame
    df = pd.DataFrame(result_data, columns=["File Name", "ACC Score"])
    sheet_name = "Entity_extraction"
    # 读取现有Excel文件（如果存在）
    if os.path.exists(excel_path):
        with pd.ExcelWriter(excel_path, mode="a", if_sheet_exists="replace", engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        with pd.ExcelWriter(excel_path, mode="w", engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
