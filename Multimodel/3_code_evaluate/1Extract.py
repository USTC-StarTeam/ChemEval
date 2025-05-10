import json
import ast
import re
import os
import shutil

def is_valid_dict_format(s):
    '''
    判断字符串是否为有效的字典格式

    参数：
    s : str
        需要判断的字符串

    返回值：
    bool
        如果是有效的字典格式返回True，否则返回False
    '''
    try:
        # 使用ast.literal_eval()尝试将字符串转换为字典
        ast.literal_eval(s)
        # 如果转换成功，返回True
        return True
    except (SyntaxError, ValueError):
        # 如果转换失败，返回False
        return False

def find_last_ABCD_letter(text):
    # 使用正则表达式在文本中查找所有匹配的字母[ABCD]
    matches = 0
    if type(text) == list:
        text = " ".join(text)
    if text != None:
        matches = re.findall(r'[ABCD]', text)
    # 如果找到了匹配项
    if matches:
        # 返回最后一个匹配项
        return matches[-1]
    else:
        # 如果没有找到匹配项，返回None
        return None

def extract_answer(file_path):
    # 不同的类型执行不同的函数
    if "评分前" in file_path  or "Before_rating" in file_path or "判断" in file_path or "True_False_Task" in file_path:
        reagent_extract_answer(file_path)
    elif "反应时间推荐" in file_path or "Reaction_Time_Recommendation" in file_path or "反应温度推荐" in file_path or "Reaction_Temperature_Recommendation" in file_path:
        avg_extract_answer(file_path)
    elif "SIDER" in file_path:
        sider_extract_answer(file_path)
    elif ("分类" in file_path and "关系分类" not in file_path and "关系抽取" not in file_path and "化学文献主题分类" not in file_path and "判断" not in file_path) or ("classification" in file_path and "化学文献主题分类" not in file_path and "relation_extraction" not in file_path and "Reaction_Type_Recognition_and_Induction" not in file_path and "True_False_Task" not in file_path and "Chemical_Literature_Topic_Classification" not in file_path):
        class_extract_answer(file_path)
    elif "regression" in file_path:
        reg_extract_answer(file_path)
    elif "试剂" in file_path  or "reagent" in file_path  or "Reagent" in file_path or "抽取" in file_path or "extraction" in file_path or "Extraction" in file_path :
        reagent_extract_answer(file_path)
    elif "基于分子结构描述分子的物理化学性质" in file_path or "简答任务" in file_path or "化学反应类型识别归纳" in file_path or "生成" in file_path or "计算任务" in file_path or "化学论文摘要生成" in file_path or "研究内容提纲生成" in file_path or "化学文献主题分类" in file_path or "Molecular_Name_Generation_from_Text_Description" in file_path or "Short_Answer_Task" in file_path or "Reaction_Type_Recognition_and_Induction" in file_path or "Generation" in file_path or "Calculation_Task" in file_path or "Chemical_Paper_Abstract_Generation" in file_path or "Research_Outline_Generation" in file_path or "Chemical_Literature_Topic_Classification" in file_path:
        gen_extract_answer(file_path)
    elif ("分子" in file_path) or ("Molecular" in file_path) or ("molecule" in file_path):
        mol_extract_answer(file_path)

# SIDER类型
def sider_extract_answer(file_path):
    # 读取文件
    all_data, gold_data = read_file(file_path)
    output = []
    print(len(all_data))
    for i in range(len(all_data)):
        if type(all_data[i]) == str:# and is_valid_dict_format(all_data[i]):
            data = (all_data[i].replace("\\n", "").replace("\\", "").replace("{{","{")
                    .replace("// Placeholder as we do not have data to predict", "")
                    .replace("// This is the specific question asked", "")
                    .replace("// Based on the specific question asked ", "")
                    .replace("// Hypothetical answer to the specific question asked", "")
                    .replace("// Speculative answer based on the original question ", "")
                    .replace("// As per the specific question asked about renal and urinary disorders.", "")
                    .replace("// Hypothetical prediction", "")
                    .replace("// Hypothetical answer based on the question", "")
                    .replace("based on the question", "")
                    .replace("```", "")
                    .replace("<end>", "")
                    .replace('}"}', '}}')
                    .replace('"{', '{'))

            # 查找json开始和结束索引
            json_start_index = data.find('{')
            json_end_index = data.find('}')
            res = data[json_start_index: json_end_index+1]
            # print(json_start_index,json_end_index)
            if json_start_index != -1 and json_end_index != -1:
                if isinstance(data, str):
                    try:
                        data = ast.literal_eval(data[json_start_index: json_end_index+1])
                    except SyntaxError as e:
                        try:
                            data = ast.literal_eval(data[json_start_index: json_end_index+1]+"}")
                        except SyntaxError as e:
                            data = data[json_start_index: json_end_index+1]

                if type(data) == dict and "answer" in data.keys():
                    question_answer_dict = {"answer": data["answer"], "gold": gold_data[i]}
                else:
                    question_answer_dict = {"answer": data, "gold": gold_data[i]}
                output.append(question_answer_dict)
        else:
            question_answer_dict = {"answer": all_data[i], "gold": gold_data[i]}
            output.append(question_answer_dict)
            print(output)
    # 写入文件
    write_file(file_path, output)

# 分类问题抽取
def class_extract_answer(file_path):
    # 读取文件
    all_data, gold_data = read_file(file_path)
    output = []
    # 遍历文件中的每一行
    for i in range(len(all_data)):
        # 如果文件路径中含有选择，则找到最后一个ABCD字母
        if "选择" in file_path or "Multiple_Choice" in file_path or "MMLU" in file_path or "SciQ" in file_path:
            all_data[i] = find_last_ABCD_letter(all_data[i])
        else:
            if isinstance(all_data[i], list):
                all_data[i] = all_data[i][0]
            else:
                all_data[i] = all_data[i]
            # 将数据转换为小写
            if all_data[i] is not None:
                data = all_data[i].lower()
            else:
                data = "null"

            if "query" in data:
                answer_end_index = all_data[i].find("query")
                data = data[:answer_end_index]
            # 如果数据中含有yes和no，则将数据置为空
            if "yes" in data and "no" in data:
                all_data[i] = None
                # change 把yes换成了数字
            elif "yes" in data or 'yes' in data:
                all_data[i] = "Yes"
            elif "no" in data or 'no' in data:
                all_data[i] = "No"
            else:
                all_data[i] = None

        # 创建一个字典，包含答案和金标准
        question_answer_dict = {"answer": all_data[i], "gold": gold_data[i]}
        # 将字典添加到输出列表中
        output.append(question_answer_dict)

    # 将输出列表写入文件
    write_file(file_path, output)


# 均值抽取
def avg_extract_answer(file_path):
    # 从文件路径中读取数据
    all_data, gold_data = read_file(file_path)
    output = []
    # 遍历数据
    for i in range(len(all_data)):
        # 处理数据，将-替换为空格
        if type(all_data[i]) != str:
            all_data[i] = None
        else:
            processed_data = all_data[i].replace("-", " ")
            # 使用正则表达式查找数字
            numbers = re.findall(r'[-+]?\d*\.?\d+', processed_data)
            # 如果找到数字
            if numbers:
                # 将数字转换为浮点数
                numbers = [float(num) for num in numbers]
                # 对数字求和，然后除以数字的数量，得到平均值
                all_data[i] = sum(numbers) / len(numbers)
            # 如果没有找到数字
            else:
                # 将该数据设置为None
                all_data[i] = None

            # 创建一个字典，包含答案和金标准
            question_answer_dict = {"answer": all_data[i], "gold": gold_data[i]}
            # 将字典添加到输出列表中
            output.append(question_answer_dict)

    # 将输出列表写入文件路径
    write_file(file_path, output)

# 回归问题抽取
def reg_extract_answer(file_path):
    # 读取文件
    all_data, gold_data = read_file(file_path)
    # 创建一个空列表，用于存储处理后的数据
    output = []
    # 遍历all_data中的每一行
    for i in range(len(all_data)):
        # print(i)
        # 使用正则表达式查找数字
        match = None
        if type(all_data[i]) ==str:
            match = re.search(r'[-+]?\d*\.?\d+(?=\D*$)', all_data[i])
        elif all_data[i] == None:
            match = 0
        # 如果找到数字
        if match:
            # 将找到的数字转换为浮点数
            all_data[i] = float(match.group())
        else:
            # 如果没找到数字，将该行数据设置为None
            all_data[i] = None

        if type(gold_data[i]) ==str:
            match = re.search(r'[-+]?\d*\.?\d+(?=\D*$)', gold_data[i])
        elif gold_data[i] == None:
            match = 0
        # 如果找到数字
        if type(gold_data[i]) ==str and match:
            # 将找到的数字转换为浮点数
            gold_data[i] = float(match.group())
        elif type(gold_data[i]) ==str:
            # 如果没找到数字，将该行数据设置为None
            gold_data[i] = None

        # 创建一个字典，用于存储该行的答案和金标答案
        question_answer_dict = {"answer": all_data[i], "gold": gold_data[i]}
        # 将该字典添加到输出列表中
        output.append(question_answer_dict)

    # 将处理后的数据写入文件
    write_file(file_path, output)

# other抽取
def reagent_extract_answer(file_path):
    if "ifly_Chem_IE" in file_path:
        return
    all_data, gold_data = read_file(file_path)
    output = []
    for i in range(len(all_data)):
        if isinstance(all_data[i], dict):
            all_data[i] = "".join(str(value) for value in all_data[i].values())
        elif isinstance(all_data[i], list):
            all_data[i] = ".".join(all_data[i])
        else:
            all_data[i] = str(all_data[i])
        # if "M_forward_reaction" in file_path:
        #     print(f'all_data_00[i] {all_data[i]}')
        if all_data[i] == None:
            question_answer_dict = {"answer": all_data[i], "gold": gold_data[i]}
            output.append(question_answer_dict)
            continue
        all_data[i] = (all_data[i].replace("\n", "").replace("\\", "").replace("'", '\"')
                       .replace("assistant", "answer").replace(": ", ":").replace('"answer":"{"answer"', '"answer"')
                       )
        if "query" in all_data[i]:
            answer_end_index = all_data[i].find("query")
            all_data[i] = all_data[i][:answer_end_index-5]
        # if "M_forward_reaction" in file_path:
        #     print(f'all_data_0[i] {all_data[i]}')
        if '"answer":"' in all_data[i]:
            answer_start_index = all_data[i].find('"answer":"') + len('"answer":"')
            answer_end_index = all_data[i].rfind('"', answer_start_index)
            if answer_start_index == -1 + len('"answer":"'):
                all_data[i] = None
            elif answer_end_index == -1:
                all_data[i] = all_data[i][answer_start_index:]
            else:
                all_data[i] = all_data[i][answer_start_index:answer_end_index]
            # if "M_forward_reaction" in file_path:
            #     print(f'all_data_1[i] {all_data[i]}')
        elif '"answer":[' in all_data[i]:
            answer_start_index = all_data[i].find('"answer":[') + len('"answer":[')
            answer_end_index = all_data[i].find(']"}', answer_start_index)
            if answer_start_index == -1 + len('"answer":['):
                all_data[i] = None
            elif answer_end_index == -1:
                all_data[i] = all_data[i][answer_start_index - 1:]
            else:
                all_data[i] = all_data[i][answer_start_index-1: answer_end_index+1]
        # 处理<SMILES> </SMILES>标签对
        if "<SMILES>" in all_data[i] and "</SMILES>" in all_data[i]:
            answer_start_index = all_data[i].find('<SMILES>') + len('<SMILES>') + 2
            answer_end_index = all_data[i].find('</SMILES>') - 2
            if answer_start_index == -1 + len('<SMILES>'):
                all_data[i] = None
            elif answer_end_index == -1:
                all_data[i] = all_data[i][answer_start_index - 1:]
            else:
                all_data[i] = all_data[i][answer_start_index - 1: answer_end_index + 1]
        if "is:" in all_data[i] and "</SMILES>" in all_data[i]:
            answer_start_index = all_data[i].find('is:') + len('is:') + 1
            answer_end_index = all_data[i].find('</SMILES>') - 2
            if answer_start_index == -1 + len('is:'):
                all_data[i] = None
            elif answer_end_index == -1:
                all_data[i] = all_data[i][answer_start_index - 1:]
            else:
                all_data[i] = all_data[i][answer_start_index - 1: answer_end_index + 1]

        try:
            all_data[i] = eval(all_data[i])
        except:
            all_data[i] = all_data[i].replace('"', '').replace("{", "").replace("}", "")
        if type(all_data[i]) == str and len(all_data[i]) != 0 and all_data[i][0] == '[' and all_data[i][-1] == ']' \
                and all_data[i].count('[') == 1 and all_data[i].count(']') == 1:
            all_data[i] = all_data[i][1:-1].replace(",", ".").split(".")
        if isinstance(all_data[i], dict):
            all_data[i] = "".join(str(value) for value in all_data[i].values())
        question_answer_dict = {"answer": all_data[i], "gold": gold_data[i]}
        output.append(question_answer_dict)

    write_file(file_path, output)

# 分子式抽取
def mol_extract_answer(file_path):
    # 读取文件
    all_data, gold_data = read_file(file_path)
    output = []
    # 遍历文件中的每一行
    for i in range(len(all_data)):
        if type(all_data[i]) == list:
            all_data[i] = "".join(all_data[i])
        if all_data[i] == None or all_data[i] =="None":
            all_data[i] = ""
        # 去除空格、换行符、单引号、双引号、反斜杠
        all_data[i] = (all_data[i].replace(" ", "").replace("\n", "").replace("'", '"').replace("{\"answer\":\"{\"answer\":\"", "{\"answer\":\"")
                       .replace("\\", "").replace(":", "").replace('""','"'))
        # 如果包含"answer"，则获取"answer"的起始和结束索引
        if '"answer"' in all_data[i]:
            answer_start_index = all_data[i].find('"answer"') + len('"answer"')
            answer_end_index = all_data[i].find('"', answer_start_index)
            # 如果起始索引为-1，则该行数据为空
            if answer_start_index == -1 + len('"answer"'):
                all_data[i] = None
            # 如果结束索引为-1，则获取该行数据的剩余部分
            elif answer_end_index == -1:
                all_data[i] = all_data[i][answer_start_index:]
            # 否则，获取该行数据的"answer"部分
            else:
                all_data[i] = all_data[i][answer_start_index:answer_end_index]
        # 如果包含"is"，则获取"is"的起始和结束索引
        elif "is" in all_data[i]:
            answer_start_index = all_data[i].find('is') + len('is')
            answer_end_index = all_data[i].find('.', answer_start_index)
            # 如果结束索引为-1，则获取该行数据的剩余部分
            if answer_end_index == -1:
                all_data[i] = all_data[i][answer_start_index:]
            # 否则，获取该行数据的"is"部分
            else:
                all_data[i] = all_data[i][answer_start_index:answer_end_index]
        # 如果可以转换为字典，则转换
        try:
            all_data[i] = eval(all_data[i])
        except:
            if all_data[i] != None:
                all_data[i] = all_data[i].replace('"', '')

        # 创建字典，存储答案和gold数据
        try:
            question_answer_dict = {"answer": all_data[i].replace("Theansweryoujudge", ""), "gold": gold_data[i]}
        except:
            question_answer_dict = {"answer": None, "gold": gold_data[i]}
        # 将字典添加到输出列表中
        output.append(question_answer_dict)

    # 将输出列表写入文件
    write_file(file_path, output)


# 生成问题抽取
def gen_extract_answer(file_path):
    # 从文件中读取数据
    all_data, gold_data = read_file(file_path)
    output = []
    for i in range(len(all_data)):
        # 处理答案为空的情况
        if all_data[i] == None:
            question_answer_dict = {"answer": all_data[i], "gold": gold_data[i]}
            output.append(question_answer_dict)
            continue
        # 判断是否是list  需要转成字符串
        if isinstance(all_data[i], list):
            all_data[i] = ".".join(all_data[i])
        all_data[i] = (all_data[i].replace("\n", "").replace("\\", "").replace("'", '"')
                       .replace("assistant", "answer").replace(": ", ":").replace("{\"answer\":\"{\"answer\":\"", "{\"answer\":\""))

        if '"answer":"' in all_data[i]:
            answer_start_index = all_data[i].find('"answer":"') + len('"answer":"')
            answer_end_index = all_data[i].find('"', answer_start_index)
            if answer_start_index != -1 + len('"answer":"') and answer_end_index == -1:
                all_data[i] = all_data[i][answer_start_index:]

            elif answer_start_index != -1 + len('"answer":"') and answer_end_index != -1:
                all_data[i] = all_data[i][answer_start_index:answer_end_index]

        elif '"answer":' in all_data[i]:
            answer_start_index = all_data[i].find('"answer":') + len('"answer":')
            if answer_start_index != -1 + len('"answer":'):
                all_data[i] = all_data[i][answer_start_index:]
        try:
            all_data[i] = str(eval(all_data[i]))
        except:
            all_data[i] = all_data[i].replace('"', '').replace("    ", "")
        # gold 处理
        if type(gold_data[i]) == list:
            gold_data[i] = '.'.join(gold_data[i])
        else:
            gold_data[i] = (gold_data[i].replace("['", "").replace("']", ""))
        question_answer_dict = {"answer": all_data[i], "gold": gold_data[i]}
        output.append(question_answer_dict)

    write_file(file_path, output)

# 写入jsonl
def write_file(file_path, output):
    with open(file_path.replace("json", "jsonl"), 'w', encoding='utf-8') as jsonl_file:
        for data_dict in output:
            # 将Python对象转换为JSON格式，并保存到文件中，便于数据交换和存储。
            json.dump(data_dict, jsonl_file, ensure_ascii=False)
            jsonl_file.write('\n')


def read_file(file_path):
    '''读取文件内容并返回答案和正确答案列表'''
    all_data = []
    gold_data = []
    with open(file_path, "r", encoding='utf-8') as file:
        '''打开文件'''
        for line in file:
            '''遍历文件每一行'''
            json_obj = json.loads(line)
            all_data.append(json_obj["answer"])
            if "gold" in json_obj.keys():
                gold_data.append(json_obj["gold"])
            elif "target" in json_obj.keys():
                gold_data.append(json_obj["target"])
            

    return all_data, gold_data

# 获取json文件路径名
def list_json_files(folder_path):
    file_list = []
    # 遍历文件夹
    for root, dirs, files in os.walk(folder_path):
        for filename in files:
            if filename.endswith(".json"):
                full_path = os.path.join(root, filename)
                file_list.append(full_path)
    return file_list

# def read_json_file(infile):



def split_jsonl_by_filename(input_file, output_folder):
    with open(input_file, 'r', encoding='utf-8') as infile:

        # zhengti = json.load(infile)

        filename = input_file.split("\\")[-1]
        # with open(infile,"r") as read:
        #     for line in read:
        #         data = json.loads(line)
        #         filename = data.get("filename")

        

        file_class = ""
        if "Chemical_Literature_Topic_Classification" in filename or   "BBBP" in filename or "ClinTox" in filename or "HIV" in filename or "判断任务" in filename or "True_False_Task" in filename or "选择任务" in filename or "Multiple_Choice" in filename or ("分子性质分类预测" in filename and "极性" in filename) or "分子性质分类预测极性" in filename or "Molecular_Property_Classification_Polarity" in input_file or "产率预测" in filename or "Product_Yield_Prediction" in filename or "Reaction_Type_Recognition_and_Induction" in filename or "化学文献主题分类" in filename:
            file_class = "classification\\"
        elif "ESOL" in filename or "HOMO" in filename or "Lipo" in filename or "LUMO" in filename or "分子性质回归预测极性" in filename or ("分子性质回归预测" in filename and  "极性" in filename) or "Molecular_Property_Regression_Polarity" in input_file or "熔点" in filename or "Melting_Point" in filename  or "沸点" in filename or "Boiling_Point" in filename or "合成难度评估" in filename  or "Synthetic_Difficulty_Evaluation" in filename or "反应温度推荐" in filename or "Reaction_Temperature_Recommendation" in filename or "反应时间推荐" in filename or "Reaction_Time_Recommendation" in filename:
            file_class = "regression\\"
        elif "SMILES_" in filename or "SMILES_Conversion" in filename or "SMILES与SELFIES互译" in filename or "SMILES转SELFIES" in filename or "2D分子图识别" in filename:
            file_class = "molecule_design_S\\"
        elif ("分子式_" in filename and "识别" not in filename ) or "Formula_test" in filename:
            file_class = "molecule_design_F\\"
        elif "IUPAC_" in filename or "IUPAC_Conversion" in filename:
            file_class = "molecule_design_I\\"
        elif "基于文本描述生成分子名称" in filename or "Molecular_Name_Generation_from_Text_Description" in filename:
            file_class = "molecule_design\\"
        elif "SIDER" in filename:
            file_class = "sider\\"
        elif "命名实体识别" in filename or "Chemical_Named_Entity_Recognition" in filename or "产率性能抽取" in filename or "Yield_Extraction" in filename or "催化类型抽取" in filename or "Catalysis_Type_Extraction" in filename or "化学反应类型识别归纳" in filename or "Reaction_Type_Recognition_and_Induction" in filename or "反应时间抽取" in filename or "Reaction_Time_Extraction" in filename or "合成反应添加剂抽取" in filename or "thetic_Reaction_Additive_Extraction" in filename or "反应温度抽取" in filename or "Reaction_Temperature_Extraction" in filename or "反应溶剂抽取" in filename or "Synthetic_Reaction_Solvent_Extraction" in filename  or "表征手段抽取" in filename or "Characterization_Method_Extraction" in filename:
            file_class = "entrty_extraction\\"
        elif "实体关系分类" in filename or "Chemical_Entity_Relationship_Classification" in filename:
            file_class = "relation_extraction\\"
        elif "合成反应产物抽取" in filename or "Reaction_Product_Extraction"  in filename or "合成反应底物抽取" in filename or "Synthetic_Reaction_Substrate_Extraction" in filename:
            file_class = "entity_recognition\\"
        elif "反应产物预测" in filename or "Reaction_Product_Prediction" in filename or "反应底物推荐" in filename or "Substrate_Recommendation" in filename or "配体推荐" in filename or "Ligand_Recommendation" in filename or "溶剂推荐" in filename or "Solvent_Recommendation" in filename or "试剂推荐" in filename or "Reagent_Recommendation" in filename  or "催化剂推荐" in filename or "Catalyst_Recommendation" in filename or "4.合成路径解析" in filename:
            file_class = "reagent_selection\\"
        elif "反应活化能刻画" in filename or "Reaction_Rate_Prediction" in filename:
            file_class = "Reaction_Rate_Prediction\\"
        elif "描述分子的物理化学性质" in filename  or "Physicochemical_Property_Prediction_from_Molecular_Structure" in filename or "填空任务" in filename or "Fill-in-the-Blank_Task" in filename or "简答任务" in filename or "Short_Answer_Task" in filename or "计算任务" in filename or "Calculation_Task" in filename or "合成路径推荐" in filename or "Synthetic_Pathway_Recommendation" in filename  or "反应中间体推导" in filename or "Intermediate_Derivation_test" in filename or "化学论文摘要生成" in filename  or "Chemical_Paper_Abstract_Generation" in filename  or "研究内容提纲生成" in filename or "Research_Outline_Generation" in filename or "名词解释" in filename or "图问答" in filename or "表问答" in filename or "谱问答" in filename or "面问答" in filename or "化学反应方程式识别" in filename or "化学反应方程式识别" in filename or "分子式识别" in filename:
            file_class = "Before_rating\\"

        # 将文件直接copy过去
        index = filename.rfind("\\") if filename.rfind("\\")!=-1 else 0
        filename = filename[index:]
        output = output_folder + file_class#.replace("D:\\ixf下载\\2024-12\\results\\机评", "D:\\ixf下载\\2024-12\\results\\机评返回")# + ".jsonl"
        os.makedirs(os.path.dirname(output), exist_ok=True)
        shutil.copy(input_file, output)



if __name__ == "__main__":
    # 文件夹路径,最多提取两级
    # 获取文件路径名[]
    file_list = list_json_files(r"D:\project\nips-chemeval\pipline\data\机评结果\按t拆分\GLM4v-3t")

    print(file_list,len(file_list))

    output_folder = r"D:\project\nips-chemeval\pipline\data\机评结果\按t拆分\GLM4v-3t\\output2_me\\"
    for file_path in file_list:
        if "original" not in file_path:
            split_jsonl_by_filename(file_path, output_folder)
    
    file_list = list_json_files(output_folder)
    for file_path in file_list:
        print(file_path)
        # 提取答案的主函数
        extract_answer(file_path)
        print(file_path, "over")