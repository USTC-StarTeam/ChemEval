'''
Author: jiancui3
Date: 2025-02-26 13:33:17
LastEditTime: 2025-03-13 21:34:35
FilePath: \测评\化学新版benchmark\code-0307\chemeval_code\Evaluate.py
'''
from Reaction_Rate_Prediction import reaction_rate
from classification import classify
from regression import reg
from molecule_design import design
from molecule_design_F import design_F
from molecule_design_I import design_I
from molecule_design_S import design_S
from entity_extraction import entity_extract
from entity_recognition import recognize
from reagent_selection import select
from sider import sider
from relation_extraction import relation_extract
from jiping import jiping
import os

if __name__ == "__main__":
    model_list = ["claude3_7t-4-1t", "gemini2.5pro-4-1t",
                  "gpt4o-4-1t", "qwen-vl-max-4-1t"]
    root_path = r"D:\project\nips-chemeval\pipline\data\机评结果\按t拆分\glm4v-4-3times-prompted_res-3t\output2_me"
    excel_path = r"D:\project\nips-chemeval\pipline\data\机评结果\按t拆分\glm4v-4-3times-prompted_res-3t\test.xlsx"

    print("========Before_rating========")
    if not any("Before_rating" in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'Before_rating' found in the root directory")
    else:
        jiping(root_path)
    print("========reaction_rate========")
    if not any("reaction_rate" in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'reaction_rate' found in the root directory")
    else:
        reaction_rate(root_path, excel_path)

    print("=======classify=========")
    if not any("classification" in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'classify' found in the root directory")
    else:
        classify(root_path,excel_path)

    print("=======sider=========")
    if not any("sider" in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'sider' found in the root directory")
    else:
        sider(root_path, excel_path)

    print("========reg========")
    if not any("regression" in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'reg' found in the root directory")
    else:
        reg(root_path)

    print("=======design=========")
    if not any("molecule_design" in folder and "molecule_design_" not in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'design' found in the root directory")
    else:
        design(root_path, excel_path)

    print("=======design_F=========")
    if not any("molecule_design_F" in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'design_F' found in the root directory")
    else:
        # pass
        design_F(root_path)

    print("=======design_I=========")
    if not any("molecule_design_I" in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'design_I' found in the root directory")
    else:
        # pass
        design_I(root_path)

    print("=======design_S=========")
    if not any("molecule_design_S" in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'design_S' found in the root directory")
    else:
        # pass
        design_S(root_path)

    print("=======entity_extract=========")
    if not any("entrty_extraction" in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'entity_extract' found in the root directory")
    else:
        entity_extract(root_path, excel_path)

    print("=======recognize========")
    if not any("entity_recognition" in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'recognize' found in the root directory")
    else:
        recognize(root_path, excel_path)

    print("=======select=========")
    if not any("reagent_selection" in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'select' found in the root directory")
    else:
        select(root_path)

    print("=======relation_extract=========")
    if not any("relation_extraction" in folder for folder in os.listdir(root_path)):
        print("Warning: No subfolder containing 'relation_extract' found in the root directory")
    else:
        relation_extract(root_path, excel_path)
    
    # 评估任务 = ["选择任务","填空任务","判断任务","简答任务","计算任务","化学命名实体识别","化学实体关系分类","合成反应底物抽取","合成反应添加剂抽取","合成反应溶剂抽取","合成反应温度抽取","合成反应时间抽取","合成反应产物抽取","表征手段抽取","催化类型抽取(热催化，电催化，光催化)","产率性能抽取","化学论文摘要生成（含方法、结论和影响力）","研究内容提纲生成","化学文献主题分类","化学反应类型识别归纳(如加成反应、消除反应、取代反应等)","基于文本描述生成分子名称","IUPAC转分子式","SMILES转分子式","IUPAC转SMILES","SMILES转IUPAC","SMILES与SELFIES互译","分子性质分类","分子性质回归","基于分子结构描述分子的物理化学性质","反应底物推荐","合成路径推荐","合成难度评估","配体推荐","试剂推荐","溶剂推荐","催化剂推荐","反应温度推荐","反应时间推荐","反应产物预测","产物产率预测","反应速率预测","反应中间体推导"]

    # 评估任务 = ["选择任务","填空任务","判断任务","简答任务","计算任务","化学命名实体识别","化学实体关系分类","合成反应底物抽取","合成反应添加剂抽取","合成反应溶剂抽取","合成反应温度抽取","合成反应时间抽取","合成反应产物抽取","表征手段抽取","催化类型抽取","产率性能抽取","化学论文摘要生成","研究内容提纲生成","化学文献主题分类","化学反应类型识别归纳","基于文本描述生成分子名称","IUPAC转分子式","SMILES转分子式","IUPAC转SMILES","SMILES转IUPAC","SMILES与SELFIES互译","分子性质分类","分子性质回归","基于分子结构描述分子的物理化学性质","反应底物推荐","合成路径推荐","合成难度评估","配体推荐","试剂推荐","溶剂推荐","催化剂推荐","反应温度推荐","反应时间推荐","反应产物预测","产物产率预测","反应速率预测","反应中间体推导"]

    # import pandas as pd
    # xls = pd.ExcelFile(excel_path)
    
    # for i, task in enumerate(评估任务):
    #     # 遍历每个 Sheet
    #     for sheet_name in xls.sheet_names:
    #         # print(f"Reading Sheet: {sheet_name}")
            
    #         # 读取当前 Sheet 的数据
    #         df = pd.read_excel(xls, sheet_name=sheet_name)
            
    #         # 遍历 DataFrame 的每一行
    #         for index, row in df.iterrows():
    #             # print(f"Row {index + 1}: {row.tolist()}")  # 转换为列表打印

    #             if task in row["File Name"]:
    #                 print(task)
    #                 break
    #             else:
    #                 print("\n")
                
            

