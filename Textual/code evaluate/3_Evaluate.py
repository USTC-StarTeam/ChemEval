'''
Author: jiancui3
Date: 2025-02-26 13:33:17
LastEditTime: 2025-03-13 21:34:35
FilePath: \测评\化学新版benchmark\code-0307\chemeval_code\Evaluate.py
'''
from Reaction_Rate_Prediction import reaction_rate
from classification import classify
from regression_new import reg
from molecule_design_F import design_F
from molecule_design_I import design_I
from molecule_design_S import design_S
from entity_extraction import entity_extract
from entity_recognition import recognize
from reagent_selection import select
from sider import sider
from relation_extraction import relation_extract
import os

if __name__ == "__main__":
    ### 可具体到t
    ### You need to modify
    root_path = r"path/to/result/dir"
    
    if root_path.endswith('times') or root_path.endswith('analysis'):
        base_depth = root_path.count(os.sep)
        for root, dirs, files in os.walk(root_path):
            if root.endswith('t') and not root.endswith('shot'):
            #if not root.endswith('t') and not root.endswith('shot') and not root.endswith('times'):
                print(root+'='*100)
                print("========reaction_rate========")
                reaction_rate(root)
                print("=======classify=========")
                classify(root)
                sider(root)
                print("========reg========")
                reg(root)
                print("=======design_F=========")
                design_F(root)
                print("=======design_I=========")
                design_I(root)
                print("=======design&design_S=========")
                design_S(root)
                print("=======entity_extract=========")
                entity_extract(root)
                print("=======recognize========")
                recognize(root)
                print("=======select=========")
                select(root)
                print("=======relation_extract=========")
                relation_extract(root)
    elif not root_path.endswith('t') and not root_path.endswith('shot'):
        for root, dirs, files in os.walk(root_path):
            if root.endswith('t') and not root.endswith('shot'):
                print(root+'='*100)
                print("========reaction_rate========")
                reaction_rate(root)
                print("=======classify=========")
                classify(root)
                sider(root)
                print("========reg========")
                reg(root)
                print("=======design_F=========")
                design_F(root)
                print("=======design_I=========")
                design_I(root)
                print("=======design&design_S=========")
                design_S(root)
                print("=======entity_extract=========")
                entity_extract(root)
                print("=======recognize========")
                recognize(root)
                print("=======select=========")
                select(root)
                print("=======relation_extract=========")
                relation_extract(root)
    elif root_path.endswith('t') and not root_path.endswith('shot'):
        print("========reaction_rate========")
        reaction_rate(root_path)
        print("=======classify=========")
        classify(root_path)
        sider(root_path)
        print("========reg========")
        reg(root_path)
        print("=======design_F=========")
        design_F(root_path)
        print("=======design_I=========")
        design_I(root_path)
        print("=======design&design_S=========")
        design_S(root_path)
        print("=======entity_extract=========")
        entity_extract(root_path)
        print("=======recognize========")
        recognize(root_path)
        print("=======select=========")
        select(root_path)
        print("=======relation_extract=========")
        relation_extract(root_path)

   
