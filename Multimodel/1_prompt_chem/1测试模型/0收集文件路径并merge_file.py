import os
import json
from pathlib import Path
from statistics import median_grouped
import os
import shutil

import langid
def has_chinese_langid(datas):
    for data in datas:
        text = data["query"]
        """检测文本语言是否为中文"""
        lang, _ = langid.classify(text)
        if lang == 'zh':
            print(True)
            print(text)
            return True
    print(False)
    return False

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

def collect_json_paths(root_folder):
    # 递归遍历所有JSON和JSONL文件路径（新增.jsonl匹配）
    json_files = list(Path(root_folder).glob("**/*.json")) + list(Path(root_folder).glob("**/*.jsonl"))  # [8](@ref)
    # 转换为字符串路径列表（绝对路径）
    file_paths = [str(file.resolve()) for file in json_files]
    # 过滤包含特定字符的路径（新增过滤条件）[1,8](@ref)
    exclude_patterns = [
        # "1.领域知识问答",
        "2.文献理解",
        # "3.分子图文理解",
        # "4.科学知识推演"
    ]
    file_paths = [
        file for file in file_paths
        if not any(pattern in file for pattern in exclude_patterns)
    ]
    return file_paths

def copy_files_to_target(source_folders, target_folder):
    """
    将多个文件夹中的文件统一复制到一个目标文件夹中

    Args:
        source_folders (list): 包含源文件夹路径的列表
        target_folder (str): 目标文件夹路径
    """
    # 确保目标文件夹存在
    os.makedirs(target_folder, exist_ok=True)

    for folder in source_folders:
        print(folder)
        # 遍历每个源文件夹
        for root, _, files in os.walk(folder):
            print(files)
            for file in files:
                source_path = os.path.join(root, file)
                target_path = os.path.join(target_folder, file)

                # 检查目标文件夹中是否已经存在同名文件
                if os.path.exists(target_path):
                    # 如果存在，重命名文件以避免覆盖
                    base, extension = os.path.splitext(file)
                    counter = 1
                    while os.path.exists(target_path):
                        new_filename = f"{base}_{counter}{extension}"
                        target_path = os.path.join(target_folder, new_filename)
                        counter += 1

                # 复制文件
                shutil.copy2(source_path, target_path)
                print(f"Copied: {source_path} -> {target_path}")


root_folder = "../../../data/t1/化学多模态数据集/"  # 替换为实际文件夹路径
path_file = "../data/文件路径-1_3_4.json"  # 输出文件路径
merge_file = "../data/merge-4.json"
file_paths = collect_json_paths(root_folder)
# file_paths = [file for file in file_paths if has_chinese_langid(read_file(file))]


# 确保输出目录存在
output_dir = os.path.dirname(path_file)
if output_dir:
    os.makedirs(output_dir, exist_ok=True)
write_file(path_file, file_paths, pattern=2)
print(f"已保存{len(file_paths)}个JSON文件路径至：{path_file}")

merge_list = []
for file in file_paths:
    a  = read_file(file)
    for data in a:
        data['file_path'] = file
        if type(data["img_path"]) == str:
            data["img_name"] = data["img_path"].split("/")[-1]
        else:
            data["img_name"] = []
            for path in data["img_path"]:
                data["img_name"].append(path.split("/")[-1])
    merge_list.extend(a)
write_file(merge_file,merge_list)

img_target_folder = "../data/images"
img_path_list = []
for file in file_paths:
    # 目标文件夹名称
    target_folder_name = "images"
    # 提取文件夹路径
    folder_path = os.path.dirname(file)
    # 将文件夹路径与目标文件夹名称组合
    target_folder_path = os.path.join(folder_path, target_folder_name)
    img_path_list.append(target_folder_path)
copy_files_to_target(img_path_list,img_target_folder)


