import json
import os


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

# 输入 JSON 文件路径
input_json_path = r"c:\Users\xiangli67\Downloads\J-3719-gpt4o0806_70b_exp6_6_3_4205条_化工大模型_export_all_hsq_gpt4o_0806_fxx_gpt4o_0806.json"

# 输出目录，用于存放按 task 拆分的 JSON 文件
output_dir = r"C:\Users\xiangli67\Downloads\J3719"

# 如果输出目录不存在，则创建它
os.makedirs(output_dir, exist_ok=True)

# 用于存储每个 task 对应的 JSON 对象
task_data = {}

# 读取输入 JSON 文件
with open(input_json_path, 'r', encoding='utf-8') as f:
    for line in f:
        data = json.loads(line)
        task = data.get("filename").split("\\")[3]
        # task = data.get("filename").split("\\0304\\")[1]
        # task = data.get("task")
        

        if task:
            # 如果 task 不在字典中，则初始化一个空列表
            if task not in task_data:
                task_data[task] = []
            # 将当前 JSON 对象添加到对应 task 的列表中
            task_data[task].append(data)

# 将每个 task 的 JSON 对象写入对应的文件
for task, data_list in task_data.items():
    output_file_path = os.path.join(output_dir, f"{task}.json")
    # new_name = str(output_file_path.split("\\")[1]).replace("\\","")
    new_name = str(output_file_path.split("\\")[-1])
    with open(output_dir+"\\"+new_name, 'w', encoding='utf-8') as f:
        for data in data_list:
            f.write(json.dumps(data, ensure_ascii=False) + '\n')

print("拆分完成！")