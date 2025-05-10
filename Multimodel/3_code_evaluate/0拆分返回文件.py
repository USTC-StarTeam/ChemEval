import json
import os

model_list = ["claude3_7t-4-3times-prompted_result","gemini2.5pro-4-3times-prompted_result",
              "gpt4o-4-3times-prompted_result","qwen-vl-max-4-3times-prompted_result"]

# 输入 JSON 文件路径
input_json_path = "../data/机评结果/按t拆分/GLM4v-3times-3.json"
# 输出目录，用于存放按 task 拆分的 JSON 文件
output_dir = r"../data/机评结果/按t拆分/GLM4v-3t/"

# 如果输出目录不存在，则创建它
os.makedirs(output_dir, exist_ok=True)

# 用于存储每个 task 对应的 JSON 对象
task_data = {}

# 读取输入 JSON 文件
with open(input_json_path, 'r', encoding='utf-8') as f:
    for line in f:
        data = json.loads(line)
        task = "\\".join(data.get("task_type").split("\\")[-3:])

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
    new_name = str(str(".".join(output_file_path.split("\\")[-5:])))
    # print(new_name)
    with open(new_name, 'w', encoding='utf-8') as f:
        for data in data_list:
            f.write(json.dumps(data, ensure_ascii=False) + '\n')

print("拆分完成！")