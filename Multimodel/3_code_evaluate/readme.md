### 化学任务自动化测评流程分为两步：  
1. 运行Extract.py格式化json文件并抽取预测值和标签值  
2. 运行Evaluate.py获得每种任务的输出指标

### 您需要做如下修改：  
- 在Extract.py文件中：  
1. file_list = list_json_files("")  # 修改为您的根目录  
output_folder = "Llama3kale\\"  # 修改为输出文件目录  
2. 在split_jsonl_by_filename(input_file, output_folder)函数中，new_data字典中的"target"键需要给你的实际标签值  

- 在Evaluate.py文件中：  
root_path = "" # 修改为之前output_folder地址
