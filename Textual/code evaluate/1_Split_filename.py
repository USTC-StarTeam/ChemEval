import json
import os,re

### 可具体到模型
### You need to modify
dir_path = r'path/to/dataset/dir'

def split_filename(dir_path):
    for dir in os.listdir(dir_path):
        if dir.endswith('.py') or dir.endswith('.json'):
            continue

        dir_sub=os.path.join(dir_path,dir)
        for file in os.listdir(dir_sub):
            if not file.endswith('.json'):
                continue

            file_path=os.path.join(dir_sub,file)
            match=re.search(r'_(.*?)t.json',file)
            times=''
            if match:
                times=match.group(1)
            else:
                print(file+'------------没有t')

            output_dir=os.path.join(dir_sub,'{}t'.format(times))
            os.makedirs(output_dir, exist_ok=True)

            # 用于存储每个 task 对应的 JSON 对象
            task_data = {}

            # 读取输入 JSON 文件
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    data = json.loads(line)
                    filename=data["filename"]
                    task = filename

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
                # new_name = str(output_file_path.split("\\")[-1])
                with open(output_dir+"\\"+f"{task.replace('\\','')}.json", 'w', encoding='utf-8') as f:
                    for data in data_list:
                        f.write(json.dumps(data, ensure_ascii=False) + '\n')

            print("{}拆分完成！".format(file_path))

def split_filename_analysis(dir_path):
    for dir in os.listdir(dir_path):
        if dir.endswith('.py') or dir.endswith('.json'):
            continue
        dir_sub=os.path.join(dir_path,dir)
        output_dir=os.path.join(dir_sub,'1t')
        os.makedirs(output_dir, exist_ok=True)

        for file in os.listdir(dir_sub):
            if file.endswith('.json'):
                # 用于存储每个 task 对应的 JSON 对象
                task_data = {}

                # 读取输入 JSON 文件
                file_path=os.path.join(dir_sub,file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    for line in f:
                        data = json.loads(line)
                        if data['answer']=='None' or data['answer']=='null':
                            data['answer']=None
                        filename=data["task_type"]
                        data['filename']=filename
                        del data['source_file']
                        del data['task_type']
                        task = filename

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
                    # new_name = str(output_file_path.split("\\")[-1])
                    with open(output_dir+"\\"+f"{task.replace('\\','')}.json", 'w', encoding='utf-8') as f:
                        for data in data_list:
                            f.write(json.dumps(data, ensure_ascii=False) + '\n')

                print("{}拆分完成！".format(file_path))

if dir_path.endswith('times'):
    for dir in os.listdir(dir_path):
        split_filename(os.path.join(dir_path,dir))
elif dir_path.endswith('analysis'):
    split_filename_analysis(dir_path)
else:
    split_filename(dir_path)

