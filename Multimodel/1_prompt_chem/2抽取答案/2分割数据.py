import json

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

data_list = read_file("../data/J-4379-gemini-gpt-claudet-qwen_export_all_fxx_gpt4o_0806_hsq_gpt4o_0806.json")
data_dic = {}
for data in data_list:
    if data['source_model'] not in data_dic.keys():
        data_dic[data['source_model']] = []
    data_dic[data['source_model']].append(data)


for key in data_dic.keys():
    write_file("../data/"+key.replace(".json","_result.json"),data_dic[key])
