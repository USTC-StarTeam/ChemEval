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

file_list = ["../../data/gemini2.5pro-4-3times-prompted.json","../../data/gpt4o-4-3times-prompted.json","../../data/claude3_7t-4-3times-prompted.json","../../data/qwen-vl-max-4-3times-prompted.json"]
data_whole = []
for file in file_list:
    datas = read_file(file)
    for data in datas:
        data['source_model'] = file.split("/")[-1]
    data_whole.extend(datas)

write_file("../../data/gemini-gpt-claudet-qwen.json",data_whole)