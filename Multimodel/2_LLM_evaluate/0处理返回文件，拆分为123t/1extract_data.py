import json
import re

res_list = []
data_list = []
model_list = ["claude3_7t-4-3times-prompted_result","gemini2.5pro-4-3times-prompted_result",
              "gpt4o-4-3times-prompted_result","qwen-vl-max-4-3times-prompted_result","glm4v-4-3times-prompted_res","gpt4o-3times"]
model_name = model_list[5]+".json"

input_file = "../../data/"+model_name
output_file_1 = "../../data/机评结果/按t拆分/"+ model_name.replace(".json","-1.json")
output_file_2 = "../../data/机评结果/按t拆分/"+ model_name.replace(".json","-2.json")
output_file_3 = "../../data/机评结果/按t拆分/"+ model_name.replace(".json","-3.json")


with open(input_file, 'r', encoding="utf-8") as f:
    for line in f:
        data = json.loads(line)
        data_list.append(data)

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
# 定义正则表达式
keys = ["<最终结果>", "<抽取结果>", "<最终选项>", "<最终得分>"]
pattern = r"<(.*?)>\n([^\n]*)"


def extract(input_str):
    # 提取字符串
    extracted_data = {}
    for key in keys:
        match = re.search(rf"{key}\**\s*(.*?)(?=\n<|$)", input_str, re.DOTALL)
        if match:
            extracted_data[key] = match.group(1).strip()
    res = {}
    if len(extracted_data.items()) == 0:
        return False
    for key, value in extracted_data.items():
        res[key] = value
    return res


def extract_chemical_and_values(text):
    if text in ['None','none']:
        return "None"
    # 提取JSON格式的反应物、条件和产物
    json_pattern = r'{"反应物": "(.*?)", "条件": "(.*?)", "产物": "(.*?)"}'
    json_matches = re.findall(json_pattern, text)
    if json_matches:
        return json_matches[0]
    # 提取分隔符后的数值（例如 +0.31, -50 等）
    value_pattern = r"([+-]?\d*\.\d+|\d+)(?=\s*([A-Za-z\(\)°\%\^]*|\s*\w+\s*|\s*MPa\^1/2)?)"
    value_matches = re.findall(value_pattern, text)
    if value_matches:
        return value_matches[0]
    # 提取包含化学式的字符串（例如 C16H28INO5Si, CCC(CC1=CC=C(C=C1)SC)N.Cl）
    chemical_pattern = r"[A-Za-z0-9\(\)\=\.\-\+]+"
    chemical_matches = re.findall(chemical_pattern, text)
    if chemical_matches:
        return chemical_matches[0]
    return False

# data = """<抽取结果>2
#  \text{Cu} + \text{O}_2 + \text{CO}_2 + \text{H}_2\text{O} \rightarrow \text{Cu}_2(\text{OH})_2\text{CO}_3
# <最终结果>
# 正确
# """
#
# print(extract(data))

for index,data in enumerate(data_list):
    res = extract(data['answer'])
    if res:
        if "<最终选项>" in res.keys():
            res = {
                "answer":res['<最终选项>']
            }
        elif "<最终得分>" in res.keys():
            res = {
                "answer":res["<最终得分>"]
            }
        elif "<最终结果>" in res.keys():
            res = {
                "answer": res["<最终结果>"]
            }
        else:
            res = {
                "answer": res["<抽取结果>"]
            }

        res['gold'] = data['gold']
        res['task_type'] = data['task_type']
        res_list.append(res)
    else:
        t = extract_chemical_and_values(data['answer'])
        if t:
            res_list.append({
                "answer": t,
                'gold':data['gold'],
                "task_type":data['task_type']
            })
        else:
            print(index,data['answer'])
            res_list.append({
                "answer": None,
                'gold': data['gold'],
                "task_type": data['task_type']
            })

# write_file("../../data/机评结果/qwen-vl-max-3times.json",res_list)

length = len(res_list)
res1 = []
res2 = []
res3 = []

if  length % 3 != 0:
    raise Exception


def deal_data_num(data_list):
    data_dic = {}
    for data in data_list:
        if data['task_type'] not in data_dic.keys():
            data_dic[data['task_type']] = []
        data_dic[data['task_type']].append(data)
    whole_res = []
    for key in data_dic.keys():
        length = len(data_dic[key])
        length = length - length%10
        data_dic[key] = data_dic[key][:length]
        whole_res.extend(data_dic[key])
        print(key,len(data_dic[key]))
    return whole_res



index = 0
while index<length:
    res1.append(res_list[index])
    res2.append(res_list[index+1])
    res3.append(res_list[index+2])
    index+=3

res1 = deal_data_num(res1)
res2 = deal_data_num(res2)
res3 = deal_data_num(res3)

# print(len)
write_file(output_file_1, res1)
write_file(output_file_2, res2)
write_file(output_file_3, res3)

deal_data_num(res1)
