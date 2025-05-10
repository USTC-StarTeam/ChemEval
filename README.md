# USAGE

## Textual

### Code Evaluate

The data after evaluating n times will be placed in the "dataset\_ntimes" folder, organized by model name, with the structure as follows:

```
dataset_ntimes/
├── model_1/
│   ├── model_1_nt.json
└── model_2/
    ├── model_2_nt.json
└── ...
```

#### STEP 1
Modify the `dir_path=path/to/dataset_ntimes` in `1_Split_filename.py` and run the file to split different task files based on the "filename" field.

#### STEP 2
Modify the `file_dir=path/to/dataset_ntimes` in `2_Extract.py` and run the file to extract the answers from the model responses and generate JSONL files, which will be placed in the `result_ntimes` folder with the following structure:

```
result_ntimes/
├── model_1/
│   ├── 0shot/
│   │   ├── 1t/
│   │   │   ├── task_dir_1/
│   │   │   │   ├──task_1.json
│   │   │   │   ├──task_1.jsonl
│   │   │   │   ├──task_2.json
│   │   │   │   ├──task_2.jsonl
│   │   │   │   ├──...
│   │   ├── 2t/
│   │   │   ├── task_dir_1/
│   │   │   │   ├──task_1.json
│   │   │   │   ├──task_1.jsonl
│   │   │   │   ├──...
│   ├── 3shot/
│   │   ├──...
└── ...
```

#### STEP 3
Modify the `root_path=path/to/result_ntimes` in `3_Evaluate.py` and run the file to obtain the evaluation results.

### LLM Evaluate
#### STEP 1
Modify the `folder_path` and `output_file` in `1_prompt_chem.py`, then run the file to submit for LLM evaluation.

#### STEP 2
Modify the `file_path` in `2_L1_task_eval.py` and run to obtain metrics for multiple-choice, true/false, fill-in-the-blank, short answer, and calculation tasks. The results will be displayed in an Excel file.

#### STEP 3
Modify the `foler_path` and `excel_path` in `3_other_task_eval.py`,then run the file to obtain metrics for abstract writing, outlining, reaction intermediates, single-step synthesis, multi-step synthesis, and physicochemical property tasks. The results will be displayed in an Excel file.

## Multimodal
### STEP 1
Run `1_prompt_chem` to obtain the evaluation data from the LLM.

### STEP 2
Run `2_LLM_evaluate` to obtain the evaluation results from the LLM.

### STEP 3
Run `3_code_evaluate` to obtain the evaluation results using code.

# Data

https://huggingface.co/datasets/Ooo1/ChemEval

# Licenses

[![CC BY-NC-SA 4.0](https://camo.githubusercontent.com/f61dcd7e9460d79b9e8e19683c964e21cc2455ff9d8860cc5ca30f35457be635/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f4c6963656e73652d434325323042592d2d4e432d2d5341253230342e302d6c69676874677265792e737667)](http://creativecommons.org/licenses/by-nc-sa/4.0/)

The ChemEval dataset is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/)

# Citation

Please cite our paper if you use our dataset.

```
@article{huang2024chemeval,
  title={ChemEval: A Comprehensive Multi-Level Chemical Evaluation for Large Language Models},
  author={Huang, Yuqing and Zhang, Rongyang and He, Xuesong and Zhi, Xuyang and Wang, Hao and Li, Xin and Xu, Feiyang and Liu, Deguang and Liang, Huadong and Li, Yi and others},
  journal={arXiv preprint arXiv:2409.13989},
  year={2024}
}
```