# !!!注意!!! 
# 遇到 genome 中有重複基因時，此 script 會出問題，重複的序列範圍會多出來在該行最後
# 為避免發生錯誤，請在使用完成後手動檢查 outfile 並利用基因範圍比對 gb file

# 使用前請先 cd 至目標 folder
# 請先統計好共同的基因並以一列一個 gene 的格式存入 txt 中
# Usage: python script.py <gbfile_path> <comparison_file_path> <target_genes_list_path>
# eg. python getGeneInfo_final.py Chloroidium_v2.gb Jaagichlorella_hainangensis_v2.gb geneList.txt
# 最終完成資料的擷取後會輸出成 final_output.txt

import sys
import re

def read_file(file_path):
    try:
        with open(file_path, "r") as file:
            data = file.read()
        return data
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        sys.exit(1)

def read_target_genes(file_path):
    data = read_file(file_path)
    target_genes = list(data.split("\n"))
    return target_genes

def buildRule(rule):
    pattern = re.compile(rule)
    return pattern

def write_to_file(output_file, gene_name, start, end, constants):
    with open(output_file, "a") as output:
        output.write(f"{gene_name}\t{constants['fileName1']}\t{constants['fileName2']}\t{constants['constant1']:.3f}\t{constants['length']}\t{constants['constant2']}\t{constants['constant3']}\t{start}\t{end}\n")

#取第一部分基因
def extract_gene_info(data, target_gene_names, rule, output_file, constants):
    pattern = buildRule(rule)
    matches = pattern.finditer(data)
    times = 0

    for match in matches:
        start = match.group("start")
        end = match.group("end")
        gene_name = match.group("gene")
        #print(gene_name)
        #length = int(end) - int(start) + 1

        if gene_name in target_gene_names: #scan gene list user input
            times += 1 # 最後可用來確認基因個數是否正確
            # if "complement" in match.group(0):  # 如果是complement基因
            #     start, end = end, start # 交換start和end
            print(f"{times}.gene: {gene_name}")
            # 將結果寫入文字檔
            write_to_file(output_file, gene_name, start, end, constants)

    print(f"gene amount : {times}")

def append_gene_info(existing_file, new_data, rule, target_gene_names, constants):
    pattern = buildRule(rule)
    matches = pattern.finditer(new_data)
    times = 0
    # 讀取現有的檔案
    with open(existing_file, "r") as existing:
        lines = existing.readlines()

    for match in matches:
        start = match.group("start")
        end = match.group("end")
        gene_name = match.group("gene")

        # 如果基因在使用者提供的列表中，則附加新的 start 和 end
        if gene_name in target_gene_names:
            
            # 避免重複搜尋 eg. ycf1 搜尋到 ycf12
            gene_name = gene_name + "\t"

            # 在相同基因的行尾附加新的 start 和 end
            for i, line in enumerate(lines):
                # print(i) # debug 
                # print(line) # debug
                
                if gene_name in line: 
                    # 取得 length，這裡已經改成由第二個檔案的 length
                    length = int(end) - int(start) + 1
                    # times += 1
                    # print(times + ". " + gene_name)
                    # 判斷是否為 complement
                    if "complement" in match.group(0):
                        start, end = end, start

                    # 將 -1 取代為實際的 length
                    line_new = line.replace("-1", str(length))

                    # 附加新的 start 和 end 以及 constants
                    line_new = line_new.rstrip() + f"\t{start}\t{end}\t{constants['constant4']}\t{constants['constant5']}\n"
                    # print(line) # debug
                    # 更新行
                    lines[i] = line_new


    # 覆寫整個檔案
    with open(existing_file, "w") as output:
        output.writelines(lines)
    print("execution completed")

def remove_gene_names(input_file):
    with open(input_file, "r") as input_data:
        lines = input_data.readlines()

    # 去除基因名稱
    new_lines = [line.split("\t")[1:] for line in lines]

    with open(input_file, "w") as output:
        # 寫入新檔案（以覆蓋的形式）
        for line in new_lines:
            output.write("\t".join(line))

def clear_file(file_path):
    with open(file_path, "w") as file:
        file.truncate(0)

def instruction():
    print("")
    print("arguments or python version error")
    print("plz read following illustration and try again")
    print("")
    print("--------------------------------------------------------")
    print("This is a tool which can parse specific genes in GeneBank via python3")
    print("Please go over to the folder which we want put result in (use cd change the dir of file )")
    print("Please provide a gb file as argument 1 (both .txt and .gb file are acceptable)")
    print("comparison file as argument 2 (both txt and gb file are acceptable)")
    print("target gene txt file as argument 3 (one row can only contain one gene)")
    print("Usage: python script.py <file_path> <comparison_file_path> <target_genes_path>")
    print("--------------------------------------------------------")
    print("")

def main():
    
    # 確保終端提供了檔案路徑
    if len(sys.argv) != 4:
        # 程式簡介
        instruction()
        sys.exit(1)
    clear_file("output.txt") 
    # 從終端參數中讀取檔案路徑
    file_path = sys.argv[1]
    comparsion_file_path = sys.argv[2]
    target_genes_path = sys.argv[3]

    # 讀取檔案
    geneBankData = read_file(file_path)

    # 定義一個要擷取的基因名稱列表
    target_gene_names = read_target_genes(target_genes_path)

    # 指定輸出文字檔的路徑
    output_file = "output.txt"

    # 寫入標頭
    # write_header(output_file)

    # 定義常數值
    constants = {
        'fileName1': "1.easyfig.fa",
        'fileName2': "2.easyfig.fa",
        'constant1': 75.000,
        'length': -1,
        'constant2': 0,
        'constant3': 0,
        'constant4': 0.0,
        'constant5': 75.2
    }

    # 建立要保留之資訊的regular expression規則 
    complementRule = r'gene\s+(?:complement\()?(?P<start>\d+)\.\.(?P<end>\d+)\)?\s+/gene="(?P<gene>\w+)"'

    # 提取並輸出基因資訊，先輸出正股再反股
    extract_gene_info(geneBankData, target_gene_names, complementRule, output_file, constants)
    
    # 讀取新檔案
    comparsion_data = read_file(comparsion_file_path)
    
    append_gene_info(output_file, comparsion_data, complementRule, target_gene_names, constants)
    remove_gene_names("output.txt")
    # clear_file("output.txt")

    # remove_header(output_file)

if __name__ == "__main__":
    main()
