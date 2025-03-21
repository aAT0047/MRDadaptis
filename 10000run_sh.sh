#!/bin/bash

# 设置工作目录
work_dir="/home/cloudam/GSDcreator/sh_files"

# 获取.sh脚本总数
script_count=100

# 使用tqdm来显示for循环的进度
echo $(seq 1 $script_count) | tqdm --unit "script" --total $script_count | while read i; do
    script_name="A_stableCallerPaperSimFlowShell_${i}.sh"
    script_path="${work_dir}/${script_name}"

    # 检查脚本是否存在
    if [ -f "$script_path" ]; then
        # 使用bash来运行.sh脚本
        cmd="bash ${script_path}"
        $cmd
        if [ $? -eq 0 ]; then
            echo "Successfully executed ${script_name}"
        else
            echo "Error in executing ${script_name}"
        fi
    else
        echo "${script_name} does not exist!"
    fi
done

echo "All scripts executed!"
