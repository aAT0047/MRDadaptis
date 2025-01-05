import os
import pandas as pd
from bayes_opt import BayesianOptimization
import concurrent.futures
from tqdm import tqdm
import subprocess
from evaf1 import evamain
from mulit_2lumpy import filter_vcf_by_evidence
from config import base_dir, input_dir, output_dir, vcf_dir, csv_file, log_folder, unique_output_folder, csv_path, last_vcf_dir, bam_dir, bayes_log, finalopt_dir
from concurrent.futures import ProcessPoolExecutor, as_completed
# 参数边界
# 定义参数范围
param_bounds = {
    "w": (1000, 1000000),
    "msw": (1, 6),
    "tt": (0, 50),
    "back_distance": (0, 150),
    "min_mapping_threshold": (0, 100),
    "min_clip": (10, 500),
    "read_length": (50, 200),
    "min_non_overlap": (50, 200),
    "discordant_z": (2,9)
}

class BayesianOptimizerF1:
    def __init__(self, data_path):
        self.data_path = data_path
        self.precision = 0.0 # 保存precision以便后续使用
        self.recall = 0.0       # 保存recall以便后续使用
    def load_data(self):
        """
        加载CSV数据，并提取优化需要的参数。
        """
        self.data = pd.read_csv(self.data_path)
        # print(self.data)
        self.data.columns = self.data.columns.str.strip()
        # print(self.data.columns)
        columns_to_keep = list(param_bounds.keys())

        # if 'w' not in self.data.columns:
        #     self.data['w'] = 10000

        self.X = self.data[columns_to_keep]

        self.y = self.data['f1_score']
     

        self.sample_name = self.data['sample_name'].iloc[0]
        self.mean = self.data['mean'].iloc[0]
        self.stdev = self.data['stdev'].iloc[0]

    def run_lumpy(self, params):
        """
        运行Lumpy并返回F1分数。
        """
        w = params["w"]
        msw = params["msw"]
        tt = params["tt"]
        back_distance = params["back_distance"]
        min_mapping_threshold = params["min_mapping_threshold"]
        min_clip = params["min_clip"]
        read_length = params["read_length"]
        min_non_overlap = params["min_non_overlap"]
        discordant_z = params["discordant_z"]

        output_vcf = os.path.join(vcf_dir, f"{self.sample_name}_optimized.vcf")
        
        cmd = [
            "lumpy",
            "-w", str(w),
            "-msw", str(msw),
            "-tt", str(tt),
            "-pe", f"id:{self.sample_name},read_group:readgroup1,bam_file:{bam_dir}/{self.sample_name}.discordants.bam,histo_file:{bam_dir}/{self.sample_name}.lib1.histo,mean:{self.mean},stdev:{self.stdev},read_length:{read_length},min_non_overlap:{min_non_overlap},discordant_z:{discordant_z},back_distance:{back_distance},weight:1,min_mapping_threshold:{min_mapping_threshold}",
            "-sr", f"id:{self.sample_name},bam_file:{bam_dir}/{self.sample_name}.splitters.bam,back_distance:{back_distance},weight:1,min_mapping_threshold:{min_mapping_threshold}",
            ">", output_vcf
        ]
        if not os.path.exists(bayes_log):
            os.makedirs(bayes_log)
        # 生成日志文件的路径
        log_file = os.path.join(bayes_log, f"{self.sample_name}_lumpy_log.txt")

        # 获取日志文件的目录路径（注意，只创建目录）
        log_dir = os.path.dirname(log_file)

        # 如果日志文件的父目录不存在，创建目录
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)  # 只创建目录，不是文件

        # 打开日志文件，并将命令的输出和错误日志写入文件
        with open(log_file, 'w') as log:
            # print(f"Logging to file {log_file}")
            process = subprocess.Popen(" ".join(cmd), shell=True, stdout=log, stderr=log)
            process.communicate()


        filter_vcf_by_evidence(output_vcf)
        output_vcf_path = os.path.join(vcf_dir, f"{self.sample_name}.vcf")
        TP, FP, FN, f1_score, precision, recall = evamain(output_vcf_path, output_vcf)
        # return f1_score
        return f1_score, precision, recall

    def objective(self, **params):
        # return self.run_lumpy(params)
        f1_score, precision, recall = self.run_lumpy(params)
        # 返回 f1_score 作为优化目标
        self.precision = precision  # 保存precision以便后续使用
        self.recall = recall        # 保存recall以便后续使用
        return f1_score  # 优化的目标仍然是f1_score

    def optimize(self):
        optimizer = BayesianOptimization(f=self.objective, pbounds=param_bounds, random_state=0)

        initial_points = self.X.to_dict(orient='records')
        # print(initial_points)
        initial_targets = self.y.tolist()
        # print(initial_points)
        # print(initial_targets)
        for point, target in zip(initial_points, initial_targets):
            optimizer.register(params=point, target=target)
            

        
        optimizer.maximize(init_points=1, n_iter=10)
        best_params = optimizer.max['params']
        best_target = optimizer.max['target']

        df_opt_params = pd.DataFrame([best_params])
        df_opt_params['f1'] = best_target
        best_precision = self.precision
        best_recall = self.recall

        df_opt_params['precision'] = best_precision
        df_opt_params['recall'] = best_recall

        return df_opt_params


    def run(self):
        self.load_data()
        return self.optimize()

    def save_optimized_parameters(self, save_path, sample_id, df_opt_params):
        df_opt_params.insert(0, 'sample_name', sample_id)
        df_opt_params.to_csv(save_path, index=False)


def optimize_for_sample(i, input_base_path, output_base_path):
    """
    处理单个样本的优化并保存结果。
    """
    input_path = os.path.join(input_base_path, f"{i}_params.csv")
    output_path = os.path.join(output_base_path, f"{i}_optparams.csv")

    optimizer = BayesianOptimizerF1(data_path=input_path)
    df_opt_params = optimizer.run()
    optimizer.save_optimized_parameters(output_path, i, df_opt_params)


def merge_optimized_csv_files(output_base_path, sample_numbers, merged_csv_path):
    """
    合并所有样本的优化结果为一个CSV文件。
    """
    dfs = []
    for sample_id in sample_numbers:
        file_path = os.path.join(output_base_path, f"{sample_id}_optparams.csv")
        if os.path.exists(file_path):
            dfs.append(pd.read_csv(file_path))
    
    if dfs:
        merged_df = pd.concat(dfs, ignore_index=True)
        merged_df.to_csv(merged_csv_path, index=False)
        print(f"Merged CSV saved to {merged_csv_path}")


def clean_up_individual_csv_files(output_base_path, sample_numbers):
    """
    删除单个样本的优化结果CSV文件。
    """
    for sample_id in sample_numbers:
        file_path = os.path.join(output_base_path, f"{sample_id}_optparams.csv")
        if os.path.exists(file_path):
            os.remove(file_path)
    print("All individual CSV files deleted.")


def run_optimization_parallel():
    """
    并行优化多个样本，并合并所有优化结果为一个文件。
    """
    input_base_path = csv_path
    output_base_path = finalopt_dir

    if not os.path.exists(output_base_path):
        os.makedirs(output_base_path)

    file_names = sorted(os.listdir(input_base_path))
    # sample_numbers = [name.split("_")[0] for name in file_names if name.endswith("_params.csv")]
    sample_numbers = [name.rsplit("_", 1)[0] for name in file_names if name.endswith("_params.csv")]
    with ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(optimize_for_sample, i, input_base_path, output_base_path): i for i in sample_numbers
        }

        # 使用 tqdm 包装 futures 的 as_completed 方法
        with tqdm(total=len(futures), desc="Optimizing samples") as pbar:
            for future in as_completed(futures):
                # try:
                future.result()  # 获取每个任务的结果，检查异常
                # except Exception as e:
                    # print(f"Error in sample {futures[future]}: {e}")
                pbar.update(1)  # 更新进度条

    merged_csv_path = os.path.join(output_base_path, "merged_dataf1.csv")
    merge_optimized_csv_files(output_base_path, sample_numbers, merged_csv_path)
    clean_up_individual_csv_files(output_base_path, sample_numbers)


def BayesianOptimizermain():
    """
    主函数：执行整个优化流程。
    """
    print("Starting optimization process...")
    
    # 初始化工作目录和配置
    # 创建输出目录（如果尚不存在）
    if not os.path.exists(finalopt_dir):
        os.makedirs(finalopt_dir)
    run_optimization_parallel()
   
    
    # # 开始并行优化过程
    # try:
    #     run_optimization_parallel()
    #     print("Optimization process completed successfully.")
    # except Exception as e:
    #     print(f"An error occurred during the optimization process: {e}")


# # 启动主函数
#   