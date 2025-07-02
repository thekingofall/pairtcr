#!/usr/bin/env python3
import pandas as pd
import argparse
import sys
import os
import re
from tqdm import tqdm

# --- 配置 ---
DEFAULT_UMI_PAIRS_FILE = os.path.join("PairTCR_results", "2_create_umi_pairs_output", "umi_pairs.tsv")
DEFAULT_TRA_EXPORT_FILE = os.path.join("PairTCR_results", "3_run_mixcr_and_export_output", "TRA_alignments_export_with_headers.tsv")
DEFAULT_TRB_EXPORT_FILE = os.path.join("PairTCR_results", "3_run_mixcr_and_export_output", "TRB_alignments_export_with_headers.tsv")
DEFAULT_OUTPUT_FILE = os.path.join("PairTCR_results", "4_pair_and_filter_clones_output", "final_paired_clones_filtered.tsv") # Changed output name

# --- 辅助函数 ---
def get_base_read_id(header):
    """
    从 FASTQ 头部字符串提取基础读段 ID。
    (去除 @, 可能的 UMI 部分, 以及 /1 或 /2)
    """
    if not header or not isinstance(header, str):
        return None
    base_id_part = header.split(maxsplit=1)[0]
    if base_id_part.startswith('@'):
        base_id_part = base_id_part[1:]
    base_id_part = re.sub(r'/[12]$', '', base_id_part)
    return base_id_part

def perform_pairing(umi_pairs_file, tra_export_file, trb_export_file, output_file):
    """
    读取输入文件, 过滤交叉比对, 并执行配对。

    Args:
        umi_pairs_file (str): UMI 配对文件路径。
        tra_export_file (str): TRA 导出文件路径。
        trb_export_file (str): TRB 导出文件路径。
        output_file (str): 输出配对结果的文件路径。
    """
    print("--- 开始执行外部配对 (带过滤) ---")

    # --- 1. 输入文件验证 ---
    required_files = [umi_pairs_file, tra_export_file, trb_export_file]
    missing_files = []
    for f in required_files:
        if not os.path.exists(f):
            missing_files.append(f)
    if missing_files:
        print(f"错误: 必需的输入文件未找到: {', '.join(missing_files)}", file=sys.stderr)
        sys.exit(1)

    # --- 2. 读取 UMI 配对文件 ---
    try:
        print(f"读取 UMI 配对文件: {umi_pairs_file}")
        umi_pairs_df = pd.read_csv(umi_pairs_file, sep='\t', usecols=['TRA_Read_ID_Base', 'TRB_Read_ID_Base'])
        umi_pairs_df = umi_pairs_df.drop_duplicates().reset_index(drop=True)
        print(f"读取了 {len(umi_pairs_df)} 条唯一的 TRA-TRB Read ID 配对关系。")
    except Exception as e:
        print(f"错误: 读取或处理 UMI 配对文件 '{umi_pairs_file}' 时出错: {e}", file=sys.stderr)
        print("请确保文件是制表符分隔，并且包含 'TRA_Read_ID_Base' 和 'TRB_Read_ID_Base' 列。", file=sys.stderr)
        sys.exit(1)

    # --- 3. 读取、过滤并处理 TRA 导出文件 ---
    try:
        print(f"读取 TRA 导出文件: {tra_export_file}")
        tra_df = pd.read_csv(tra_export_file, sep='\t', low_memory=False) # low_memory=False may help with mixed types
        print(f"读取了 {len(tra_df)} 条 TRA 原始比对记录。")

        # *** 新增过滤步骤 ***
        # 确保基因列是字符串类型，处理可能的 NaN 值
        tra_df['bestVGene'] = tra_df['bestVGene'].astype(str)
        tra_df['bestJGene'] = tra_df['bestJGene'].astype(str)
        # 过滤：只保留 V 或 J 基因以 'TRA' 开头的行
        tra_df_filtered = tra_df[
            tra_df['bestVGene'].str.startswith('TRA', na=False) |
            tra_df['bestJGene'].str.startswith('TRA', na=False)
        ].copy() # 使用 .copy() 避免 SettingWithCopyWarning
        print(f"过滤后剩下 {len(tra_df_filtered)} 条 TRA 比对记录。")
        del tra_df # 释放内存

        # 从 descrsR1 列提取基础 Read ID
        print("提取 TRA 基础 Read ID...")
        tqdm.pandas(desc="处理 TRA Headers")
        tra_df_filtered['base_read_id'] = tra_df_filtered['descrsR1'].progress_apply(get_base_read_id)
        tra_df_filtered = tra_df_filtered.dropna(subset=['base_read_id'])

        # 为合并做准备：重命名列
        tra_df_renamed = tra_df_filtered.rename(columns={
            'bestVGene': 'TRA_VGene',
            'bestJGene': 'TRA_JGene',
            'nSeqCDR3': 'TRA_nCDR3',
            'aaSeqCDR3': 'TRA_aCDR3',
            'descrsR1': 'TRA_Header'
        })

        # 选择列，并处理重复的 base_read_id
        tra_df_unique = tra_df_renamed[['base_read_id', 'TRA_Header', 'TRA_VGene', 'TRA_JGene', 'TRA_nCDR3', 'TRA_aCDR3']].drop_duplicates(subset=['base_read_id'], keep='first').reset_index(drop=True)
        print(f"处理后得到 {len(tra_df_unique)} 条唯一的 TRA Read ID 及其比对信息。")
        del tra_df_filtered, tra_df_renamed # 释放内存

    except Exception as e:
        print(f"错误: 读取或处理 TRA 导出文件 '{tra_export_file}' 时出错: {e}", file=sys.stderr)
        sys.exit(1)

    # --- 4. 读取、过滤并处理 TRB 导出文件 ---
    try:
        print(f"读取 TRB 导出文件: {trb_export_file}")
        trb_df = pd.read_csv(trb_export_file, sep='\t', low_memory=False)
        print(f"读取了 {len(trb_df)} 条 TRB 原始比对记录。")

        # *** 新增过滤步骤 ***
        trb_df['bestVGene'] = trb_df['bestVGene'].astype(str)
        trb_df['bestJGene'] = trb_df['bestJGene'].astype(str)
        # 过滤：只保留 V 或 J 基因以 'TRB' 开头的行
        trb_df_filtered = trb_df[
            trb_df['bestVGene'].str.startswith('TRB', na=False) |
            trb_df['bestJGene'].str.startswith('TRB', na=False)
        ].copy()
        print(f"过滤后剩下 {len(trb_df_filtered)} 条 TRB 比对记录。")
        del trb_df

        # 提取 TRB 基础 Read ID
        print("提取 TRB 基础 Read ID...")
        tqdm.pandas(desc="处理 TRB Headers")
        trb_df_filtered['base_read_id'] = trb_df_filtered['descrsR1'].progress_apply(get_base_read_id)
        trb_df_filtered = trb_df_filtered.dropna(subset=['base_read_id'])

        # 重命名列
        trb_df_renamed = trb_df_filtered.rename(columns={
            'bestVGene': 'TRB_VGene',
            'bestJGene': 'TRB_JGene',
            'nSeqCDR3': 'TRB_nCDR3',
            'aaSeqCDR3': 'TRB_aCDR3',
            'descrsR1': 'TRB_Header'
        })

        # 选择列并去重
        trb_df_unique = trb_df_renamed[['base_read_id', 'TRB_Header', 'TRB_VGene', 'TRB_JGene', 'TRB_nCDR3', 'TRB_aCDR3']].drop_duplicates(subset=['base_read_id'], keep='first').reset_index(drop=True)
        print(f"处理后得到 {len(trb_df_unique)} 条唯一的 TRB Read ID 及其比对信息。")
        del trb_df_filtered, trb_df_renamed

    except Exception as e:
        print(f"错误: 读取或处理 TRB 导出文件 '{trb_export_file}' 时出错: {e}", file=sys.stderr)
        sys.exit(1)

    # --- 5. 执行合并配对 ---
    print("开始合并过滤后的数据进行配对...")
    # 合并 umi_pairs 和 过滤后的 TRA 数据
    paired_df = pd.merge(
        umi_pairs_df,
        tra_df_unique,
        left_on='TRA_Read_ID_Base',
        right_on='base_read_id',
        how='inner'
    )
    paired_df = paired_df.drop(columns=['base_read_id'])
    print(f"与过滤后的 TRA 数据合并后，有 {len(paired_df)} 条配对记录。")
    del tra_df_unique # 释放内存

    # 合并结果和 过滤后的 TRB 数据
    final_paired_df = pd.merge(
        paired_df,
        trb_df_unique,
        left_on='TRB_Read_ID_Base',
        right_on='base_read_id',
        how='inner'
    )
    final_paired_df = final_paired_df.drop(columns=['base_read_id'])
    print(f"与过滤后的 TRB 数据合并后，最终得到 {len(final_paired_df)} 条完整的 TRA-TRB 配对记录。")
    del paired_df, trb_df_unique # 释放内存

    # --- 6. 输出结果 ---
    if not final_paired_df.empty:
        output_columns = [
            'TRA_Read_ID_Base', 'TRB_Read_ID_Base',
            'TRA_VGene', 'TRA_JGene', 'TRA_nCDR3', 'TRA_aCDR3',
            'TRB_VGene', 'TRB_JGene', 'TRB_nCDR3', 'TRB_aCDR3',
            'TRA_Header', 'TRB_Header'
        ]
        final_paired_df = final_paired_df[[col for col in output_columns if col in final_paired_df.columns]]

        try:
            print(f"将过滤后的配对结果写入: {output_file}")
            final_paired_df.to_csv(output_file, sep='\t', index=False)
            print("写入完成。")
        except Exception as e:
            print(f"错误: 写入输出文件 '{output_file}' 时出错: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print("警告: 过滤后未找到任何可以配对的记录，输出文件将为空或不创建。")
        try:
            with open(output_file, 'w') as f_out:
                 header = ['TRA_Read_ID_Base', 'TRB_Read_ID_Base', 'TRA_VGene', 'TRA_JGene', 'TRA_nCDR3', 'TRA_aCDR3', 'TRB_VGene', 'TRB_JGene', 'TRB_nCDR3', 'TRB_aCDR3', 'TRA_Header', 'TRB_Header']
                 f_out.write('\t'.join(header) + '\n')
            print(f"已创建空的输出文件（带表头）: {output_file}")
        except Exception as e:
             print(f"错误: 尝试创建空输出文件 '{output_file}' 时出错: {e}", file=sys.stderr)

    print("--- 外部配对 (带过滤) 完成 ---")

# --- 主执行逻辑 ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="根据 MiXCR 导出的比对文件和 UMI 配对文件，过滤交叉比对后执行外部 TRA-TRB 配对。",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--umi-pairs", default=DEFAULT_UMI_PAIRS_FILE,
                        help="包含 TRA_Read_ID_Base 和 TRB_Read_ID_Base 列的 UMI 配对 TSV 文件路径。")
    parser.add_argument("--tra-export", default=DEFAULT_TRA_EXPORT_FILE,
                        help="包含 TRA 比对信息和 descrsR1 列的 TSV 文件路径 (由 exportAlignments 生成)。")
    parser.add_argument("--trb-export", default=DEFAULT_TRB_EXPORT_FILE,
                        help="包含 TRB 比对信息和 descrsR1 列的 TSV 文件路径 (由 exportAlignments 生成)。")
    parser.add_argument("-o", "--output", default=DEFAULT_OUTPUT_FILE,
                        help="输出最终过滤后配对结果的 TSV 文件路径。")

    args = parser.parse_args()

    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.isdir(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
            print(f"已创建输出目录: {output_dir}")
        except OSError as e:
            print(f"错误: 无法创建输出目录 '{output_dir}': {e}", file=sys.stderr)
            sys.exit(1)

    perform_pairing(args.umi_pairs, args.tra_export, args.trb_export, args.output)

