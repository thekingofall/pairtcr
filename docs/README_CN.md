# PairTCR 分析流程文档

## 第一章：引言

### 1.1 流程的目标

欢迎使用 PairTCR 分析流程。本软件的主要目标是从双端高通量测序数据中鉴定和定量配对的T细胞受体（TCR）α链（TRA）和β链（TRB）。在免疫学中，TCR是T淋巴细胞（或T细胞）表面的一个蛋白质复合体，负责识别与主要组织相容性复合体（MHC）分子结合的抗原肽段。

这种识别的特异性由TRA链和TRB链的独特组合决定。了解构成单个T细胞受体的确切TRA/TRB配对，对于理解其抗原特异性、追踪疫苗或疾病引起的克隆反应以及开发新型免疫疗法至关重要。

本流程旨在从原始测序读段（reads）中重建这些关键的TRA-TRB配对，将海量的基因组数据转化为具有实际意义的生物学见解。

### 1.2 设计哲学

本流程基于几个关键原则构建，以确保其健壮性、用户友好性和高效率：

1.  **模块化与关注点分离**：整个工作流程被分解为多个独立的、顺序执行的步骤，每个主要任务都由一个专门的脚本处理（例如 `1_preprocess_and_trim.py`、`2_create_umi_pairs.py` 等）。这种模块化设计使得流程更易于理解、调试和修改。它隔离了复杂的处理过程，例如基于UMI的读段连接和生物序列注释，从而实现了专注的开发和验证。

2.  **自动化与流程编排**：`5_runpipeline.py` 脚本充当主控制器或编排器。它能自动按正确顺序执行所有步骤，管理中间文件的数据流，处理日志记录，并为用户提供一个统一、简单的界面。这对于确保可重复性和易用性至关重要。

3.  **状态管理与可恢复性**：本流程被设计为容错的。在执行一个步骤之前，编排脚本会检查该步骤预期的输出文件是否存在。如果文件存在且不为空，则跳过该步骤。这种智能的状态检查机制使得流程能够在发生错误后从中断点恢复，从而在长时间运行中途出错时节省大量的计算时间。同时提供了一个 `--force` 标志来覆盖此行为，强制从头开始整个过程。

4.  **性能优化**：工作流程中内置了一个关键但非显而易见的优化。在基于UMI连接识别出潜在的TRA-TRB配对后（步骤2），流程会创建一个新的、经过筛选的FASTQ文件集，其中*仅包含*属于有效配对的读段（步骤2.5）。这显著减少了后续注释步骤中计算密集型工具MiXCR需要处理的读段数量，从而大幅提升性能并缩短运行时间。

5.  **清晰性与易于调试**：本流程为透明度而设计。所有中间文件都存储在与每个步骤相对应的、命名清晰的专用目录中。此外，还为整个流程和每个独立步骤生成了全面的日志。这种结构化的输出使得追踪数据流、检查中间结果以及在分析的任何阶段诊断问题都变得简单明了。

## 第二章：流程工作流与数据流

本流程执行一系列脚本，将数据从原始读段处理成最终的配对克隆。下图展示了整体的工作流程。

```mermaid
graph TD
    subgraph "原始数据"
        A[R1/R2 FASTQ 文件]
    end

    subgraph "步骤 1: 预处理"
        B(1_preprocess_and_trim.py)
        A --> B
        B --> C[TRA 特异性 FASTQ 文件<br/>(头部含 UMI)]
        B --> D[TRB 特异性 FASTQ 文件<br/>(头部含 UMI)]
    end

    subgraph "步骤 2: UMI 配对"
        E(2_create_umi_pairs.py)
        C --> E
        D --> E
        E --> F[umi_pairs.tsv<br/>(TRA_ReadID <-> TRB_ReadID)]
    end

    subgraph "步骤 2.5: FASTQ 筛选 (优化)"
        G(主流程脚本筛选FASTQ)
        F --> G
        C --> G
        D --> G
        G --> H[匹配的 TRA FASTQ<br/>(仅含配对中的读段)]
        G --> I[匹配的 TRB FASTQ<br/>(仅含配对中的读段)]
    end

    subgraph "步骤 3: TCR 注释"
        J(3_run_mixcr_and_export.sh)
        H --> J
        I --> J
        J --> K[TRA 比对结果 TSV<br/>(MiXCR 输出)]
        J --> L[TRB 比对结果 TSV<br/>(MiXCR 输出)]
    end

    subgraph "步骤 4: 最终组装与筛选"
        M(4_pair_and_filter_clones.py)
        F --> M
        K --> M
        L --> M
        M --> N{最终配对克隆表<br/>(TRA_CDR3, TRB_CDR3)}
    end

    style A fill:#e6e6e6,stroke:#333,stroke-width:2px
    style N fill:#b4e7b4,stroke:#333,stroke-width:2px
```

### 2.1 步骤详解

#### **步骤 1: 预处理与剪切 (`1_preprocess_and_trim.py`)**

*   **目的**: 识别源自TRA或TRB转录本的读段，提取其UMI，并清理序列。
*   **输入**: 原始的双端FASTQ文件 (`*_1.fq.gz`, `*_2.fq.gz`)。
*   **过程**:
    1.  脚本扫描每个读段对，寻找对TRA和TRB建库特有的特定DNA接头和连接子序列。
    2.  如果找到特定结构，它会从中提取两个唯一分子标识符（UMI）。
    3.  UMI信息被连接起来并附加到读段对中两条读段的FASTQ头部（例如 `@READ_ID UMI:TRA:GATTACA_CATTGTC`）。
    4.  读段序列被剪切，以去除所有非生物学碱基（接头、连接子、UMI），仅留下内源性的TCR转录本序列用于下游分析。
*   **输出**: 四个新的、按链类型分开的gzipped FASTQ文件：`PREFIX_TRA_1.fq.gz`, `PREFIX_TRA_2.fq.gz`, `PREFIX_TRB_1.fq.gz`, 和 `PREFIX_TRB_2.fq.gz`。

#### **步骤 2: 创建UMI配对 (`2_create_umi_pairs.py`)**

*   **目的**: 将源自同一个单细胞的TRA和TRB读段连接起来。
*   **过程**: 此步骤利用了文库设计的关键特性：TRB读段的UMI序列是其对应TRA读段UMI的反向互补序列。
    1.  首先，它读取所有 `_TRA_1.fq.gz` 记录，并构建一个字典，将每个唯一的TRA UMI映射到携带它的读段标识符集合。
    2.  然后，它计算每个TRA UMI的反向互补序列，并创建第二个映射：`reverse_complement(TRA_UMI) -> original_TRA_UMI`。
    3.  最后，它遍历 `_TRB_1.fq.gz` 文件。对于每个TRB读段，它提取UMI并检查此UMI是否存在于反向互补映射中。
    4.  如果找到匹配，意味着已识别出一个配对。脚本记录下TRB读段的基础Read ID与相应TRA读段之间的关联。
*   **输出**: 一个制表符分隔的文件 (`umi_pairs.tsv`)，包含两列：`TRA_Read_ID_Base` 和 `TRB_Read_ID_Base`，明确列出了源自单个原始RNA分子的每个读段配对。

#### **步骤 2.5: 创建匹配的FASTQ (在 `5_runpipeline.py` 内部执行)**

*   **目的**: 通过减小计算最密集步骤（MiXCR）的输入大小来优化流程。
*   **过程**:
    1.  流程编排器读取 `umi_pairs.tsv` 文件，并编译一个包含所有成功配对的TRA和TRB读段ID的唯一集合。
    2.  然后，它筛选来自步骤1的四个FASTQ文件，写入新的FASTQ文件，其中*仅包含*ID在此"已配对"集合中的读段。
*   **输出**: 一组新的、更小的、"匹配的"FASTQ文件，位于 `matched_fastq_output/` 目录中。

#### **步骤 3: 运行MiXCR并导出 (`3_run_mixcr_and_export.sh`)**

*   **目的**: 使用MiXCR工具进行TCR序列注释。
*   **输入**: 来自步骤2.5的、经过筛选的"匹配的"FASTQ文件。
*   **过程**: 这个shell脚本自动化了MiXCR的工作流程。
    1.  它分别对TRA和TRB的FASTQ配对文件运行 `mixcr analyze amplicon`。此步骤将读段与V、D、J基因的参考库进行比对。至关重要的是，它使用了 `--align "-OsaveOriginalReads=true"` 选项，以确保保留原始的FASTQ头部（带有我们的UMI标签）。
    2.  然后，它运行 `mixcr assemble` 从比对结果中构建克隆。
    3.  最后，它使用 `mixcr exportAlignments` 生成详细的、人类可读的报告。导出格式被配置为包括V/J基因的鉴定结果、CDR3区域的核苷酸和氨基酸序列，以及原始读段的描述信息（`-descrsR1`），其中包含了最终配对所需的、至关重要的读段ID。
*   **输出**: 两个TSV文件：`TRA_alignments_export_with_headers.tsv` 和 `TRB_alignments_export_with_headers.tsv`。

#### **步骤 4: 配对与筛选克隆 (`4_pair_and_filter_clones.py`)**

*   **目的**: 将所有处理过的数据整合到一个最终的、全面的配对TRA-TRB克隆表中。
*   **输入**: 来自步骤2的 `umi_pairs.tsv` 文件和来自步骤3的两个MiXCR比对报告。
*   **过程**: 这个Python脚本使用 `pandas` 库执行一系列数据库风格的连接操作。
    1.  首先加载两个MiXCR比对表，并进行筛选，以移除跨链的错误比对（例如，来自TRB文库文件的读段被MiXCR错误地注释为TRA基因）。
    2.  它从每个表的头部信息列中提取基础读段ID。
    3.  它将 `umi_pairs.tsv` 表与筛选后的TRA比对表合并，连接键为 `TRA_Read_ID_Base`。
    4.  然后，它将合并后的结果与筛选后的TRB比对表再次合并，连接键为 `TRB_Read_ID_Base`。
*   **输出**: 整个流程的最终结果：`final_paired_clones_filtered.tsv`。该文件中的每一行代表一个成功测序、组装和注释的TRA-TRB配对，包含了两条链的V/J基因使用情况和CDR3序列。

## 第三章：安装与使用

### 3.1 依赖项

- Python 3.x
- `pandas` 库
- `tqdm` 库
- MiXCR (v3.x 或更高版本)

### 3.2 运行流程

整个流程通过 `5_runpipeline.py` 脚本执行。

**基本用法:**

```bash
python3 scripts/5_runpipeline.py /path/to/raw/fastq/dir
```

**常用参数:**

*   `input_dir`: (位置参数) 包含原始 `_1.fq.gz` 和 `_2.fq.gz` 文件的目录路径。默认为 `./raw`。
*   `-o, --output-root`: 用于存储所有结果的根目录。默认为 `./PairTCR_results`。
*   `-p, --prefix`: 所有生成文件的统一前缀。默认为 `TCR_TSO_18`。
*   `-n, --read-limit`: 从输入FASTQ中处理的最大读段对数量。便于测试。默认为 100,000。
*   `-t, --threads`: MiXCR使用的线程数。默认为 90。
*   `--mixcr-jar`: 你的 `mixcr.jar` 文件的路径。默认为 `scripts/mixcr.jar`。
*   `--force`: 强制流程从头开始，删除所有先前结果。
*   `--use-c`: 使用预编译的C语言版预处理器脚本，可显著加快步骤1的速度。

**示例:**

```bash
python3 scripts/5_runpipeline.py MyRawData \
    --output-root MyProject_Results \
    --prefix PatientA \
    --read-limit 5000000 \
    --threads 16 \
    --mixcr-jar /path/to/my/mixcr.jar
```
此命令将:
- 处理来自 `MyRawData` 目录的前500万个读段对。
- 将所有结果存储在一个名为 `MyProject_Results` 的顶级文件夹中。
- 用前缀 `PatientA` 命名所有中间文件。
- 在MiXCR分析步骤中使用16个线程。
- 使用位于 `/path/to/my/mixcr.jar` 的MiXCR可执行文件。

## 第四章：输出结构

成功完成后，输出根目录将包含以下子目录：

*   `logs/`: 包含主流程和每个独立步骤的详细日志文件。对调试至关重要。
*   `1_preprocess_and_trim_output/`: 包含来自步骤1的四个中间FASTQ文件。
*   `2_create_umi_pairs_output/`: 包含来自步骤2的 `umi_pairs.tsv` 文件。
*   `matched_fastq_output/`: 包含来自步骤2.5优化的、经过筛选的FASTQ文件。
*   `3_run_mixcr_and_export_output/`: 包含MiXCR报告、比对文件 (`.vdjca`)、克隆文件 (`.clns`) 和来自步骤3的导出的TSV比对文件。
*   `4_pair_and_filter_clones_output/`: 包含流程的最终主要输出文件：`final_paired_clones_filtered.tsv`。

---
最终的 `final_paired_clones_filtered.tsv` 文件是主要结果，提供了样本中所有已识别的TRA-TRB配对的全面列表。 