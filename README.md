# DNABarcodeDiscovery

通过组合 Pi 和 sequence clusters的方法实现药用植物鉴别 DNA 条形码方法。


### 依赖
本项目依赖numpy, plotly 和 pymafft包。pymafft包基于Rust重新编写了 MAFFT 算法，因此需要提前安装 Rust 工具链才能完成编译。同时要求你的python版本不低于3.10 且不大于 3.14。

### 使用方法

总共分四步骤

第一步首先执行多序列比对

```bash
uv run barcode-discovery align raw.fasta -o aligned.fasta --strategy auto -t 10
```

第二步首先计算序列多样性，会生成一个交互式点图，辅助选择最优区间用于合并：

```bash
# 第一步：计算窗口 Pi，生成 window_pi.tsv 和可交互 HTML 点图
uv run barcode-discovery scan aligned.fasta -o results --pi-threshold 0.05
```

第三步生成配置文件模板

```bash
# 生成配置模板
uv run barcode-discovery init-config -o config.toml
```

第四步填好 config.toml， 主要将需要的区间发送给下一步脚本。

```bash
# 在点图上悬停读取窗口起点，填好 config.toml 后运行第二步
uv run barcode-discovery evaluate config.toml
```

### 配置文件格式

```toml
# barcode-discovery 第二步配置文件
# 运行: barcode-discovery evaluate config.toml

# 第一步使用的 MAFFT 比对 FASTA
fasta = "data/salvia_cp_26spe_align.fas"

# 第一步输出的窗口 Pi 表
window_pi = "results/window_pi.tsv"

# 用于注释真实坐标的物种 ID (与 FASTA 头 '>' 后的 ID 完全一致, 或其第一个空白前的部分)
reference_id = "NC_020431.1"

# (可选) 参考基因组的 GFF3 注释; 必须与 reference_id 对应同一个基因组。
# 填写后 evaluate 会把每个区间注释为 geneA(10bp)+geneA-geneB(180bp)+geneB(20bp)。
# 不需要注释时留空。
gff3 = "data/S.miltiorrhiza_chloroplast.gff3"
# (可选) GFF3 含多条序列时, 指定与参考基因组对应的 seqid
gff3_seqid = ""

# 输出目录
out = "results"

# 序列聚类比阈值 (%), 大于该值判定为 PASS
ratio_threshold = 80.0

# ---- 候选条形码 ----
# 每个 [[candidate]] 是一个候选条形码, 可以写任意多个。
# intervals: 一个或多个区间, 每个区间 = [起始窗口的起点坐标, 结束窗口的起点坐标]
#   起点坐标取自 window_pi.tsv 的 start 列 (或交互图悬停显示的 window start)。
#   区间内的窗口会被融合为 [起始窗口 start, 结束窗口 end]。
#   多个区间会拼接为同一个条形码 (重叠/相邻的区间自动合并)。
[[candidate]]
name = "candidate_1"
intervals = [
    #[84, 246],
]

#[[candidate]]
#name = "candidate_2"
#intervals = [
    #[111976, 112056],
#]

```

### 结果文件

- `window_pi.tsv` : 窗口 Pi 值表
- `pi_sliding_window.html` : 窗口 Pi 值交互图
- `barcode_candidates.tsv` : 候选条形码表
