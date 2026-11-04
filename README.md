# DNABarcodeDiscovery

## usage

```bash
# 第一步：计算窗口 Pi，生成 window_pi.tsv 和可交互 HTML 点图
uv run barcode-discovery scan aligned.fasta -o results --pi-threshold 0.05

# 生成配置模板
uv run barcode-discovery init-config -o config.toml

# 在点图上悬停读取窗口起点，填好 config.toml 后运行第二步
uv run barcode-discovery evaluate config.toml
```
