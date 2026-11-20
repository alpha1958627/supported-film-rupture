# 运行与绘图工具

从仓库根目录执行。两套方法各用独立私有工作目录，真实sites.json不上传。

| 文件 | 用途 |
|---|---|
| be/cluster.py | BE的SIPF优先调度、验证门槛和候选清理 |
| basilisk_midpoint/cluster.py | 中点版的MPI编译、调度和验证门槛 |
| plot_rupture.py | 读取私有核验快照，只导出首页一张结果图 |

BE：先在be/编译，将wall_native.so和站点配置放到私有workspace。

```bash
python tools/be/cluster.py --workspace <private-BE-workspace> prepare
python tools/be/cluster.py --workspace <private-BE-workspace> watch
```

中点版：在私有workspace配置各站MPI、Python和Basilisk路径，并复制中点campaign.json。

```bash
python tools/basilisk_midpoint/cluster.py --workspace <private-midpoint-workspace> prepare
python tools/basilisk_midpoint/cluster.py --workspace <private-midpoint-workspace> watch
```

watch会实际提交和管理任务；prepare只打包本地运行文件。SIPF先投，排队时竞速XH与WUZH，通过运行门槛后保留赢家，其他副本校验归档后清理。只处理私有回执登记的本批次作业。

从含25条已完成中点轨迹、参数、首达括区和来源SHA的私有核验快照重画：

```bash
python tools/plot_rupture.py --snapshot <private-path>/bm_snapshot.json
```

脚本只向results/写一张图的PDF和600 dpi PNG；逐条事件与均值核验表写在输入快照旁边，不在仓库公开。不提交快照、采样曲线、CSV、验证数据或旧看板。
