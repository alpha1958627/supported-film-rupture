# 工具入口

从仓库根目录执行命令。工具本身不是新的数值版本；核心分别在be/和basilisk_midpoint/。

| 方法 | 文件（相对tools/） | 职责 |
|---|---|---|
| BE | be/cluster.py | SIPF优先、排队竞速、验证后提交、重复作业归档清理 |
| BE | be/collect.py | 核验作业和已保存历史，采集私有快照 |
| BE | be/render.py | 全部种子、TCW+mode和HTML看板 |
| 中点 | basilisk_midpoint/cluster.py | MPI编译、验证门槛和调度 |
| 中点 | basilisk_midpoint/report.py | 核验SHA，采集并生成中点版看板 |
| 中点 | basilisk_midpoint/mean.py | 完整历史按共同物理时间对齐，计算每厚度5条停止路径的算术平均 |
| 统计 | analysis/compare_events.py | 初始厚度1%事件与旧Basilisk对照 |
| 统计 | analysis/compare_h03.py | 固定参考厚度0.3或0.4首次事件 |
| 原图 | analysis/review_c5_image_reference.py | C=5及原图数据；初始厚度40%判据 |

只更新结果：

```bash
python tools/be/collect.py --workspace <BE-workspace>
python tools/be/render.py --workspace <BE-workspace>
python tools/basilisk_midpoint/report.py --workspace <midpoint-workspace>
```

统计比较：

```bash
python tools/analysis/compare_events.py --workspace <BE-workspace>
python tools/analysis/compare_h03.py --workspace <BE-workspace>
python tools/analysis/compare_h03.py --workspace <BE-workspace> --threshold 0.4
python tools/analysis/review_c5_image_reference.py --workspace <BE-workspace> --source-root <historical-data-directory>
```

`cluster.py ... prepare`只打包本地运行文件；`watch`会实际提交和管理任务。两套方法使用不同私有工作目录。真实sites.json不上传；[BE站点模板](../config/README.md)。
