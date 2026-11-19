# 历史版本与恢复

这些ZIP是冻结原件。当前入口：[BE–MG](../be/README.md)、[无滑移随机中点](../basilisk_midpoint/README.md)。

| 归档 | 内容 | 对应方法 |
|---|---|---|
| [legacy_code_20261119.zip](legacy_code_20261119.zip) | 原仓库源码、原README/Notebook及既有补充代码 | JAM、冻结迁移率C–FFT |
| [legacy_validation_20261119.zip](legacy_validation_20261119.zip) | 上述旧方法的验证、图像和动画 | 只验证对应旧方法 |
| [be_pilot_and_checks_20261119.zip](be_pilot_and_checks_20261119.zip) | 0.7 BE先导试跑及换单位检查 | 当前BE路线的先导记录 |
| [basilisk_midpoint_original_20261025.zip](basilisk_midpoint_original_20261025.zip) | 原中点源码、Python参考、Slurm及源清单 | 含滑移0.2；旧“.3”指原生最大步长 |

恢复原仓库：将legacy_code和legacy_validation两个ZIP解压到一个新的空目录。其他两个ZIP各自单独解压，避免同名文件互相覆盖。

校验：[旧仓库逐文件与ZIP SHA](manifest.json) · [BE先导ZIP SHA](extra_checksums.json) · [原中点逐文件与ZIP SHA](../basilisk_midpoint/campaign.json)。目录整理不重打包或修改原件。
