# 附壁液膜破碎：代码与结果导航

**目前使用两套独立求解器：BE–MG 与 Basilisk 随机中点。** 两者均计算无滑移固壁上的二维高度场，对应三维液膜自由表面。先选方法，再进入该方法的源码、验证与看板。

## 版本一览

| 版本 | 状态与用途 | 方法 | 源码 / 验证 / 结果 |
|---|---|---|---|
| **BE–MG，无滑移** | 当前25条轨迹；初厚0.5–0.9各5条 | 局部h³；随机后向欧拉；多重网格与非线性迭代 | [源码](be/README.md) · [验证](validation/README.md) · [看板](dashboard/index.html) |
| **Basilisk 随机中点，无滑移** | 新增25条测试；初厚0.5–0.9各5条 | 局部h³；随机中点；Newton–GMRES，FFT预条件 | [源码](basilisk_midpoint/README.md) · [验证](dashboard/basilisk_noslip/validation.json) · [看板](dashboard/basilisk_noslip/README.md) |
| Basilisk 随机中点，滑移0.2 | 历史对照；旧图“Basilisk .3” | h³+0.2h²；随机中点；原16 MPI核 | [代码ZIP](archives/basilisk_midpoint_original_20261025.zip) · [数据来源](reference/README.md) |
| BE 0.7先导试跑 | 归档；换单位及试跑检查 | 当前BE路线的早期独立试跑 | [代码与检查ZIP](archives/be_pilot_and_checks_20261119.zip) |
| 原仓库 JAM / C–FFT | 冻结归档 | 冻结迁移率、加性噪声；Fourier中点推进 | [代码ZIP](archives/legacy_code_20261119.zip) · [旧验证ZIP](archives/legacy_validation_20261119.zip) |

**“Basilisk .3”中的 `.3` 是原生最大时间步，不是版本号或破碎阈值。** BE也使用Basilisk基础设施；两套名称区分时间推进与非线性求解方法。

## 两套当前工况的区别

| 设置 | BE–MG | Basilisk 随机中点 |
|---|---|---|
| 网格 / 噪声 | 64² / white6 | 64² / 面白噪声 |
| 强度 / 滑移 | φ=0.01 / 无滑移 | φ=0.01 / 无滑移 |
| 共同单位中的吸引系数 | A=4π²/3 | A=13.159，保留历史系数 |
| 吸引项离散 | 面厚度上的链式近似 | 直接化学势面差分 |
| 每条资源 | 4物理OpenMP核 | 4 MPI进程，每进程1核 |
| **实际停止阈值** | **hmin/h0=0.005** | **hmin/h0=0.4** |
| 工况与种子 | [campaign.json](campaign.json) | [campaign.json](basilisk_midpoint/campaign.json) |
| TCW+mode图线 | 旧主图C=3.5；另有C=5版；256²谱和插值传递 | C=5；64²连续谱，无插值滤波 |

`h0`是固定参考厚度，`hi`是每组初始厚度。**固定hmin/h0=0.4与hmin/hi=0.4是两个判据。** [原图C=5核对版](dashboard/review_C5_image_reference/README.md) 使用后者；各图分别标明。

## 目录分布

```text
be/                    当前BE–MG内核、噪声、理论、运行端验证
basilisk_midpoint/     当前无滑移中点内核、参考实现、理论、运行端验证
tools/
  be/                  BE调度、采集、看板生成
  basilisk_midpoint/    中点版调度、采集和看板生成
  analysis/            事件、阈值、理论和历史数据比较
config/                站点配置模板
validation/            BE实际短测回执及失败版本记录
reference/             历史Basilisk对照数据与来源
dashboard/             已发布结果；basilisk_noslip/为新中点版看板
postprocess/           两版共享绘图样式，不是第三套求解器
archives/              原件ZIP、恢复说明、SHA校验
docs/                  方法和结果的详细记录
campaign.json          BE的25条轨迹工况表
```

## 常用入口

从仓库根目录运行；`<workspace>`是含站点配置、作业回执和私有检查点的本地工作目录。

| 操作 | BE–MG | Basilisk 随机中点 |
|---|---|---|
| 调度 | `python tools/be/cluster.py --workspace <workspace> watch` | `python tools/basilisk_midpoint/cluster.py --workspace <workspace> watch` |
| 采集 | `python tools/be/collect.py --workspace <workspace>` | 与下一步合并 |
| 看板 | `python tools/be/render.py --workspace <workspace>` | `python tools/basilisk_midpoint/report.py --workspace <workspace>` |

部署前的编译和prepare见各方法说明。**更新看板只需采集和作图，不必再次部署。** [工具导航](tools/README.md) · [旧入口迁移表](docs/layout.md)

[BE方法与结果详记](docs/be_method_and_results.md) · [中点版方法](basilisk_midpoint/README.md) · [归档与恢复](archives/README.md) · [结果导航](dashboard/README.md)

数值核心与历史ZIP保持原字节；服务器任务与检查点不移动。原始场、完整历史和真实账号配置不放入公开仓库。
