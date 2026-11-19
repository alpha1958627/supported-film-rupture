# BE–MG：当前无滑移后向欧拉版

迁移率m=h³，保守乘性面噪声；确定性通量和噪声在末态计算，以非线性外迭代和Basilisk多重网格求解完整更新。空间吸引项采用面链式近似。

| 文件 | 用途 |
|---|---|
| wall_native.c | C/Basilisk数值内核 |
| noise.py、colored.py | 白噪声与相关噪声实现；本批次使用white6 |
| worker.py、launch.py | 运行端自适应推进、检查点与启动 |
| validate.py | 提交生产前的运行端验证 |
| tcw_mode.py | 原登记TCW+mode模型，不随种子拟合 |
| build.sh | Linux/Basilisk编译入口 |

本批次初厚0.5–0.9各5条，φ=.01、无滑移、64²、每条4物理OpenMP核；实际停止阈值是固定hmin/h0=.005。0.3与0.4图另外提取统计事件，不更改运行端目标。

[工况与种子](../campaign.json) · [验证回执](../validation/README.md) · [结果](../dashboard/README.md) · [方程、验证和编译步骤](../docs/be_method_and_results.md) · [调度与作图](../tools/README.md)
