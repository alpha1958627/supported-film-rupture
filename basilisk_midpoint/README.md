# Basilisk 无滑移随机中点版

这是一套独立于 `be/` 的实现：Basilisk 周期 multigrid 网格，随机中点时间推进，Newton–GMRES 解完整非线性更新，FFT 只作线性预条件。不是旧的冻结迁移率 FFT 时间推进，也不是 BE 时间推进。历史 Basilisk 原件单独保存在 `archives/basilisk_midpoint_original_20261025.zip`，不覆盖旧验证。

| 文件 | 用途 |
|---|---|
| `basilisk_midpoint.c` | C/Basilisk随机中点数值内核 |
| `basilisk_midpoint_candidate.py` | C内核的ctypes接口与指纹检查点；文件名沿用来源 |
| `midpoint_candidate.py` | NumPy参考实现，用于逐项对照 |
| `adaptive_runner.py` | Brownian路径、自适应拒步与事件细化 |
| `driver.py` | MPI运行端、历史输出和资源窗保存 |
| `validate.py` | 算子、重启与步长筛查 |
| `tcw_mode.py` | 本版对应的TCW+mode曲线 |
| `build.sh`、`campaign.json` | 编译与25条工况表 |

本页工具命令从仓库根目录执行。[工具导航](../tools/README.md) · [结果导航](../dashboard/README.md) · [历史归档](../archives/README.md)

## 本次工况

初厚 `r=h_i/h_0=.5,.6,.7,.8,.9`，各5个独立种子；64²网格，两个方向独立面白噪声，`φ=.01`，无排斥、无滑移；每条4 MPI进程，每进程1核。全部种子列在 [campaign.json](campaign.json)。这25条是新的测试，不代替或暂停原BE任务。

**停止判据是固定参考厚度 `u_min=h_min/h_0=.4`。** 它不等于先前图片红点采用的 `h_min/h_i=.4`。后者在共同单位中是 `.4*r`，不能混合比较破碎时间。

原生变量为 `H=h/h_i`，初值1，盒长 `L=100/r`，`A_native=.0013159/r²`，`φ_native=.01/r²`。共同变量 `u=rH`、`τ=r*t_native/10^8`；共同盒长1、`A=13.159`、`m(u)=u³`。BE使用 `A=4π²/3=13.1594725348`，两者相差约0.0036%，本次保留历史Basilisk的A以追溯来源。

无滑移修改同时进入 C 内核和 Python 参考：原 `slip_native=.2/r` 改为0；面迁移率和 `sqrt(2φ*m_face)` 随之共同改变。化学势采用 `μ=-ΔH+A/H³` 的面差分；当前 BE 对吸引项使用面链式近似，有限网格上不完全相同。不要把两者差异全部归因于时间方法或CPU速度。

## 验证与运行

`validate.py` 检查 C/Python 力、噪声及方向导数、完整一步、质量、无滑移值和移除滑移后的噪声幅度关系；还检查同一 Brownian 路径检查点重启。每厚度两个独立验证种子，用参考原生步长0.009375，对比32个共同时间、至原生时间9.6。候选 `.3,.15,.075,.0375,.01875` 由大到小筛查，最大场RMS相对差≤2%、最大共同 `u_min` 差≤.005才可选用；不是长期首达收敛证明。

验证只为避免中途截断而使用原生阈值 `.01`；**生产阈值始终是 `.4/r`**。生产完整非线性残差≤1e-10、累计相对质量误差<1e-9、正有限厚度；不裁剪、不补质量、不拒步重抽噪声。阈值细化沿同一 Brownian 路径，保存首次跨越的接受步时间括区，原生括区≤1e-5。残差、正值和质量失效会保留数值失败证据，不当物理破碎。

`driver.py` 每接受步记录hmin，每300秒保存指纹检查点。6小时资源窗提前保存并原参数重排；MPI集体决策保存/退出，避免进程进入不同通信分支。重启要求源码、共享库、Python运行环境及MPI进程数指纹一致。

在配置好 `BASILISK` 和 MPI `CC99` 后执行 `bash build.sh`。SIPF使用qcc直接编译；WUZH旧qcc转换器触发内部断言，部署器在SIPF用 `qcc -source -D_GNU_SOURCE -D_MPI=1 -O3 -shared -fPIC basilisk_midpoint.c` 生成标准C，在WUZH由当地mpicc编译相同数值内核。生成物留在本地/服务器，不随仓库上传。每个运行站点仍须重新通过C/Python算子门槛。

本地部署入口：`python tools/basilisk_midpoint/cluster.py --workspace <私有工作目录> prepare` 后 `watch`。私有目录包含 sites.json、完整工况、候选及赢家路由。SIPF先提交，排队才同时竞速XH/WUZH；检查实际运行及算子门槛后保留赢家，取消其他候选并校验归档后清理。站点配置与原始场不公开。

## 图与理论

[实时派生看板](../dashboard/basilisk_noslip/README.md) 由 `python tools/basilisk_midpoint/report.py --workspace <私有工作目录>` 刷新。实线画所有已保存种子；方形末点为未到阈值，圆形末点为实际跨越。公开仅抽样hmin、门槛验证和脱敏状态，完整历史与场留在服务器。

`tcw_mode.py` 用本实现的64²连续Fourier白噪声谱、A=13.159、无滑移和 `C_rms=5`，不加BE旧理论中的双线性插值传递函数；值及斜率连续拼接至最快增长模。每个厚度只有一条确定模型线，所有种子共用。它是线性涨落包络与模式增长预测，不是随机首达平均；与BE旧256²、sinc⁴过滤的模型线有定义差异。
