# 无滑移 Basilisk 随机中点

Basilisk周期网格上的守恒随机中点法，完整非线性更新由Newton–GMRES求解，FFT只作预条件。局部面迁移率与sqrt(2φm_face)随机通量一致，无滑移m=u³；吸引项采用直接化学势面差分。

| 文件 | 用途 |
|---|---|
| basilisk_midpoint.c | C/Basilisk数值内核 |
| basilisk_midpoint_candidate.py | ctypes接口与指纹检查点，文件名沿用实现 |
| midpoint_candidate.py | NumPy参考实现 |
| adaptive_runner.py | Brownian路径、拒步和事件细化 |
| driver.py | MPI运行端、检查点和历史输出 |
| validate.py | 算子、重启和步长筛查 |
| tcw_mode.py | 本版64² TCW+mode模型 |
| build.sh、campaign.json | 编译与工况模板 |

初厚r=hi/h0=.5/.6/.7/.8/.9各5个种子，64²、φ=.01、无滑移、无排斥、每条4 MPI进程。固定共同阈值u_min=hmin/h0=.4，不是初始厚度的40%。配置与种子见 [campaign.json](campaign.json)。

原生变量H=h/hi初值1、L=100/r、A_native=.0013159/r²、φ_native=.01/r²。共同变量u=rH、τ=r*t_native/10^8，盒长1、A=13.159、m=u³。BE使用A=4π²/3，且吸引项有限网格离散不同。

设置BASILISK与MPI CC99后执行`bash build.sh`。部分旧qcc站点可用支持的qcc生成标准C，再由当地mpicc编译；运行前仍须通过本地算子门槛。[调度说明](../tools/README.md)

validate.py检查C/Python算子、导数和单步、质量、滑移值、同路径检查点重放。每厚度2验证种子，在原生时间9.6内32个共同时间对比参考步长.009375；候选.3/.15/.075/.0375/.01875，场RMS相对差≤2%、共同hmin差≤.005。仅为短测，不是首达统计收敛验证。

生产残差≤1e-10、相对质量误差<1e-9、正有限厚度；拒步不重抽噪声、不补质量或裁剪。阈值沿同一路径细化，原生时间括区≤1e-5。6小时资源窗保存后原参数续排，数值失败不当破碎。

首页结果仅画各厚度5条实际首达时间的算术平均与一条C=5理论线。理论使用64²连续Fourier谱，无插值滤波，值和斜率连续拼接到最快增长模。图是有限样本对照，不将包络交点称为随机首达平均。
