# 附壁液膜

本仓库用二维长波模型推进刚性固壁上液膜的三维自由表面 $`z=h(x,y,t)`$。壁面无滑移，气相剪切忽略。表面张力抹平起伏，范德华吸引可使长波扰动增长，热噪声则驱动稳定平膜的厚度涨落。算例覆盖近乎平膜阶段。

## 模型与时间推进

记无量纲厚度为 $`H=h/h_0`$，平膜对应 $`H=1`$。长度以 $`h_0`$ 归一化，时间以 $`\mu h_0/\gamma`$ 归一化。迁移率在平膜上取值，压力保留非线性的范德华项。无热噪声时，方程为

```math
\partial_{t^*}H=-\mathcal M P(H),\qquad
P(H)=-\nabla_*^2 H+a(H^{-3}-1).
```

其中 $`a`$ 为无量纲吸引系数，$`\nabla_*`$ 对无量纲膜面坐标求导。对模态 $`(m_x,m_y)`$，无量纲波数向量为 $`\boldsymbol q=2\pi(m_x/L_x^*,m_y/L_y^*)`$，其中 $`L_i^*=L_i/h_0`$，并记 $`q=|\boldsymbol q|`$。无滑移固壁与零切应力自由表面给出迁移率算子的 Fourier 乘子

```math
M(q)=\frac{q^2}{3},\quad q>0,\qquad M(0)=0.
```

代码取无量纲液体黏度和表面张力为 1。将平膜线性化，连续理论给出模态率和振幅

```math
s^*(q)=\frac{q^2}{3}\left(3a-q^2\right),\qquad
A(t^*)=A(0)\,e^{s^*(q)t^*}.
```

纯表面张力下，所有非零模态衰减。有吸引时，$`0<q<\sqrt{3a}`$ 的模态增长。正则模对照通过拟合 $`\log[A(t^*)/A(0)]`$ 的斜率得到数值模态率。

两套程序在周期等距网格上使用中心差分拉普拉斯算子，在 Fourier 空间作用迁移率，并按隐式中点推进

```math
H^{n+1}-H^n=-\Delta t^*\mathcal M_d P_d\!\left(\tfrac{H^{n+1}+H^n}{2}\right)
 + \frac{\sqrt{2\theta\Delta t^*\mathcal M_d}}{\Delta x^*}\,\xi^n.
```

$`P_d`$ 和 $`\mathcal M_d`$ 分别为离散压力和迁移率。设 $`q_{d,i}=2\sin(q_i\Delta x^*/2)/\Delta x^*`$，离散拉普拉斯的乘子为 $`-q_d^2=-(q_{d,x}^2+q_{d,y}^2)`$，离散迁移率为 $`M_d=q_d^2/3`$。平方根按 Fourier 乘子逐模态计算。$`\theta=k_BT/(\gamma h_0^2)`$ 为无量纲噪声强度，$`\xi^n`$ 是该步的标准高斯网格场。持续加噪每步输入新的噪声场，仅首步加噪从第二步起将噪声项置零。

模型与 TCW 理论参考 [Zhang 等，Phys. Rev. E 102, 053105](https://doi.org/10.1103/PhysRevE.102.053105) 和 [Zhang 等，J. Fluid Mech. 915, A135](https://doi.org/10.1017/jfm.2021.164)。

## 二维无涨落正则模

两组均采用 $`64\times64`$ 周期网格，取 $`L_x/h_0=L_y/h_0=128`$、$`\Delta t^*=8`$。初始扰动为斜向 $`(m_x,m_y)=(2,2)`$ 的余弦波，振幅 $`10^{-3}h_0`$，热噪声为零。

| 工况 | JAM 与 Basilisk 拟合的模态率 | 连续长波理论 |
| --- | --- | --- |
| 有吸引，$`a=0.0084375`$ | $`3.90562\times10^{-5}`$ | $`3.87840\times10^{-5}`$ |
| 无吸引，$`a=0`$ | $`-1.23068\times10^{-4}`$ | $`-1.23862\times10^{-4}`$ |

左图只展示有吸引时 Basilisk 在 $`t^*=640`$ 的厚度场。右图在同一坐标框内比较增长和衰减，黑线为连续理论。图的横纵轴、线性尺度、色标和范围沿用 [自由液膜仓库](https://github.com/alpha1958627/free-film-rupture)。附壁膜的演化较慢，到 $`t^*=640`$，两组振幅分别为初值的 102.5% 和 92.4%。

![有吸引的末态厚度场，以及有无吸引的模态振幅对照](postprocess/figures/normal_mode_2d.png)

[正则模.ipynb](postprocess/正则模.ipynb) 用保存的厚度场和振幅历史重画合图。[合图 PDF](postprocess/figures/normal_mode_2d.pdf)

## 二维热涨落厚度场

[tcw.ipynb](postprocess/tcw.ipynb) 比较持续加噪与首步加噪的二维壳平均谱。对无吸引的稳定平膜，记非零 Fourier 模态的功率为 $`S(q,t^*)=\langle|\widehat{H-1}|^2\rangle`$。在边长为 $`L^*`$ 的正方形周期域，连续理论给出

```math
\begin{aligned}
S_{\rm on}(q,t^*)&=\frac{\theta}{(L^*)^2q^2}\left[1-e^{2s^*(q)t^*}\right],\\
S_{\rm off}(q,t^*)&=S(q,t_0^*)\,e^{2s^*(q)(t^*-t_0^*)},\qquad t_0^*=\Delta t^*.
\end{aligned}
```

持续加噪从平膜起步，首步加噪以实测的首步功率谱为初值。图中对每个波数壳的功率取平均后开平方。横轴为 $`|\boldsymbol q|h_0`$，纵轴为 $`10^3|\widehat{\delta h}/h_0|_{\rm rms}`$，均采用与自由液膜仓库相同的线性尺度和范围。

下图来自 SIPF 作业 `62083`，使用 $`512\times512`$ 网格。持续加噪和仅首步加噪各计算 30 对 JAM/Basilisk 轨迹，保留 84 个完整波数壳，并与连续长波理论同图对照。

![SIPF 512×512 网格的 TCW 谱，比较 JAM、Basilisk 与连续理论](postprocess/figures/tcw_sipf_512.png)

[SIPF 512×512 TCW 谱 PDF](postprocess/figures/tcw_sipf_512.pdf)

动图取自 $512\times512$ 周期网格上种子 `61100` 的配对轨迹，左侧为 JAM，右侧为 Basilisk。两组显示相同的 $256h_0\times256h_0$ 区域，使用与自由液膜仓库相同的 82 个时刻、视角、纵轴范围和色标。

### 持续加噪

![持续加噪下 JAM 与 Basilisk 的三维厚度场并排动画](postprocess/figures/tcw_continuous_3d_paired.gif)

### 仅首步加噪

![仅首步加噪下 JAM 与 Basilisk 的三维厚度场并排动画](postprocess/figures/tcw_first_step_3d_paired.gif)

## 附壁液膜破碎与预测对照

![附壁液膜初始厚度与破碎时间、最低厚度历史及原预测的对照](postprocess/figures/wall_film_rupture_prediction_20261112.png)

## 代码入口

- [JAM 二维求解器](jam/wall_film_2d.py) 的 `WallFilm2D.advance(h, white)` 推进一步。线性项在 Fourier 空间消元，有吸引时用解析雅可比修正非线性残差。
- [Basilisk 二维求解器](basilisk/wall_film_midpoint_2d.c) 用 `qcc` 编译，调用 [二维 FFT](basilisk/gas_fft_2d.h)，用定点迭代求解非线性中点压力。程序从标准输入读取初始厚度和逐步高斯场，从标准输出返回保存的厚度帧。
- [配对运行与谱统计](postprocess/run_tcw.py) 给两个求解器传入相同的初态和随机场，再计算 Fourier 壳平均谱。[后处理目录](postprocess/) 包含正则模与 TCW 的 Notebook、绘图脚本和 SIPF 的逐种子壳功率。图轴参考值保存在 [figure_reference.json](postprocess/data/figure_reference.json)。

## 复算

需要 NumPy、Matplotlib、Jupyter Python 内核、Pillow 和用于图中文字渲染的 LaTeX。使用 Basilisk `qcc` 编译二维中点程序。

```sh
python -m pip install -r requirements.txt
mkdir -p build
qcc -O3 -disable-dimensions -I"$PWD/basilisk" basilisk/wall_film_midpoint_2d.c -o build/wall_film_midpoint_2d -lm
python postprocess/run_tcw.py build/wall_film_midpoint_2d runs/tcw_512 --grid 512 --replicas 30 --workers 6
python postprocess/figures.py
python postprocess/replay_surfaces.py build/wall_film_midpoint_2d
python postprocess/render_surfaces.py
```

$512^2$ 算例使用 $L_x/h_0=L_y/h_0=1024$，保留 $\Delta x/h_0=\Delta y/h_0=2$。在相同波数上限 $|\boldsymbol q|h_0\leq\pi/6$ 内，共有 84 个完整的非零 Fourier 壳。终点 $t^*=1280$，每条轨迹 160 步。两种协议各有 30 对轨迹，每对包含 JAM 和 Basilisk。`run_tcw.py` 把结果写入指定目录，Notebook 和绘图脚本默认读取 `postprocess/data` 中发布的数据。

SIPF `node2` 于北京时间 2026-11-11 02:28:10–02:33:54 完成原批次，实测 5 分 44 秒。持续加噪与首步加噪谱相对离散中点理论的壳 RMS 偏差中位数分别为 0.56% 和 0.55%。首步加噪在长时间高波数端的离散衰减与连续理论可见差异，CSV 同时保留两种理论。原批次 $128^2$ 网格到 $t^*=16384$ 的发展率检查保存在 [长时间验证数据](postprocess/data/validation_long_window.json)。
