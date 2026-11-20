# 无滑移 BE–MG

局部迁移率m=h³，末态保守乘性面噪声；随机后向欧拉，以非线性外迭代和Basilisk多重网格求解。共同单位盒长1、A=4π²/3，吸引项使用面厚度链式近似。

| 文件 | 用途 |
|---|---|
| wall_native.c | C/Basilisk内核 |
| noise.py、colored.py | 白噪声与相关噪声实现 |
| worker.py、launch.py | 自适应推进、检查点和启动 |
| validate.py | 算子、离散模态、步长和重启验证 |
| tcw_mode.py | 此路线的原登记理论模型 |
| build.sh | Linux编译 |

初厚.5–.9各5条的配置见 [campaign.json](../campaign.json)，φ=.01、64²、每条4物理OpenMP核，运行端停止于固定hmin/h0=.005。首页图使用另一套中点版的.4事件，不能混用阈值或事件时间。

编译：设置BASILISK到Basilisk src目录，在本目录运行`bash build.sh`。将共享库复制到私有工作目录，填写sites.json，按 [工具说明](../tools/README.md) 运行prepare与watch。完整更新残差<1e-8、相对质量误差<1e-9、正有限厚度；不裁剪、不补质量。拒步沿同一Brownian路径细分。

验证通过后才允许提交相应厚度的生产轨迹。短窗步长对照不等于长期首达统计收敛。实际结果、验证回执、场与检查点只存私有目录。
