# 历史对照数据

这里是历史数据，不是求解器或新无滑移批次。

- [basilisk_cap03.csv](basilisk_cap03.csv)：2026-10-25原生最大步长0.3的历史摘要；滑移0.2，原16 MPI核。
- [provenance.json](provenance.json)：定义、来源与SHA。
- [原图C=5核对](../dashboard/review_C5_image_reference/README.md)：包含后来补充的0.5组，覆盖范围与上述原批次不同。

原生H=h/hi；共同u=(hi/h0)H，τ=(hi/h0)t_native/10^8。原图Hmin=.4对应u_min=.4*(hi/h0)，不是固定u_min=.4。稀疏保存上界、插值估计与新轨迹接受步首达括区按各图说明区分。

[历史源码](../archives/basilisk_midpoint_original_20261025.zip) · [当前无滑移中点版](../basilisk_midpoint/README.md)
