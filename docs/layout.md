# 目录整理与旧入口迁移

此次整理工具与文档，不改变数值算法、种子、接受门槛或服务器部署目录。数值内核和历史ZIP逐字节核验不变。

| 旧根目录入口 | 新入口（相对仓库根目录） |
|---|---|
| cluster.py | tools/be/cluster.py |
| collect.py | tools/be/collect.py |
| render.py | tools/be/render.py |
| bm_cluster.py | tools/basilisk_midpoint/cluster.py |
| bm_report.py | tools/basilisk_midpoint/report.py |
| compare_events.py | tools/analysis/compare_events.py |
| compare_h03.py | tools/analysis/compare_h03.py |
| review_c5_image_reference.py | tools/analysis/review_c5_image_reference.py |
| sites.example.json | config/be_sites.example.json |

参数不变，只更换脚本路径。工具从自身位置向上两级解析仓库根目录，不依赖调用目录。campaign.json仍专属BE；中点工况在basilisk_midpoint/campaign.json。

已启动的本地守护继续使用已加载代码，无需重启；后续守护与12小时报告使用新路径。旧服务器任务、检查点、场与完整历史不移动。两套数值核心保留原路径，便于核对运行指纹。

原首页详记保存在 [be_method_and_results.md](be_method_and_results.md)，新首页只保留版本与入口导航。
