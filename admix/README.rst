Admix
=====

Admix 是一个用于计算祖先成分（混合比例）的简单工具，支持来自各种DNA测试供应商（如23andme和AncestryDNA）的SNP原始数据。

功能特性
--------

- **多模型支持**: 支持30多种不同的祖先分析模型，包括K7b、K12b、K36、globe13等
- **多数据格式**: 支持23andme、AncestryDNA等主流DNA测试公司的原始数据格式
- **基于位置的分析**: 提供更高精度的基于位置的祖先分析方法
- **利用率显示**: 显示数据利用率，帮助评估分析质量
- **多语言支持**: 支持英文和中文人口名称显示
- **灵活输出**: 支持控制台输出和文件保存
- **优化选项**: 可调整优化容差、排序和过滤零值比例

安装
----

使用pip安装:

.. code-block:: bash

    pip install admix

或者从源码安装:

.. code-block:: bash

    git clone https://github.com/TheYtree/admix.git
    cd admix
    python setup.py install

使用方法
--------

基本用法:

.. code-block:: bash

    # 使用默认设置分析23andme数据
    admix -f your_genome.txt

    # 指定特定模型
    admix -f your_genome.txt -m K7b K12b K36

    # 使用基于位置的分析（更高精度）
    admix -f your_genome.txt --position-based

    # 显示中文人口名称
    admix -f your_genome.txt -z

    # 查看数据利用率（标准分析和基于位置的分析都支持）
    admix -f your_genome.txt -m K7b K12b

    # 保存结果到文件
    admix -f your_genome.txt -o results.txt

    # 排序并忽略零值
    admix -f your_genome.txt --sort --ignore-zeros

命令行参数
----------

- ``-f, --file``: 原始基因组数据文件路径
- ``-m, --models``: 指定要使用的混合模型（默认使用所有可用模型）
- ``-v, --vendor``: 设置DNA测试供应商（默认: 23andme）
- ``-o, --output``: 将结果保存到文件
- ``-z, --zhongwen``: 显示中文人口名称
- ``-t, --tolerance``: 设置优化容差（默认: 1e-3）
- ``--sort``: 按比例排序混合结果
- ``--ignore-zeros``: 仅显示非零比例
- ``--position-based``: 使用基于位置的分析方法

可用模型
--------

项目包含30多种不同的祖先分析模型:

- **K7b**: 7个祖先成分（南亚、西亚、西伯利亚、非洲、地中海-中东、大西洋波罗的海、东亚）
- **K12b**: 12个祖先成分（格德罗西亚、西伯利亚、西北非、东南亚等）
- **K36**: 36个详细祖先成分（美洲印第安、阿拉伯、亚美尼亚、巴斯克等）
- **globe13**: 13个全球祖先成分
- **globe10**: 10个全球祖先成分
- **E11**: 11个东亚特定祖先成分
- **HarappaWorld**: 16个南亚特定祖先成分
- **AncientNearEast13**: 13个古代近东祖先成分
- **EUtest13**: 13个欧洲测试祖先成分
- **Jtest14**: 14个犹太测试祖先成分
- **TurkicK11**: 11个突厥祖先成分
- **KurdishK10**: 10个库尔德祖先成分
- **Africa9**: 9个非洲特定祖先成分
- **EastSeaK12**: 12个东海祖先成分
- **ProjectLiK11/14**: 11/14个李项目祖先成分
- **ProLi14**: 14个李项目祖先成分
- **MichalK25**: 25个Michal祖先成分
- **MDLPK27**: 27个MDLP祖先成分
- **K47**: 47个详细祖先成分
- **K25R1**: 25个R1祖先成分
- **K18M4**: 18个M4祖先成分
- **K14M1**: 14个M1祖先成分
- **K13M2**: 13个M2祖先成分
- **K8AMI**: 8个AMI祖先成分
- **K7AMI**: 7个AMI祖先成分
- **K7M1**: 7个M1祖先成分
- **puntDNAL**: 12个puntDNAL祖先成分
- **weac2**: 7个weac2祖先成分
- **world9**: 9个世界祖先成分
- **Eurasia7**: 7个欧亚祖先成分

数据格式支持
------------

- **23andme**: 默认支持的格式
- **AncestryDNA**: 支持AncestryDNA的原始数据格式
- **其他格式**: 可通过修改代码支持其他DNA测试公司的数据格式

示例输出
--------

.. code-block:: text

    K7b (88%)
    Atlantic Baltic: 45.23%
    Southern: 28.67%
    West Asian: 15.34%
    East Asian: 8.12%
    Siberian: 2.64%

    K12b (92%)
    North European: 38.45%
    Atlantic Med: 25.67%
    Caucasus: 18.23%
    Gedrosia: 12.34%
    Southwest Asian: 5.31%

利用率说明
----------

利用率表示成功匹配的SNP数量占总SNP数量的百分比，显示在模型名称后的括号中（如K7b (88%)）。

- **高利用率（80%以上）**: 数据质量较好，分析结果更可靠
- **中等利用率（50-80%）**: 数据质量一般，结果仅供参考
- **低利用率（50%以下）**: 数据质量较差，建议检查数据格式或使用其他数据

标准分析和基于位置的分析都支持利用率显示，帮助用户评估分析质量。

依赖项
------

- Python 2.7+ 或 Python 3.5+
- numpy
- scipy

许可证
------

本项目采用 GNU General Public License v3.0 许可证。

更多信息
--------

更多详细信息，请访问: https://github.com/TheYtree/admix

贡献
----

欢迎提交问题报告和功能请求。如果您想贡献代码，请fork项目并提交pull request。
