# Admix - 祖先成分分析工具

[![Build Status](https://travis-ci.org/stevenliuyi/admix.svg?branch=master)](https://travis-ci.org/stevenliuyi/admix)
[![PyPI version](https://badge.fury.io/py/admix.svg)](https://badge.fury.io/py/admix)

Admix 是一个用于计算祖先成分（混合比例）的简单工具，支持来自各种DNA测试供应商（如 [23andme](https://www.23andme.com/) 和 [AncestryDNA](https://www.ancestry.com/dna/)）的SNP原始数据。

## 功能特性

- **多供应商支持**: 支持23andme、AncestryDNA、FTDNA、WeGene、MyHeritage等主流DNA测试供应商
- **多种分析模型**: 支持30+种公开的混合分析模型
- **基于位置的分析**: 新增基于染色体位置的高精度分析功能
- **批量分析**: 支持同时分析多个模型
- **中英文支持**: 支持中英文种群名称显示
- **结果排序**: 支持按比例降序排列结果
- **文件输出**: 支持将结果保存到文件

## 目录

- [安装](#安装)
  - [从GitHub安装](#从github安装)
  - [从PyPI安装](#从pypi安装)
- [使用方法](#使用方法)
  - [基本用法](#基本用法)
  - [基于位置的分析](#基于位置的分析)
  - [高级选项](#高级选项)
- [输出示例](#输出示例)
- [常见问题](#常见问题)
- [原始数据格式](#原始数据格式)
- [分析模型](#分析模型)
- [实现原理](#实现原理)

## 安装

### 从GitHub安装

您可以使用 `pip` 直接从GitHub仓库安装Admix：

```bash
pip install git+https://github.com/TheYtree/admix
```

### 从PyPI安装

您也可以从 [PyPI](https://pypi.python.org/pypi/admix) 安装Admix：

```bash
pip install admix
```

注意：由于大小限制，PyPI上的包只包含五个模型（`K7b`、`K12b`、`globe13`、`world9` 和 `E11`）。如果您需要所有模型，请从GitHub仓库安装。

## 使用方法

### 基本用法

假设您已经下载了23andme原始数据并将其放在当前目录中，文件名为 `my_raw_data.txt`。您可以通过指定计算模型（本例中为 `K7b`）来执行混合计算：

```bash
admix -f my_raw_data.txt -v 23andme -m K7b
```

您也可以设置多个模型进行计算：

```bash
admix -f my_raw_data.txt -v 23andme -m K7b K12b
```

如果未设置模型，程序将应用所有可用模型：

```bash
admix -f my_raw_data.txt -v 23andme
```

您可以通过更改 `-v` 或 `--vendor` 参数来选择原始数据格式。支持的值列在[原始数据格式](#原始数据格式)部分。

您也可以设置 `-o` 或 `--output` 参数将祖先成分结果写入文件：

```bash
admix -f my_raw_data.txt -v 23andme -o result.txt
```

如果您还没有原始数据，也可以使用程序提供的演示23andme数据文件来测试程序：

```bash
admix -m world9
```

中文用户可以通过开启 `-z` 标志来显示中文种群名称：

```bash
admix -z -m E11
```

此外，您可以使用 `--sort` 标志对比例进行排序，使用 `--ignore-zeros` 标志只显示非零比例。

### 基于位置的分析

新增的基于位置的分析功能提供更高精度的分析结果。使用 `--position-based` 参数启用此功能：

```bash
admix -f my_raw_data.txt -m E11 --position-based
```

基于位置的分析会显示匹配的SNP统计信息，提供更详细的分析过程：

```bash
admix -f my_raw_data.txt -m E11 --position-based -z --sort
```

### 高级选项

- `-t, --tolerance`: 设置优化容差（默认：1e-3）
- `--sort`: 对结果进行降序排序
- `--ignore-zeros`: 只显示非零比例
- `-z, --zhongwen`: 显示中文种群名称
- `--position-based`: 使用基于位置的分析方法（更高精度）

获取更多帮助信息：

```bash
admix -h
```

## 输出示例

### 英文输出

命令：`admix -m K12b`

输出：
```
Gedrosia: 0.06%
Siberian: 3.71%
Northwest African: 0.00%
Southeast Asian: 33.43%
Atlantic Med: 0.07%
North European: 0.00%
South Asian: 0.00%
East African: 0.00%
Southwest Asian: 0.01%
East Asian: 62.72%
Caucasus: 0.00%
Sub Saharan: 0.00%
```

### 中文输出

命令：`admix -m K12b -z`

输出：
```
格德罗西亚: 0.06%
西伯利亚: 3.71%
西北非: 0.00%
东南亚: 33.43%
大西洋地中海: 0.07%
北欧: 0.00%
南亚: 0.00%
东非: 0.00%
西南亚: 0.01%
东亚: 62.72%
高加索: 0.00%
撒哈拉以南非洲: 0.00%
```

### 基于位置的分析输出

命令：`admix -f my_raw_data.txt -m E11 --position-based -z`

输出：
```
Using position-based analysis for higher precision...

Admixture calculation models: E11

Calcuation is started...

读取到 739868 个有效SNP
位置匹配统计: 160026/170288 SNPs 匹配
使用基于位置的分析方法
E11
华东 East Chinese: 47.56%
傣族 South Chinese Dai: 19.87%
彝族 Southwest Chinese Yi: 16.44%
鄂伦春 North Chinese Oroqen: 8.33%
日本 Japanese: 7.25%
非洲 African: 0.41%
美洲 American: 0.16%
```

## 常见问题

**问题：** *为什么我得到的每个种群的估计比例都相同？*

**答案：** 这个包使用 [scipy.optimize.minimize](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html) 优化函数，它有一个参数 `tol` 来控制优化器的终止容差。默认容差设置为 `1e-3`。这在大多数情况下都能正常工作，但有时 `1e-3` 太大，导致过早终止。您可以手动设置更小的容差（比如 `1e-4`）来获得正确的结果，尽管优化器运行时间会更长。您可以使用 `-t` 或 `--tolerance` 标志来做到这一点，例如：

```bash
admix -f my_raw_data.txt -t 1e-4
```

## 原始数据格式

Admix 支持以下DNA测试供应商的原始数据格式，使用 `-v` 或 `--vendor` 参数：

| 参数值 | 供应商 |
| ------ | ------ |
| 23andme | [23andme](https://www.23andme.com/) |
| ancestry | [AncestryDNA](https://www.ancestry.com/dna/) |
| ftdna | [FamilyTreeDNA Family Finder](https://www.familytreedna.com/products/family-finder) |
| ftdna2 | [FamilyTreeDNA Family Finder](https://www.familytreedna.com/products/family-finder) (新格式) |
| wegene | [WeGene](https://www.wegene.com/en/) |
| myheritage | [MyHeritageDNA](https://www.myheritage.com/dna) |

## 分析模型

Admix 支持许多公开可用的混合分析模型。所有计算器文件都是其作者的财产，不受本程序许可证的约束。提供了包含每个模型更多信息的链接。

| 模型值 | 模型名称 | 来源 |
| ------ | -------- | ---- |
| `K7b` | Dodecad K7b | [链接](http://dodecad.blogspot.com/2012/01/k12b-and-k7b-calculators.html) |
| `K12b` | Dodecad K12b | [链接](http://dodecad.blogspot.com/2012/01/k12b-and-k7b-calculators.html) |
| `globe13` | Dodecad globe13 | [链接](http://dodecad.blogspot.com/2012/10/globe13-calculator.html) |
| `globe10` | Dodecad globe10 | [链接](http://dodecad.blogspot.com/2012/10/globe10-calculator.html) |
| `world9` | Dodecad world9 | [链接](http://dodecad.blogspot.com/2011/12/world9-calculator.html) |
| `Eurasia7` | Dodecad Eurasia7 | [链接](http://dodecad.blogspot.com/2011/10/eurasia7-calculator.html) |
| `Africa9` | Dodecad Africa9 | [链接](http://dodecad.blogspot.com/2011/09/africa9-calculator.html) |
| `weac2` | Dodecad weac (West Eurasian cline) 2 | [链接](http://dodecad.blogspot.com/2012/06/weac2-calculator.html) |
| `E11` | E11 | [链接](http://www.ranhaer.com/thread-32241-1-1.html) |
| `K36` | Eurogenes K36 | [链接](http://bga101.blogspot.com/2013/03/eurogenes-k36-at-gedmatch.html) |
| `EUtest13` | Eurogenes EUtest K13 | [链接](http://bga101.blogspot.com/2013/11/updated-eurogenes-k13-at-gedmatch.html) |
| `Jtest14` | Eurogenes Jtest K14 | [链接](http://bga101.blogspot.com/2012/09/eurogenes-ashkenazim-ancestry-test-files.html) |
| `HarappaWorld` | HarappaWorld | [链接](http://www.harappadna.org/2012/05/harappaworld-admixture/) |
| `TurkicK11` | Turkic K11 | [链接](http://www.anthrogenica.com/showthread.php?8817-Turkic-K11-Admixture-Calculator) |
| `KurdishK10` | Kurdish K10 | [链接](http://www.anthrogenica.com/showthread.php?8571-K10-Kurdish-Calculator-Version-1/page6) |
| `AncientNearEast13` | Ancient Near East K13 | [链接](http://www.anthrogenica.com/showthread.php?8193-ancient-DNA-in-the-Gedrosia-Near-East-Neolithic-K13) |
| `K7AMI` | Eurogenes K7 AMI | [链接](http://www.anthrogenica.com/showthread.php?4548-Upcoming-DIY-Eurogenes-K7-amp-K8-Calculator-amp-Oracles-for-tracking-E-Asian-amp-ASI) |
| `K8AMI` | Eurogenes K8 AMI | [链接](http://www.anthrogenica.com/showthread.php?4548-Upcoming-DIY-Eurogenes-K7-amp-K8-Calculator-amp-Oracles-for-tracking-E-Asian-amp-ASI) |
| `MDLPK27` | MDLP K27 | [链接](http://www.anthrogenica.com/showthread.php?4557-Post-MDLP-K27-Results) |
| `puntDNAL` | puntDNAL K12 Ancient World | [链接](http://www.anthrogenica.com/showthread.php?8034-PuntDNAL-K12-Ancient-World-Results) |
| `K47` | LM Genetics K47 | [链接](https://anthrogenica.com/showthread.php?12788-New-K30-K47-World-Calculator) |
| `K7M1` | Tolan K7M1 | [链接](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `K13M2` | Tolan K13M2 | [链接](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `K14M1` | Tolan K14M1 | [链接](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `K18M4` | Tolan K18M4 | [链接](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `K25R1` | Tolan K25R1 | [链接](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `MichalK25`| Michal World K25 | [链接](https://anthrogenica.com/showthread.php?13359-Michal-s-World-K25-calculator) |
| `EastSeaK12` | East Sea K12 | [链接](https://anthrogenica.com/showthread.php?13359-Michal-s-World-K25-calculator) |
| `ProjectLiK11` | Project Li K11 | [链接](https://anthrogenica.com/showthread.php?13359-Michal-s-World-K25-calculator) |
| `ProLi14` | Pro Li 14 | [链接](https://anthrogenica.com/showthread.php?13359-Michal-s-World-K25-calculator) |

## 实现原理

使用最大似然估计（MLE）算法进行祖先成分计算，实现相对简单直接。

设 $F_{nk}$ 为种群 $k$ 的SNP标记 $n$ 的次要等位基因频率，
$l_n^{minor}$ 和 $l_n^{major}$ 分别为标记 $n$ 的次要和主要等位基因，
$G_{ni}$ 为我们感兴趣的个体在标记 $n$ 处的等位基因（$i=1,2$）。

我们的目标是找到个体的混合分数 $q_k$，它使对数似然函数最大化

$$\chi_{{l^{minor}_n}}(G_{ni})j_i\log(F_{nk}q_k)+\chi_{{l^{major}_n}}(G_{ni})j_i\log((J_{nk}-F_{nk})q_k)$$

其中 $\chi$ 是指示函数，$J$ 和 $j$ 是全1矩阵/向量。注意这里隐含了爱因斯坦求和约定。在约束条件 $0 \leq q_k \leq 1$ 和 $\sum {q_k} = 1$ 下，我们可以通过应用优化技术获得混合比例 $q_k$。

### 基于位置的分析原理

基于位置的分析使用染色体位置信息来提高匹配精度：

1. **位置匹配**: 使用染色体位置（chromosome:position）作为唯一标识符
2. **精确匹配**: 只匹配具有相同染色体位置的SNP
3. **统计报告**: 显示匹配的SNP数量和总SNP数量
4. **优化算法**: 使用相同的MLE算法但基于更精确的数据

这种方法减少了由于SNP名称不一致导致的匹配错误，提高了分析结果的准确性。
