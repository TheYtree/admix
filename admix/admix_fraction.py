import numpy as np
import scipy.optimize as optimize
from raw_data_processing import read_raw_data, read_model
import admix_models


# genotype matches
def genotype_matches(genome_data, snp, major_alleles, minor_alleles):
    g_major1 = [(genome_data.get(snp[i], '-')[0] == a)
                for (i, a) in enumerate(major_alleles)]
    g_major2 = [(genome_data.get(snp[i], '-')[-1] == a)
                for (i, a) in enumerate(major_alleles)]
    g_major = np.array(g_major1).astype(int) + np.array(g_major2).astype(int)

    g_minor1 = [(genome_data.get(snp[i], '-')[0] == a)
                for (i, a) in enumerate(minor_alleles)]
    g_minor2 = [(genome_data.get(snp[i], '-')[-1] == a)
                for (i, a) in enumerate(minor_alleles)]
    g_minor = np.array(g_minor1).astype(int) + np.array(g_minor2).astype(int)

    # 计算匹配的SNP数量
    matched_count = sum(1 for i, s in enumerate(snp) if s in genome_data)
    total_snps = len(snp)

    return g_major, g_minor, matched_count, total_snps


# numerical guard against log(0) -> -inf / nan at the boundaries
_EPS = 1e-12


# log likelihood function
def likelihood(g_major, g_minor, frequency, admixture_fraction):
    # p = 每个SNP上 major 等位基因在混血群体中的频率
    # 裁剪到 (_EPS, 1-_EPS) 避免边界 log(0)；这是原实现的数学缺陷之一：
    # 某成分被顶到 1 时 p 可能出现 0，log(0)=-inf 使似然无定义、优化退化。
    p = np.clip(np.dot(frequency, admixture_fraction), _EPS, 1.0 - _EPS)
    q = 1.0 - p
    l1 = np.dot(g_major, np.log(p))
    l2 = np.dot(g_minor, np.log(q))
    # add minus sign so we can minimize the function instead of maximize
    return -(l1 + l2)


# obtain admixture fraction by applying the maximum likelihood estimation
def admix_fraction(model, raw_data_format, raw_data_file=None, tolerance=1e-3):
    genome_data = read_raw_data(raw_data_format, raw_data_file)
    snp, minor_alleles, major_alleles, frequency = read_model(model)
    g_major, g_minor, matched_count, total_snps = genotype_matches(genome_data, snp, major_alleles,
                                        minor_alleles)

    n_pop = admix_models.n_populations(model)

    # constraints of admixture fractions
    bounds = tuple((0, 1) for i in range(n_pop))
    constraints = ({'type': 'eq', 'fun': lambda af: np.sum(af) - 1})

    # MLE
    # 原实现：单起点(uniform) + 松容差 tol=1e-3（默认）-> SLSQP 过早满足收敛条件，
    # 常停在次优/退化局部解（例如把单一成分顶到 100%），导致同一祖先群体的样本
    # 结果不一致、且 >98% 的退化解还会被上层丢弃。修复：
    #   (1) 收紧 ftol（<=1e-6），让优化真正收敛；
    #   (2) 若仍落到退化解(最大成分 >99%)或未收敛/非有限，再用少量确定起点重试取最优。
    # 注意：本次改动只影响数值求解，模型/频率数据与概率模型本身不变。
    likelihood_func = lambda af: likelihood(g_major, g_minor, frequency, af)
    ftol = min(abs(float(tolerance)), 1e-6)

    def _solve(x0):
        return optimize.minimize(
            likelihood_func, x0, bounds=bounds, constraints=constraints,
            method='SLSQP', options={'maxiter': 800, 'ftol': ftol})

    res = _solve(np.ones(n_pop) / n_pop)

    def _degenerate(r):
        return (r is None or (not np.isfinite(r.fun)) or (not r.success)
                or float(np.max(r.x)) > 0.99)

    if _degenerate(res):
        rng = np.random.default_rng(12345)
        starts = [rng.dirichlet(np.ones(n_pop)) for _ in range(5)]
        for k in range(n_pop):  # 各轴顶点附近
            v = np.full(n_pop, 0.05 / max(n_pop - 1, 1))
            v[k] = 0.95
            starts.append(v)
        for x0 in starts:
            try:
                r = _solve(x0)
            except Exception:
                continue
            if np.isfinite(r.fun) and (not np.isfinite(res.fun) or r.fun < res.fun):
                res = r

    admix_frac = np.clip(res.x, 0, 1)
    s = admix_frac.sum()
    if s > 0:
        admix_frac = admix_frac / s
    else:
        admix_frac = np.ones(n_pop) / n_pop

    # 计算利用率
    utilization_rate = matched_count / total_snps if total_snps > 0 else 0

    return admix_frac, utilization_rate
