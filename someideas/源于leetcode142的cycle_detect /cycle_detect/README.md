# 🧭 环检测算法探索手记：从力扣 142 到「已知最优」

> 源于力扣142
> 闭式地同时定位环长 λ 与入口节点 μ
> 利用kmpnextval数组优化策略跳过重复向前探测环的步骤

## 一、缘起

- 出发点：**力扣 142「环形链表 II」**
- 核心直觉：在一趟遍历里**同时推进并维护环长 λ 与入口 μ**，用几何间距（2 的幂）的锚点，配合「失配后复用已匹配前缀」（KMP `nextval` 思想），避免重复探索环外（tail）与环内节点。

## 二、核心想法

用 O(1) 内存，在确定性迭代轨道上，仅靠**局部探测 + 等值判断**，闭式地定位隐藏的全局参数 μ 与 λ。

- 锚点按 $2^k$ 间距布设：既能框定未知尺度，又能在失配时跳到最长可复用前缀；
- **不重启即可定位入口 μ**——这是和经典 Floyd / Brent 最不一样的一点。

## 三、关键结论（已验证）

**闭式表达**（求值次数 $A$）：

$$r_0=\lceil\log_2(\mu+1)\rceil,\quad r_\lambda=\lceil\log_2\lceil\lambda/3\rceil\rceil,\quad r_{\min}=\max(r_0,r_\lambda)$$

$$A=2(\mu+\lambda)+R,\qquad R=\begin{cases}1,&\lambda\le 2^{r_{\min}}\\ 2^{r_{\min}}-1,&\text{otherwise}\end{cases}$$

（自环 $(0,1)$ 是唯一例外，此时 $A=2$。）

**平均常数**：在随机映射模型下，$\bar c_A\approx 2.180337$，带有约 $1.2\times10^{-3}$ 的微小周期涟漪。

**与经典算法对比**（实验中 100% 不劣于二者）：

| 算法 | 平均常数 c | 相对本文 |
|---|---|---|
| 本文做法 | ≈ 2.18 | 1.00× |
| Brent (1980) | ≈ 3.08 | 0.71× |
| Floyd (1967) | ≈ 3.47 | 0.63× |

## 四、撞壁记录（诚实版）

查新后发现的前人工作，以及它们如何「封住」了新颖性：

- **Floyd 1967 / Brent 1980**：经典 O(1) 环检测；Brent 还能直接、一趟求出 λ。
- **Sedgewick–Szymanski–Yao (SSY) 1982**：时空权衡的最优算法 $n(1+O(1/\sqrt M))$，并明确指出经典算法「除了回到起点重来，没有直接找入口的办法」。
- **Fich 1981/1983**：存储 ≤ M 个值的算法，求值次数下界为 $(\lambda+\mu)\left(1+\frac{1}{M-1}\right)$。
- **Nivasch 2004**：对数内存的栈式算法。

**致命的一击**：把 Fich 下界代入 $M=2$，得到 $2(\lambda+\mu)$——而我的做法正好是 $\approx 2(\mu+\lambda)$、O(1) 内存。

> ⚠️ **我没有「打败」任何东西——我恰好坐在了 M=2 的已知理论最优线上。
> 独立摸到最优，是能力的证明；但它不是新颖性。**

## 五、文件结构
cycle-detection/
├── README.md                # 本文件（完整记录，建议先读这个）
├── nextval_unified.tex      # 论文源码（英文版）
├── nextval_unified_zh.tex   # 论文源码（中文版）
├── gen_grid.cpp             # 生成验证用的 (μ, λ) 网格
├── factor_bench.cpp         # Pollard-ρ 因子分解对比实验
├── mkfigs.py / mkfig_avg.py # 画图脚本
├── verify2.py               # 闭式公式验证脚本
└── figures_nextval/         # 图（PDF）


## 六、参考文献

- R. W. Floyd, 1967 — 龟兔算法（Cycle detection）。
- R. P. Brent, 1980 — *An improved Monte Carlo factorization algorithm*, BIT 20.
- R. Sedgewick, T. G. Szymanski, A. C. Yao, 1982 — *The complexity of finding cycles in periodic functions*, SIAM J. Comput. 11(2):376–390.
- F. E. Fich, 1981/1983 — 周期函数求值次数的下界。
- G. Nivasch, 2004 — 对数内存栈式环检测。
- D. E. Knuth, *TAOCP* Vol. 2, §3.1（习题 3.1-6, 3.1-7）。

