**Nolting Ch.6 笔记：量子多体基础 —— 有限温微扰论 (Part III：两粒子松原函数——密度关联、极化传播子、有效相互作用、顶角函数与横向自旋磁化率)**

> **目标读者**：已完成 Ch.6 Part I（Matsubara 方法与松原函数、广义 Wick 定理、图规则、环图、单粒子松原函数的一阶图）与 Part II（Dyson 方程与骨架图、Hartree-Fock 近似、二阶骨架图、Hubbard 与 Jellium 模型）的读者。
> **撰写原则**：按物理主题重组，保留全部公式与推导逻辑链；引用内容时直接重写；公式一律不编号，后文引用以名称或大意指代；教材正文提及习题处，以（推导如下）引用框补全结论。
> **对应教材**：Wolfgang Nolting, *Fundamentals of Many-Body Physics* (Springer, 2009), 2nd ed., translated by William D. Brewer, Chapter 6, Section 6.3.
> **本文采用 AI 进行文献整理与内容梳理，所有核心内容在作者指导下完成。**

> **缩写说明：** M — Matsubara / 松原；ret — retarded / 推迟；conn — connected / 连通；IR — irreducible / 不可约；RPA — Random Phase Approximation / 无规相近似；LA — ladder approximation / 梯近似；EoM — Equation of Motion / 运动方程；$\varepsilon$ — 统计因子（$+1$ 玻色、$-1$ 费米）；$f_\varepsilon$ — 广义 Fermi-Dirac/Bose-Einstein 函数（$f_-(x)$ 为费米分布）；$G^{0,M}$ — 自由单粒子松原函数；$G^M$ — 完整（敷饰）单粒子松原函数；$S_{\mathbf{p}\sigma}(E)$ — 谱密度.

本文承接 Part II 的单粒子自能概念，把"部分求和"的思想推进到**两粒子**层次。叙事线是：**两粒子松原函数（密度关联）的定义与图规则 → 通过"极化传播子"块化（Dyson 方程与 RPA）→ 通过"有效相互作用"块化（动态屏蔽）→ 通过"顶角函数"块化（梯近似）→ 一个具体应用（横向自旋磁化率）**。教材中所有公式与推导均完整收录；习题 6.3.1–6.3.5 在正文提及处补全。需要特别指出：$T\neq 0$ 的松原图与 $T=0$ 的因果格林函数图**结构完全相同**（这是 Matsubara 形式体系最重要的性质之一），因此本篇每个概念都可以与 Ch.5 §5.6 中 $T=0$ 的对应物（极化传播子、有效相互作用、顶角函数）对照——那套图在本篇中几乎原样复用，只是传播子换成松原函数。

---

# 绪论：两粒子函数的"图块化"

单粒子自能概念带来巨大简化：把无穷级数中的不可约部分提取出来，再通过 Dyson 方程求和。但部分求和并非只有自能这一种方式——Ch.5 §5.6 在 $T=0$ 情形下已经示范了另外三种"图块"：**极化传播子**（把粒子-空穴气泡链求和）、**有效相互作用**（把极化插入到相互作用线里求和）、**顶角函数**（把粒子-空穴散射过程块化）。这些思想可以几乎原封不动地搬进 $T\neq 0$ 的 Matsubara 形式体系，因为：

- $T\neq 0$ 的松原图与 $T=0$ 的因果格林函数图结构相同；
- 松原函数的 Dyson 型方程在结构上与 $T=0$ 完全相同，只需把自由传播子换成自由松原函数、自能换成不可约贡献之和。

本篇以**密度关联函数** $\langle\!\langle \rho_{\mathbf{q}};\rho_{\mathbf{q}}^+\rangle\!\rangle$ 为主线展开，它通过（普通 Fourier 变换引入的）库仑矩阵元 $v(\mathbf{q})$ 与介电函数直接相连，是处理屏蔽、集体激发（等离激元）的标准工具。三种图块各自对应一类无穷部分和的求和，并互为工具：极化传播子进入有效相互作用与密度关联的 Dyson 方程，顶角函数则给出极化传播子的另一种表示（梯近似）。最后，**横向自旋磁化率**展示了两粒子松原函数的一个特殊情形——当外接传播子携带不同自旋时，原本复杂的 Dyson 求和自动坍缩为极化传播子本身，Hubbard 模型的梯近似可以精确求解。

---

# 主题一：密度关联函数与介电函数

## 密度算符、推迟密度关联与介电函数

密度关联是第三章引入过的推迟格林函数，其定义基于（伴随）**密度算符**：

$$
\rho_{\mathbf{q}} = \sum_{\mathbf{k}\sigma} a_{\mathbf{k}\sigma}^{+}a_{\mathbf{k}+\mathbf{q}\sigma};\qquad \rho_{\mathbf{q}}^{+} \equiv \rho_{-\mathbf{q}},
$$

即"密度涨落 $\rho_{\mathbf{q}}$ 对密度涨落 $\rho_{\mathbf{q}}^+$ 的响应"。以 Jellium 模型为例，第三章展示了它与物理上重要的**介电函数**的密切联系：

$$
\frac{1}{\varepsilon(\mathbf{q},E)} = 1 + \frac{1}{\hbar}\,v(\mathbf{q})\,\langle\!\langle \rho_{\mathbf{q}};\rho_{\mathbf{q}}^{+}\rangle\!\rangle_E^{\text{ret}},
$$

其中 $v(\mathbf{q}) = e^{2}/\varepsilon_{0}Vq^{2}$ 是 Jellium 模型的相关库仑矩阵元。第三章的推导表明，这个因子来自某些正规 Fourier 变换，与 Jellium 模型中的相互作用无关，因此上式应当**普遍成立**。唯一的前提是：平衡态下电子与离子系统的电荷密度精确相互补偿，且外部扰动电荷只作用于（响应更快的）电子子系统。

> **外加扰动电荷与诱导电荷**（第三章，内容转写）：外部"扰动电荷" $\rho_{\text{ext}}(\mathbf{q},E)$ 与其诱导的电荷密度 $\rho_{\text{ind}}(\mathbf{q},E)$ 之间有关系
> $$
> \rho_{\text{ind}}(\mathbf{q},E) = \left(\frac{1}{\varepsilon(\mathbf{q},E)} - 1\right)\rho_{\text{ext}}(\mathbf{q},E),
> $$
> 由此可以区分出若干有趣的极限情形：
> - $\varepsilon(\mathbf{q},E) \gg 1$：扰动电荷被**几乎完全屏蔽**；
> - $\varepsilon(\mathbf{q},E) \to 0 \iff \langle\!\langle \rho_{\mathbf{q}};\rho_{\mathbf{q}}^{+}\rangle\!\rangle_E^{\text{ret}}$ 出现奇异性：任意小的扰动电荷都引起电荷密度的**有限涨落** $\implies$ 系统的元激发，即 **plasmons**（等离激元）$E = E(\mathbf{q})$。

推迟密度关联（连同介电函数）在第四章已经用运动方程方法近似算过；本节的任务是展示如何用**图解 Matsubara 形式体系**来计算它。为此先固定图的记号约定（Fig 6.26）：

![Fig 6.26：顶点记号](./fig有限温微扰论3两粒子松原函数/fig626.png)

> **顶点记号**（Fig 6.26，内容转写）：
> 1. 每个顶点处**自旋守恒**：$\sigma_k = \sigma_n$；$\sigma_l = \sigma_m$；
> 2. 顶点处**动量守恒**：$\mathbf{k} + \mathbf{l} = \mathbf{m} + \mathbf{n}$；
> 3. 相互作用矩阵元至多依赖于**转移动量**：$\mathbf{q} \equiv \mathbf{k} - \mathbf{n} = \mathbf{m} - \mathbf{l}$。

这些假设对我们关心的绝大多数模型都成立，至少对 Hubbard 模型与 Jellium 模型成立。在此约定下，一般的相互作用矩阵元化为

$$
v(kl;nm) \rightarrow v_{\sigma_k\sigma_l}(\mathbf{q} = \mathbf{k}-\mathbf{n})\,\delta_{\mathbf{k}+\mathbf{l},\mathbf{m}+\mathbf{n}}\,\delta_{\sigma_k\sigma_n}\,\delta_{\sigma_l\sigma_m},
$$

其中两个具体模型的 $v_{\sigma_k\sigma_l}(\mathbf{q})$ 分别为

$$
v_{\sigma_k\sigma_l}(\mathbf{q}) = \left\{\begin{array}{ll} v(\mathbf{q}) & \text{(Jellium)}\\[2mm] \dfrac{U}{N}\,\delta_{\sigma_k,-\sigma_l} & \text{(Hubbard)} \end{array}\right.
$$

## 两粒子松原函数及其图展开

我们的出发点是如下的**两粒子松原函数**：

$$
D_{\mathbf{q}}(E_0) = \langle\!\langle \rho_{\mathbf{q}};\rho_{\mathbf{q}}^{+}\rangle\!\rangle_{E_0}^{M} = \int_0^{\hbar\beta} d\tau\, e^{\frac{\mathrm{i}}{\hbar}E_0(\tau-\tau')}\,D_{\mathbf{q}}(\tau-\tau'),
$$

$$
D_{\mathbf{q}}(\tau-\tau') = -\left\langle T_{\tau}\left(\rho_{\mathbf{q}}(\tau)\,\rho_{\mathbf{q}}^{+}(\tau')\right)\right\rangle .
$$

算符仍处于修正 Heisenberg 表示。过渡到修正 Dirac 表示的方式与单粒子松原函数的情形完全相同（与编时算符 $T_\tau$、演化算符 $U(\hbar\beta,0)$ 的标准做法一致），得到

$$
D_{\mathbf{q}}(\tau-\tau') = -\sum_{\substack{\mathbf{k}\mathbf{p}\\ \sigma\sigma'}}\frac{\left\langle T_{\tau}\left(U(\hbar\beta,0)\,a_{\mathbf{k}\sigma}^{+}(\tau)a_{\mathbf{k}+\mathbf{q}\sigma}(\tau)\,a_{\mathbf{p}\sigma'}^{+}(\tau')a_{\mathbf{p}-\mathbf{q}\sigma'}(\tau')\right)\right\rangle^{(0)}}{\langle U(\hbar\beta,0)\rangle^{(0)}} \equiv \sum_{\sigma\sigma'}D_{\mathbf{q}\sigma\sigma'}(\tau-\tau').
$$

此后所有算符都理解为 Dirac 表示，下标 D 一律省略。对后面要讨论的部分求和而言，先对 **自旋分辨的密度关联** $D_{\mathbf{q}\sigma\sigma'}(\tau-\tau')$ 做图解分析是方便的；最终只需对 $\sigma$ 与 $\sigma'$ 求和即得真正关心的密度关联 $D_{\mathbf{q}}$。

> **连通图定律**（Part I，内容转写）：由于"连通图定律"，上式分母 $\langle U(\hbar\beta,0)\rangle^{(0)}$ 恰好抵消，求值时只需对所有**连通开放图**求和：
> $$
> D_{\mathbf{q}\sigma\sigma'}(\tau-\tau') = -\sum_{\mathbf{k}\mathbf{p}}\left\langle T_{\tau}\left(U(\hbar\beta,0)\,a_{\mathbf{k}\sigma}^{+}(\tau)a_{\mathbf{k}+\mathbf{q}\sigma}(\tau)\,a_{\mathbf{p}\sigma'}^{+}(\tau')a_{\mathbf{p}-\mathbf{q}\sigma'}(\tau')\right)\right\rangle_{\text{conn. open}}^{(0)}.
> $$

每个被加项

$$
\widehat{D}_{\mathbf{k}\mathbf{p}\mathbf{q}\sigma\sigma'}(\tau-\tau') = -\Big\langle T_{\tau}\left(U(\hbar\beta,0)\,a_{\mathbf{k}\sigma}^{+}(\tau)a_{\mathbf{k}+\mathbf{q}\sigma}(\tau)\,a_{\mathbf{p}\sigma'}^{+}(\tau')a_{\mathbf{p}-\mathbf{q}\sigma'}(\tau')\right)\Big\rangle_{\substack{\text{conn.}\\ \text{open}}}^{(0)}
$$

对应一簇开放、连通的图，在 $\tau$ 和 $\tau'$ 处各带两条外线（一进一出），如图 6.27 所示。对密度关联而言（见 $\rho_{\mathbf{q}}$ 的定义），这些外传播子携带**相同的自旋指标** $\sigma$；在 §主题五 中将见到外传播子携带不同自旋的例子。

![Fig 6.27–6.29：密度关联图的一般结构（6.27）、开放的不可连通图结构（6.28）与零阶、一阶密度关联图的例子（6.29）](./fig有限温微扰论3两粒子松原函数/fig627fig628fig629.png)

与 §5.6.1 一样，可以立即确信：**所有开放图都自动是连通的**。由于顶点处动量守恒的假设，像图 6.28 那样的非连通图结构只有在 $\mathbf{q} = 0$ 时才可能出现；在 Jellium 模型中这类图由于 $v(\mathbf{0}) = 0$ 而不作贡献，况且它们本身也乏味——$\mathbf{q} = \mathbf{0}$ 时密度算符就是粒子数算符 $\widehat{N} = \sum_{\mathbf{k}\sigma}a_{\mathbf{k}\sigma}^{+}a_{\mathbf{k}\sigma}$。图 6.29 给出零阶与一阶图的一些例子（阶数仍由相互作用线数目决定），它们与 §5.6.1 的 $T=0$ 图在结构上自然完全相同。

## 能域表示与图规则

求值时采用能域表示是方便的：

$$
\widehat{D}_{\mathbf{k}\mathbf{p}\mathbf{q}\sigma\sigma'}(\tau-\tau') = \frac{1}{\hbar\beta}\sum_{E_0}e^{-\frac{\mathrm{i}}{\hbar}E_0(\tau-\tau')}\,\widehat{D}_{\mathbf{k}\mathbf{p}\mathbf{q}\sigma\sigma'}(E_0).
$$

![Fig 6.30：密度关联图的能域表示](./fig有限温微扰论3两粒子松原函数/fig630.png)

与 Part II 中单粒子松原函数的情形相同，这一变换首先带来顶点处的**能量守恒**。外线需要特殊处理：它们对图的贡献可分解为三个因子（与单粒子情形中导出的分解方式相同）——第一因子进入顶点处能量守恒；第二、第三因子给出如下形式的贡献：

$$
\begin{array}{c}
(\mathbf{k},E_1):\quad \left(-G_{\mathbf{k}\sigma}^{0,M}(E_1)\right)\left(\frac{1}{\sqrt{\hbar\beta}}\exp\left(\frac{\mathrm{i}}{\hbar}E_1\tau\right)\right)\\[2mm]
(\mathbf{k}+\mathbf{q},E_2):\quad \left(-G_{\mathbf{k}+\mathbf{q}\sigma}^{0,M}(E_2)\right)\left(\frac{1}{\sqrt{\hbar\beta}}\exp\left(-\frac{\mathrm{i}}{\hbar}E_2\tau\right)\right)\\[2mm]
(\mathbf{p},E_3):\quad \left(-G_{\mathbf{p}\sigma'}^{0,M}(E_3)\right)\left(\frac{1}{\sqrt{\hbar\beta}}\exp\left(\frac{\mathrm{i}}{\hbar}E_3\tau'\right)\right)\\[2mm]
(\mathbf{p}-\mathbf{q},E_4):\quad \left(-G_{\mathbf{p}-\mathbf{q}\sigma'}^{0,M}(E_4)\right)\left(\frac{1}{\sqrt{\hbar\beta}}\exp\left(-\frac{\mathrm{i}}{\hbar}E_4\tau'\right)\right).
\end{array}
$$

把图核的贡献记为 $A_{\mathbf{k}\mathbf{p}\mathbf{q}\sigma\sigma'}(E_1\ldots E_4)$（图 6.30），合起来得到

$$
-\widehat{D}_{\mathbf{k}\mathbf{p}\mathbf{q}\sigma\sigma'}(\tau-\tau') = \frac{\varepsilon}{\hbar^2\beta^2}\sum_{E_1\ldots E_4}\left(-G_{\mathbf{k}\sigma}^{0,M}(E_1)\right)\left(-G_{\mathbf{k}+\mathbf{q}\sigma}^{0,M}(E_2)\right)\left(-G_{\mathbf{p}\sigma'}^{0,M}(E_3)\right)\left(-G_{\mathbf{p}-\mathbf{q}\sigma'}^{0,M}(E_4)\right)A_{\mathbf{k}\mathbf{p}\mathbf{q}\sigma\sigma'}(E_1\ldots E_4)\,e^{-\frac{\mathrm{i}}{\hbar}\left((E_2-E_1)\tau-(E_3-E_4)\tau'\right)},
$$

因子 $\varepsilon$ 来自环规则。由于所考虑的哈密顿量不含显式时间依赖，上述表达式只能依赖于时间差 $\tau - \tau'$（"时间均匀性"），这要求

$$
E_2 - E_1 \stackrel{!}{=} E_3 - E_4 \equiv E_0.
$$

作为两个松原能量之差，$E_0$ 无论如何都是**玻色型**的（对费米子：$e^{-\mathrm{i}\beta E_1} = e^{-\mathrm{i}\beta E_2} = -1 \Rightarrow e^{-\mathrm{i}\beta E_0} = +1$）。作代换

$$
E_1 = E;\qquad E_2 = E + E_0;\qquad E_3 = E';\qquad E_4 = E' - E_0,
$$

再对上面的能域表示作 Fourier 变换，得到

$$
-\widehat{D}_{\mathbf{k}\mathbf{p}\mathbf{q}\sigma\sigma'}(E_0) = \frac{\varepsilon}{\hbar\beta}\sum_{E,E'}\left(-G_{\mathbf{k}\sigma}^{0,M}(E)\right)\left(-G_{\mathbf{k}+\mathbf{q}\sigma}^{0,M}(E+E_0)\right)\left(-G_{\mathbf{p}\sigma'}^{0,M}(E')\right)\left(-G_{\mathbf{p}-\mathbf{q}\sigma'}^{0,M}(E'-E_0)\right)A_{\mathbf{k}\mathbf{p}\mathbf{q}\sigma\sigma'}(E,E',E_0),
$$

并利用

$$
D_{\mathbf{q}\sigma\sigma'}(E_0) = \sum_{\mathbf{k}\mathbf{p}}\widehat{D}_{\mathbf{k}\mathbf{p}\mathbf{q}\sigma\sigma'}(E_0),
$$

即可表述自旋分辨密度关联 $-D_{\mathbf{q}\sigma\sigma'}(E_0)$ 的**图规则**。

![Fig 6.31：能域中的自旋分辨密度关联图](./fig有限温微扰论3两粒子松原函数/fig631.png)

> **图规则（$-D_{\mathbf{q}\sigma\sigma'}(E_0)$，内容转写）**：求和对象为所有带四条外连续线（图 6.31）的开放连通图。$n$ 阶图（$n$ 个顶点）按下列规则求值：
> 1. 顶点 $\Longleftrightarrow \dfrac{1}{\hbar\beta}\,v_{\sigma_k\sigma_l}(\mathbf{q})\,\delta_{E_k+E_l,E_m+E_n}\,\delta_{\mathbf{k}+\mathbf{l},\mathbf{m}+\mathbf{n}}\,\delta_{\sigma_k\sigma_n}\,\delta_{\sigma_l\sigma_m}$；$(\mathbf{q} = \mathbf{k}-\mathbf{n})$（见 $v_{\sigma_k\sigma_l}(\mathbf{q})$ 的两个模型表达式）；
> 2. 内部连续线（传播或非传播）$\Longleftrightarrow -G_{\mathbf{n}\sigma_n}^{0,M}(E_n) = \dfrac{-\hbar}{\mathrm{i}E_n - \varepsilon(\mathbf{n}) + \mu}$；
> 3. 非传播线附加因子 $\Longleftrightarrow \exp\left(\dfrac{\mathrm{i}}{\hbar}E_n\cdot 0^{+}\right)$；
> 4. 外接（传播子）$\Longleftrightarrow$ 左侧 $\left(-G_{\mathbf{k}\sigma}^{0,M}(E)\right)\left(-G_{\mathbf{k}+\mathbf{q}\sigma}^{0,M}(E+E_0)\right)$；右侧 $\left(-G_{\mathbf{p}\sigma'}^{0,M}(E')\right)\left(-G_{\mathbf{p}-\mathbf{q}\sigma'}^{0,M}(E'-E_0)\right)$；
> 5. 对所有"内部"波数、自旋与松原能量求和，即对 $\mathbf{k},\mathbf{p},E,E'$ 求和，**不对** $\mathbf{q},E_0,\sigma,\sigma'$ 求和；
> 6. 因子 $\dfrac{1}{\hbar\beta}\left(-\dfrac{1}{\hbar}\right)^n\varepsilon^{S}$；$S$ = 环数（loop number）。

规则 6 中多出的因子 $1/\hbar\beta$ 来自**四条外接**（单粒子松原函数只有两条外接）。要得到真正的密度关联 $-D_{\mathbf{q}}(E_0)$，只需再把 $-D_{\mathbf{q}\sigma\sigma'}(E_0)$ 对 $\sigma$ 和 $\sigma'$ 求和。

## 零阶近似：粒子-空穴气泡

作为形式体系的第一项应用，显式计算最低阶（零阶）的密度关联。求值图 6.32 中的图，此时必须有 $\sigma = \sigma'$：

![Fig 6.32：最低（零）阶密度关联图](./fig有限温微扰论3两粒子松原函数/fig632.png)

$$
-\hbar\Lambda_{\mathbf{q}}^{(0)}(E_0) \equiv -D_{\mathbf{q}}^{(n=0)}(E_0) = -\sum_{\sigma\sigma'}D_{\mathbf{q}\sigma\sigma'}^{(n=0)}(E_0)\,\delta_{\sigma\sigma'} = \frac{\varepsilon}{\hbar\beta}\sum_{\mathbf{k},E,\sigma}G_{\mathbf{k}+\mathbf{q}}^{0,M}(E+E_0)\,G_{\mathbf{k}}^{0,M}(E).
$$

这里引入记号 $\hbar\Lambda_{\mathbf{q}}^{(0)}(E_0)$，是为了预告"极化传播子"（下一主题）——零阶密度关联正是零阶极化传播子。自由松原函数与自旋无关，对 $\sigma$ 求和只给出因子 $2$：

$$
\hbar\Lambda_{\mathbf{q}}^{(0)}(E_0) = -2\varepsilon\hbar^2\sum_{\mathbf{k}}I_{\mathbf{k}}(\mathbf{q}) = -2\varepsilon\hbar^2\cdot\frac{1}{\hbar\beta}\sum_{\mathbf{k}}\sum_{E}\frac{1}{\mathrm{i}E - \varepsilon(\mathbf{k}) + \mu}\cdot\frac{1}{\mathrm{i}(E+E_0) - \varepsilon(\mathbf{k}+\mathbf{q}) + \mu}.
$$

> **松原频率求和（围道积分方法）**（Part I，内容转写）：对松原能量 $E_n$ 的求和可用（Part I 建立的）围道积分公式
> $$
> \frac{1}{\hbar\beta}\sum_{E_n}H(\mathrm{i}E_n) = \frac{\varepsilon}{2\pi\mathrm{i}\hbar}\oint_{C'}dE\,\frac{H(E)}{e^{\beta E}-\varepsilon},
> $$
> 其中 $C'$ 为 Part I 图 6.3 所示的路径，数学上沿负方向绕行，$\varepsilon$ 为统计因子。其来源是：$1/(e^{\beta z}-\varepsilon)$ 在费米/玻色松原能量处具有一阶极点（残数 $\varepsilon/\beta$），于是被积函数的极点即 $H(E)$ 在复平面上的极点，残数定理给出结果。

对本问题，围道积分给出

$$
I_{\mathbf{k}}(\mathbf{q}) = \frac{\varepsilon}{2\pi\mathrm{i}\hbar}\oint_{C'}dE\,\frac{1}{e^{\beta E}-\varepsilon}\cdot\frac{1}{(E-\varepsilon(\mathbf{k})+\mu)(E+\mathrm{i}E_0-\varepsilon(\mathbf{k}+\mathbf{q})+\mu)}.
$$

路径 $C'$ 沿数学负方向绕行。被积函数有两个极点 $E_1 = \varepsilon(\mathbf{k}) - \mu$ 与 $E_2 = \varepsilon(\mathbf{k}+\mathbf{q}) - \mu - \mathrm{i}E_0$，残数定理给出

$$
I_{\mathbf{k}}(\mathbf{q}) = \frac{-2\varepsilon\pi\mathrm{i}}{2\pi\mathrm{i}\hbar}\left(\frac{1}{e^{\beta(\varepsilon(\mathbf{k})-\mu)}-\varepsilon}\cdot\frac{1}{\varepsilon(\mathbf{k})-\mu+\mathrm{i}E_0-\varepsilon(\mathbf{k}+\mathbf{q})+\mu} + \frac{1}{e^{\beta(\varepsilon(\mathbf{k}+\mathbf{q})-\mu-\mathrm{i}E_0)}-\varepsilon}\cdot\frac{1}{\varepsilon(\mathbf{k}+\mathbf{q})-\mu-\mathrm{i}E_0-\varepsilon(\mathbf{k})+\mu}\right).
$$

如前所述 $E_0$ 是玻色型的，故 $\exp(-\mathrm{i}\beta E_0) = +1$，从而 $f_\varepsilon$ 出现：

$$
I_{\mathbf{k}}(\mathbf{q}) = -\frac{\varepsilon}{\hbar}\frac{1}{\mathrm{i}E_0 + \varepsilon(\mathbf{k}) - \varepsilon(\mathbf{k}+\mathbf{q})}\left(f_\varepsilon(\varepsilon(\mathbf{k})-\mu) - f_\varepsilon(\varepsilon(\mathbf{k}+\mathbf{q})-\mu)\right) = -\frac{\varepsilon}{\hbar}\frac{\langle n_{\mathbf{k}}\rangle^{(0)} - \langle n_{\mathbf{k}+\mathbf{q}}\rangle^{(0)}}{\mathrm{i}E_0 + \varepsilon(\mathbf{k}) - \varepsilon(\mathbf{k}+\mathbf{q})},
$$

其中 $f_\varepsilon(x) = 1/(e^{\beta x}-\varepsilon)$ 是广义 Fermi-Dirac/Bose-Einstein 函数，$\langle n_{\mathbf{k}}\rangle^{(0)}$ 为无相互作用系统的占位数。最终得到

$$
\boxed{\;\Lambda_{\mathbf{q}}^{(0)}(E_0) = 2\sum_{\mathbf{k}}\frac{\langle n_{\mathbf{k}}\rangle^{(0)} - \langle n_{\mathbf{k}+\mathbf{q}}\rangle^{(0)}}{\mathrm{i}E_0 + \varepsilon(\mathbf{k}) - \varepsilon(\mathbf{k}+\mathbf{q})}\;}
$$

用这个结果，密度关联在最低级近似下已经确定。把它代入介电函数关系式，并做解析延拓 $\mathrm{i}E_0 \to E + \mathrm{i}0^{+}$（推迟函数的过渡），得到介电函数的近似表达式：

$$
\frac{1}{\varepsilon^{(0)}(\mathbf{q},E)} = 1 + 2v(\mathbf{q})\sum_{\mathbf{k}}\frac{\langle n_{\mathbf{k}}\rangle^{(0)} - \langle n_{\mathbf{k}+\mathbf{q}}\rangle^{(0)}}{E + \mathrm{i}0^{+} + \varepsilon(\mathbf{k}) - \varepsilon(\mathbf{k}+\mathbf{q})}.
$$

介电函数的零点代表系统的元激发。把这个表达式用于 Jellium 模型：其零点恰好对应**粒子-空穴激发**，不出现其他零点——例如**完全没有集体激发（等离激元）的迹象**。这正是下一主题要补救的：通过对无穷部分级数求和（RPA），等离激元会从介电函数的额外零点中涌现出来。

---

# 主题二：极化传播子与无规相近似

## 极化贡献与不可约极化贡献

与单粒子松原函数的 Dyson 方程（$G^M = G^{0,M} + G^{0,M}\frac{1}{\hbar}\Sigma G^M$，内容见 Part II）相平行，对密度关联也可以分出无穷部分级数。类比 $T=0$ 的 §5.6.1：

> **定义（自旋分辨极化贡献）**：自旋分辨极化贡献 = $-D_{\mathbf{q}\sigma\sigma'}(E_0)$ 的图贡献，它具有两条外接的相互作用线，并且此外各有一条入射、一条出射的传播子（图 6.33）。

![Fig 6.33–6.34：自旋分辨极化贡献的一般图结构（6.33）与（不可约及可约）极化贡献的例子（6.34）](./fig有限温微扰论3两粒子松原函数/fig633fig634.png)

可以立刻看出，$-D_{\mathbf{q}\sigma\sigma'}(E_0)$ 展开中的全部图都是极化贡献（对比图 6.31 与图 6.33）；例子见图 6.34。下一步定义

> **定义（不可约自旋分辨极化贡献）**：不可约自旋分辨极化贡献 = 不能通过**切断一条相互作用线**而分解为两个低阶独立极化贡献图的极化贡献。

图 6.34 中第三个图显然是可约的，而前两个图不可约。可约图的一般形式如图 6.35 所示，由三个结构单元组成：部分 (a) 代表某个不可约自旋分辨极化贡献，部分 (b) 是一条相互作用线，部分 (c) 代表某个低阶的（可约或不可约）自旋分辨密度关联图。显然，若在 (a) 中遍历所有不可约自旋分辨极化贡献、在 (c) 中遍历所有自旋分辨密度关联图，就得到全部这类图。这引出定义：

![Fig 6.35–6.37：可约极化贡献的三单元结构（6.35）、极化传播子记号（6.36）与自旋分辨密度关联的 Dyson 方程（6.37）](./fig有限温微扰论3两粒子松原函数/fig635fig636fig637.png)

> **定义（自旋分辨极化传播子）**：自旋分辨极化传播子 $-\hbar\Lambda_{\mathbf{q}\sigma\sigma'}(E_0)$ = 所有不可约自旋分辨极化贡献之和。图式上用图 6.36 的记号表示。

## 自旋分辨密度关联的 Dyson 方程

借助自旋分辨极化传播子，可以对真正关心的自旋分辨密度关联写下 Dyson 方程（图 6.37），它与单粒子松原函数的 Dyson 方程等价：

$$
\boxed{\;D_{\mathbf{q}\sigma\sigma'}(E_0) = \hbar\Lambda_{\mathbf{q}\sigma\sigma'}(E_0) + \sum_{\sigma''\sigma'''}\Lambda_{\mathbf{q}\sigma\sigma''}(E_0)\,v_{\sigma''\sigma'''}(\mathbf{q})\,D_{\mathbf{q}\sigma'''\sigma'}(E_0)\;}
$$

图 6.37 中顶点的记号仍需说明。按图规则 1，顶点"通常"携带因子

$$
\frac{1}{\hbar\beta}\,v_{\sigma\sigma'}(\mathbf{q})\,\delta_{E_1+E_2,E_3+E_4}\,\delta_{\mathbf{k}+\mathbf{l},\mathbf{m}+\mathbf{n}}\,\delta_{\sigma_k\sigma_n}\,\delta_{\sigma_l\sigma_m},
$$

Kronecker 符号可以略去，因为能量、自旋与动量守恒已经直接写进图的记号。若把 $-\hbar\Lambda_{\mathbf{q}\sigma\sigma''}$ 的右接与 $\left(-D_{\mathbf{q}\sigma'''\sigma'}(E_0)\right)$ 的左接归入"外接"，则按规则 5 它们贡献一个 $1/\hbar\beta$ 因子（"内部"传播子不携带这个因子），于是这个因子不再需要由所画出的顶点提供。最后，图的阶 $n$ 由顶点数给出，带来因子 $(-1/\hbar)^n$；如果按图 6.37 的 Dyson 方程那样抽出一个特殊顶点，它必须伴随因子 $(-1/\hbar)$。于是剩下的是上述 Dyson 方程。

把剩余项在自旋空间中组合为 $2\times 2$ 矩阵，模型相关的相互作用矩阵 $\widetilde{V}(\mathbf{q})$ 的元素即 $v_{\sigma\sigma'}(\mathbf{q})$：

$$
\widetilde{V}(\mathbf{q}) = \Big(v_{\sigma\sigma'}(\mathbf{q})\Big)_{\substack{\sigma=\uparrow,\downarrow\\ \sigma'=\uparrow,\downarrow}},
$$

则 Dyson 方程可以写成矩阵方程

$$
\widetilde{D}_{\mathbf{q}}(E_0) = \hbar\widetilde{\Lambda}_{\mathbf{q}}(E_0) + \widetilde{\Lambda}_{\mathbf{q}}(E_0)\,\widetilde{V}(\mathbf{q})\,\widetilde{D}_{\mathbf{q}}(E_0),
$$

其解为

$$
\boxed{\;\widetilde{D}_{\mathbf{q}}(E_0) = \frac{\hbar\,\widetilde{\Lambda}_{\mathbf{q}}(E_0)}{\mathbb{1} - \widetilde{\Lambda}_{\mathbf{q}}(E_0)\,\widetilde{V}(\mathbf{q})}\;}
$$

这样，密度关联就完全由极化传播子确定；但注意：解出矩阵方程后，还必须对矩阵 $\widetilde{D}_{\mathbf{q}}(E_0)$ 的所有元素 $D_{\mathbf{q}\sigma\sigma'}(E_0)$ 求和。因此要严格区分 $D_{\mathbf{q}}(E_0)$ 与 $\widetilde{D}_{\mathbf{q}}(E_0)$：

$$
D_{\mathbf{q}}(E_0) = \langle\!\langle \rho_{\mathbf{q}};\rho_{\mathbf{q}}^{+}\rangle\!\rangle_{E_0}^{M} = \sum_{\sigma\sigma'}D_{\mathbf{q}\sigma\sigma'}(E_0).
$$

## Jellium 模型：归约为总极化传播子（推导如下）

由于 Jellium 模型相互作用矩阵的特殊形式

$$
\widetilde{V}(\mathbf{q}) \equiv v(\mathbf{q})\begin{pmatrix}1 & 1\\ 1 & 1\end{pmatrix},
$$

密度关联可以直接用实际的（即非自旋分辨的）极化传播子

$$
\Lambda_{\mathbf{q}}(E_0) = \sum_{\sigma\sigma'}\Lambda_{\mathbf{q}\sigma\sigma'}(E_0)
$$

表达：

$$
D_{\mathbf{q}}(E_0) = \langle\!\langle \rho_{\mathbf{q}};\rho_{\mathbf{q}}^{+}\rangle\!\rangle_{E_0}^{M} = \frac{\hbar\,\Lambda_{\mathbf{q}}(E_0)}{1 - v(\mathbf{q})\Lambda_{\mathbf{q}}(E_0)}.
$$

> **（推导如下：教材习题 6.3.1）** —— 上式有两种证法。
>
> 方法一（在 Dyson 方程中使用特殊相互作用形式）：Jellium 模型中 $v_{\sigma''\sigma'''}(\mathbf{q}) \equiv v(\mathbf{q})$，于是自旋分辨 Dyson 方程为
> $$
> D_{\mathbf{q}\sigma\sigma'}(E_0) = \hbar\Lambda_{\mathbf{q}\sigma\sigma'}(E_0) + v(\mathbf{q})\sum_{\sigma''\sigma'''}\Lambda_{\mathbf{q}\sigma\sigma''}(E_0)\,D_{\mathbf{q}\sigma'''\sigma'}(E_0).
> $$
> 对 $\sigma\sigma'$ 求和，两个求和号立即因子化：
> $$
> \begin{aligned}
> D_{\mathbf{q}}(E_0) &= \hbar\sum_{\sigma\sigma'}\Lambda_{\mathbf{q}\sigma\sigma'}(E_0) + v(\mathbf{q})\left(\sum_{\sigma\sigma''}\Lambda_{\mathbf{q}\sigma\sigma''}(E_0)\right)\left(\sum_{\sigma'''\sigma'}D_{\mathbf{q}\sigma'''\sigma'}(E_0)\right)\\
> &= \hbar\Lambda_{\mathbf{q}}(E_0) + v(\mathbf{q})\,\Lambda_{\mathbf{q}}(E_0)\,D_{\mathbf{q}}(E_0),
> \end{aligned}
> $$
> 解出 $D_{\mathbf{q}}(E_0) = \hbar\Lambda_{\mathbf{q}}(E_0)/\left(1 - v(\mathbf{q})\Lambda_{\mathbf{q}}(E_0)\right)$。
>
> 方法二（矩阵方程直接相乘）：由于
> $$
> \widetilde{V}(\mathbf{q})\widetilde{D}_{\mathbf{q}}(E_0) = v(\mathbf{q})\begin{pmatrix}1&1\\1&1\end{pmatrix}\begin{pmatrix}D_{\uparrow\uparrow}&D_{\uparrow\downarrow}\\ D_{\downarrow\uparrow}&D_{\downarrow\downarrow}\end{pmatrix} = v(\mathbf{q})\begin{pmatrix}D_{\uparrow\uparrow}+D_{\downarrow\uparrow} & D_{\uparrow\downarrow}+D_{\downarrow\downarrow}\\ D_{\uparrow\uparrow}+D_{\downarrow\uparrow} & D_{\uparrow\downarrow}+D_{\downarrow\downarrow}\end{pmatrix},
> $$
> 矩阵 $\widetilde{\Lambda}_{\mathbf{q}}(E_0)\widetilde{V}(\mathbf{q})\widetilde{D}_{\mathbf{q}}(E_0)$ 的四个元素各有同样的求和结构，把所有元素对 $\sigma\sigma'$ 相加后，各项都因子化为总极化传播子与总密度关联的乘积：
> $$
> \sum_{\sigma\sigma'}\left(\widetilde{\Lambda}_{\mathbf{q}}(E_0)\widetilde{V}(\mathbf{q})\widetilde{D}_{\mathbf{q}}(E_0)\right)_{\sigma\sigma'} = v(\mathbf{q})\left(\sum_{\sigma\sigma'}\Lambda_{\mathbf{q}\sigma\sigma'}(E_0)\right)\left(\sum_{\sigma\sigma'}D_{\mathbf{q}\sigma\sigma'}(E_0)\right) = v(\mathbf{q})\Lambda_{\mathbf{q}}(E_0)D_{\mathbf{q}}(E_0),
> $$
> 对矩阵 Dyson 方程逐元素求和即得同一结果。

对物理上重要的介电函数，这意味着

$$
\boxed{\;\varepsilon(\mathbf{q},E_0) = 1 - v(\mathbf{q})\Lambda_{\mathbf{q}}(E_0)\;}
$$

## RPA 介电函数与等离激元

在第一近似中，Jellium 模型的极化传播子 $\Lambda_{\mathbf{q}}(E_0)$ 应替换为零阶表达式 $\Lambda_{\mathbf{q}}^{(0)}(E_0)$。这对 $D_{\mathbf{q}}(E_0)$ 已经意味着对**无穷部分级数**求和（即"无规相近似"，RPA）：

$$
\varepsilon_{RPA}(\mathbf{q},E) = 1 - 2v(\mathbf{q})\sum_{\mathbf{k}}\frac{\langle n_{\mathbf{k}}\rangle^{(0)} - \langle n_{\mathbf{k}+\mathbf{q}}\rangle^{(0)}}{E + \mathrm{i}0^{+} + \varepsilon(\mathbf{k}) - \varepsilon(\mathbf{k}+\mathbf{q})},
$$

其中已经完成到推迟函数的过渡。

> **注意：** 尽管形式上相似，$\varepsilon_{RPA}(\mathbf{q},E)$ 与 $\varepsilon^{(0)}(\mathbf{q},E)$ 并不是同一个东西。后者的倒数展开
> $$
> \frac{1}{\varepsilon} = \frac{1}{1 - v\Lambda} = \sum_{n=0}^{\infty}(v\Lambda)^n = \underbrace{1 + v\Lambda}_{\text{（零阶近似）}} + \cdots
> $$
> 只对应头两项。而对无穷部分级数求和后，$\varepsilon_{RPA}(\mathbf{q},E)$ 会出现一个物理上意义重大的**额外零点**——它正可识别为**等离激元激发**。

在第四章中，我们已经用格林函数运动方程方法（内容：含 (4.143) 的推导）导出了与 RPA 完全相同的结果，那里的物理解释（特别是图解说明）可以完整搬用，无需重复。

## Hubbard 模型：矩阵解（推导如下）

对 Hubbard 模型，矩阵方程必须以相互作用

$$
\widetilde{V} \equiv \frac{U}{N}\begin{pmatrix}0 & 1\\ 1 & 0\end{pmatrix}
$$

求解，它不允许像 Jellium 模型那样的直接简化。元素 $D_{\mathbf{q}\sigma\sigma'}(E_0)$ 的显式计算如下。

> **（推导如下：教材习题 6.3.2）** —— Hubbard 模型中 $v_{\sigma\sigma'}(\mathbf{q}) \equiv \frac{U}{N}\delta_{\sigma,-\sigma'}$，代入自旋分辨 Dyson 方程得
> $$
> D_{\mathbf{q}\sigma\sigma'}(E_0) = \hbar\Lambda_{\mathbf{q}\sigma\sigma'}(E_0) + \frac{U}{N}\sum_{\sigma''}\Lambda_{\mathbf{q}\sigma\sigma''}(E_0)\,D_{\mathbf{q},-\sigma'',\sigma'}(E_0),
> $$
> 逐元素写出：
> $$
> D_{\uparrow\uparrow} = \hbar\Lambda_{\uparrow\uparrow} + \frac{U}{N}\left(\Lambda_{\uparrow\uparrow}D_{\downarrow\uparrow} + \Lambda_{\uparrow\downarrow}D_{\uparrow\uparrow}\right),
> $$
> $$
> D_{\uparrow\downarrow} = \hbar\Lambda_{\uparrow\downarrow} + \frac{U}{N}\left(\Lambda_{\uparrow\uparrow}D_{\downarrow\downarrow} + \Lambda_{\uparrow\downarrow}D_{\uparrow\downarrow}\right),
> $$
> $$
> D_{\downarrow\uparrow} = \hbar\Lambda_{\downarrow\uparrow} + \frac{U}{N}\left(\Lambda_{\downarrow\uparrow}D_{\downarrow\uparrow} + \Lambda_{\downarrow\downarrow}D_{\uparrow\uparrow}\right),
> $$
> $$
> D_{\downarrow\downarrow} = \hbar\Lambda_{\downarrow\downarrow} + \frac{U}{N}\left(\Lambda_{\downarrow\uparrow}D_{\downarrow\downarrow} + \Lambda_{\downarrow\downarrow}D_{\uparrow\downarrow}\right),
> $$
> 其中能量与波矢自变量 $\mathbf{q},E_0$ 均略去。由第一与第三式解出耦联方程组：
> $$
> D_{\uparrow\uparrow} = \frac{\hbar\Lambda_{\uparrow\uparrow}}{1 - \frac{U}{N}\Lambda_{\uparrow\downarrow}} + \frac{\frac{U}{N}\Lambda_{\uparrow\uparrow}}{1 - \frac{U}{N}\Lambda_{\uparrow\downarrow}}\,D_{\downarrow\uparrow},\qquad
> D_{\downarrow\uparrow} = \frac{\hbar\Lambda_{\downarrow\uparrow}}{1 - \frac{U}{N}\Lambda_{\downarrow\uparrow}} + \frac{\frac{U}{N}\Lambda_{\downarrow\downarrow}}{1 - \frac{U}{N}\Lambda_{\downarrow\uparrow}}\,D_{\uparrow\uparrow},
> $$
> 代入并化简，最终得到
> $$
> \boxed{\;D_{\mathbf{q}\uparrow\uparrow}(E_0) = \hbar\frac{\Lambda_{\mathbf{q}\uparrow\uparrow}(E_0)}{1 - \frac{U}{N}\left(\Lambda_{\mathbf{q}\uparrow\downarrow}(E_0) + \Lambda_{\mathbf{q}\downarrow\uparrow}(E_0)\right) + \frac{U^2}{N^2}\left(\Lambda_{\mathbf{q}\downarrow\uparrow}(E_0)\Lambda_{\mathbf{q}\uparrow\downarrow}(E_0) - \Lambda_{\mathbf{q}\uparrow\uparrow}(E_0)\Lambda_{\mathbf{q}\downarrow\downarrow}(E_0)\right)}\;}
> $$
> 由此 $D_{\mathbf{q}\downarrow\uparrow}$ 也随之确定；另外两个矩阵元 $D_{\mathbf{q}\downarrow\downarrow}$ 与 $D_{\mathbf{q}\uparrow\downarrow}$ 只需在方程两侧同时做自旋翻转（$\uparrow\leftrightarrow\downarrow$）即得。
>
> **（参数情形，顺磁性系统）** 对参数磁性电子系统可取
> $$
> \Lambda_{\mathbf{q}\sigma\sigma'}(E_0) = \Lambda_{\mathbf{q},-\sigma,-\sigma'}(E_0),
> $$
> 于是
> $$
> \Lambda_{\mathbf{q}\uparrow\uparrow} + \Lambda_{\mathbf{q}\downarrow\uparrow} \equiv \Lambda_{\mathbf{q}\downarrow\downarrow} + \Lambda_{\mathbf{q}\uparrow\downarrow} \equiv \frac{1}{2}\Lambda_{\mathbf{q}}(E_0),
> $$
> 密度关联化为与 Jellium 模型结构相似、大幅简化的形式：
> $$
> \boxed{\;D_{\mathbf{q}}(E_0) = \frac{\hbar\,\Lambda_{\mathbf{q}}(E_0)}{1 - \frac{U}{2N}\Lambda_{\mathbf{q}}(E_0)}\;}
> $$
>
> **（零阶极化传播子，RPA）** 取极化传播子零阶结果 $\Lambda_{\mathbf{q}}^{(0)}(E_0)$（即无相互作用、顺磁系统的结果），就得到 RPA 密度关联：
> $$
> D_{\mathbf{q}}(E_0) = \frac{\hbar\,\Lambda_{\mathbf{q}}^{(0)}(E_0)}{1 - \frac{U}{2N}\Lambda_{\mathbf{q}}^{(0)}(E_0)}.
> $$
>
> **（饱和铁磁系统）** 系统只含 $\uparrow$ 电子时，Hubbard 模型中不存在相互作用（没有自旋相反的电子对），因此
> $$
> D_{\mathbf{q}}(E_0) \equiv D_{\mathbf{q}\uparrow\uparrow}(E_0) = \hbar\Lambda_{\mathbf{q}}(E_0) = \hbar\Lambda_{\mathbf{q}\uparrow\uparrow}(E_0) = \hbar\Lambda_{\mathbf{q}}^{(0)}(E_0).
> $$

---

# 主题三：有效相互作用与双重重正化

## 裸相互作用与有效相互作用

独立于描述密度关联这一原始目标，极化传播子还有重要的应用。在 $T=0$ 的 §5.6.2 中，我们示范了如何借助极化传播子发展**有效相互作用**的概念；这几乎可以直接搬进 $T\neq 0$ 的 Matsubara 形式体系。

单粒子松原函数的自能图（Part II）中，有些在一条相互作用线里含有可约或不可约的极化贡献（例子见图 6.38）。引入"有效相互作用" $v_{\text{eff},\sigma\sigma'}$ 可以概括全部这类图。图式上，我们把"裸"相互作用与"有效"相互作用区分如下：

![Fig 6.38–6.39：带极化贡献的自能图（6.38）与有效相互作用图的一般结构（6.39）](./fig有限温微扰论3两粒子松原函数/fig638fig639.png)

![裸相互作用与有效相互作用的记号对应](./fig有限温微扰论3两粒子松原函数/ExprEffIntAndBareInt.png)


一般图的结构如图 6.39 所示：由两条"裸"相互作用线与某个可约或不可约（自旋分辨）极化贡献构成。所有这些图之和就给出有效相互作用。合理地为有效相互作用也定义"阶"：$v_{\text{eff},\sigma\sigma'}$ 的 **$n$ 阶包含 $(n+1)$ 条相互作用线**。零、一、二、三阶的例子见图 6.40。

![Fig 6.40–6.41：有效相互作用 $v_{\mathrm{eff},\sigma\sigma'}(\mathbf{q},E)$ 的显式图贡献（6.40，未计自旋）与其示意图结构（6.41）](./fig有限温微扰论3两粒子松原函数/fig640fig641.png)

除零阶外，每个图都具有图 6.41 所示的结构：部分 (a) 为"裸"相互作用；部分 (b) 为"任意"不可约极化贡献；部分 (c) 为 $v_{\text{eff},\sigma\sigma'}(\mathbf{q},E)$ 的"任意"图。显然，若在 (b) 中遍历所有不可约极化贡献（即极化传播子）、在 (c) 中遍历有效相互作用的全部图，并加上零阶图，就得到所有相互作用图。这又可以写成 Dyson 方程：

$$
-\frac{1}{\hbar}\,v_{\text{eff},\sigma\sigma'}(\mathbf{q},E) = -\frac{1}{\hbar}\,v_{\sigma\sigma'}(\mathbf{q}) + \sum_{\sigma''\sigma'''}\left(-\frac{1}{\hbar}\,v_{\sigma\sigma''}(\mathbf{q})\right)\left(-\hbar\Lambda_{\mathbf{q}\sigma''\sigma'''}(E)\right)\left(-\frac{1}{\hbar}\,v_{\text{eff},\sigma'''\sigma'}(\mathbf{q},E)\right),
$$

矩阵表述更易理解：

$$
\boxed{\;\widetilde{v}_{\text{eff}}(\mathbf{q},E) = \widetilde{V}(\mathbf{q}) + \widetilde{V}(\mathbf{q})\,\widetilde{\Lambda}_{\mathbf{q}}(E)\,\widetilde{v}_{\text{eff}}(\mathbf{q},E) = \frac{\widetilde{V}(\mathbf{q})}{\mathbb{1} - \widetilde{V}(\mathbf{q})\,\widetilde{\Lambda}_{\mathbf{q}}(E)}\;}
$$

这样，"有效"相互作用就完全由极化传播子确定。

## Jellium 模型：动态屏蔽（推导如下）

与密度关联的情形类似，Jellium 模型"裸"相互作用矩阵的特殊形式允许用实际的（非自旋分辨的）极化传播子来表达有效相互作用，即 Jellium 模型的有效相互作用不显含自旋依赖：

$$
\boxed{\;v_{\text{eff},\sigma\sigma'}(\mathbf{q},E) \equiv v_{\text{eff}}(\mathbf{q},E) = \frac{v(\mathbf{q})}{1 - v(\mathbf{q})\Lambda_{\mathbf{q}}(E)} = \frac{v(\mathbf{q})}{\varepsilon(\mathbf{q},E)}\;}
$$

> **（推导如下：教材习题 6.3.3）** —— 两种证法。
>
> 方法一（在 Dyson 方程中使用特殊相互作用形式）：Jellium 模型中 $v_{\sigma\sigma''}(\mathbf{q}) \equiv v(\mathbf{q})$ 与自旋无关，于是
> $$
> v_{\text{eff},\sigma\sigma'}(\mathbf{q},E) = v(\mathbf{q}) + v(\mathbf{q})\sum_{\sigma''\sigma'''}\Lambda_{\mathbf{q}\sigma''\sigma'''}(E)\,v_{\text{eff},\sigma'''\sigma'}(\mathbf{q},E).
> $$
> 右侧不依赖 $\sigma$，因此有效相互作用至少不依赖其第一个自旋指标；可把右侧的有效相互作用作为公共因子提出：
> $$
> v_{\text{eff},\sigma\sigma'}(\mathbf{q},E_0) = v(\mathbf{q}) + v(\mathbf{q})\,v_{\text{eff},\sigma\sigma'}(\mathbf{q},E)\,\Lambda_{\mathbf{q}}(E),
> $$
> 即 Jellium 模型中有效相互作用与相互作用伙伴的自旋无关：
> $$
> v_{\text{eff},\sigma\sigma'}(\mathbf{q},E_0) \equiv v_{\text{eff}}(\mathbf{q},E_0) = \frac{v(\mathbf{q})}{1 - v(\mathbf{q})\Lambda_{\mathbf{q}}(E_0)}.
> $$
>
> 方法二（矩阵方程直接相乘）：由于 $\widetilde{V}(\mathbf{q})\widetilde{\Lambda}_{\mathbf{q}}(E_0)$ 的两行相同，矩阵 Dyson 方程的前两行（$v_{\text{eff},\uparrow\uparrow}$ 与 $v_{\text{eff},\downarrow\uparrow}$ 的方程）完全相同，逐元素解之得
> $$
> v_{\text{eff},\uparrow\uparrow}(\mathbf{q},E_0) = v(\mathbf{q}) + v(\mathbf{q})\Big(\Lambda_{\mathbf{q}\uparrow\uparrow} + \Lambda_{\mathbf{q}\downarrow\uparrow} + \Lambda_{\mathbf{q}\uparrow\downarrow} + \Lambda_{\mathbf{q}\downarrow\downarrow}\Big)v_{\text{eff},\uparrow\uparrow}(\mathbf{q},E_0) = v(\mathbf{q}) + v(\mathbf{q})\Lambda_{\mathbf{q}}(E_0)v_{\text{eff},\uparrow\uparrow}(\mathbf{q},E_0),
> $$
> $$
> \Rightarrow v_{\text{eff},\uparrow\uparrow}(\mathbf{q},E_0) = \frac{v(\mathbf{q})}{1 - v(\mathbf{q})\Lambda_{\mathbf{q}}(E_0)} = v_{\text{eff},\downarrow\uparrow}(\mathbf{q},E_0),
> $$
> 其余两个矩阵元完全类似，再次验证 Jellium 模型的有效相互作用与所涉及电子的自旋无关。

最后一步引入了介电函数 $\varepsilon(\mathbf{q},E) = 1 - v(\mathbf{q})\Lambda_{\mathbf{q}}(E)$。显然，$\varepsilon(\mathbf{q},E)$ 描述的是"裸"相互作用被关联的 Jellium 粒子系统的极化所**动态屏蔽**。在 $T=0$ 情形下，我们曾通过因果格林函数得到过类似结果（(5.189) 的内容）。

## Hubbard 模型：自旋依赖的有效相互作用（推导如下）

Hubbard 模型的特殊相互作用矩阵不允许像 Jellium 模型那样的简化，必须直接使用矩阵方程。注意：与"裸"相互作用相反，有效相互作用的上、下顶点**不一定**需要携带不同自旋，即一般情形下 Hubbard 模型中也有 $v_{\text{eff},\sigma\sigma}(\mathbf{q},E) \neq 0$。

> **（推导如下：教材习题 6.3.4）** —— Hubbard 模型中 $v_{\sigma\sigma''}(\mathbf{q}) \equiv \frac{U}{N}\delta_{\sigma,-\sigma''}$，代入 Dyson 方程：
> $$
> v_{\text{eff},\sigma\sigma'}(\mathbf{q},E_0) = \frac{U}{N}\delta_{\sigma,-\sigma'} + \frac{U}{N}\sum_{\sigma'''}\Lambda_{\mathbf{q},-\sigma,\sigma'''}(E_0)\,v_{\text{eff},\sigma'''\sigma'}(\mathbf{q},E_0),
> $$
> 即
> $$
> v_{\text{eff},\uparrow\uparrow} = \frac{U}{N}\left(\Lambda_{\downarrow\uparrow}v_{\text{eff},\uparrow\uparrow} + \Lambda_{\downarrow\downarrow}v_{\text{eff},\downarrow\uparrow}\right),\qquad
> v_{\text{eff},\downarrow\uparrow} = \frac{U}{N} + \frac{U}{N}\left(\Lambda_{\uparrow\uparrow}v_{\text{eff},\uparrow\uparrow} + \Lambda_{\uparrow\downarrow}v_{\text{eff},\downarrow\uparrow}\right),
> $$
> 其中自变量 $(\mathbf{q},E_0)$ 均略去。由第一式解出 $v_{\text{eff},\uparrow\uparrow} = \frac{U}{N}\Lambda_{\downarrow\downarrow}v_{\text{eff},\downarrow\uparrow}/\left(1 - \frac{U}{N}\Lambda_{\downarrow\uparrow}\right)$，代入第二式并化简，得到
> $$
> \boxed{\;v_{\text{eff},\uparrow\uparrow}(\mathbf{q},E_0) = \frac{\frac{U^2}{N^2}\Lambda_{\mathbf{q}\downarrow\downarrow}(E_0)}{\left(1 - \frac{U}{N}\Lambda_{\mathbf{q}\downarrow\uparrow}(E_0)\right)\left(1 - \frac{U}{N}\Lambda_{\mathbf{q}\uparrow\downarrow}(E_0)\right) - \frac{U^2}{N^2}\Lambda_{\mathbf{q}\downarrow\downarrow}(E_0)\Lambda_{\mathbf{q}\uparrow\uparrow}(E_0)}\;}
> $$
> 再代回非对角元素的表达式：
> $$
> \boxed{\;v_{\text{eff},\downarrow\uparrow}(\mathbf{q},E_0) = \frac{\frac{U}{N}\left(1 - \frac{U}{N}\Lambda_{\mathbf{q}\downarrow\uparrow}(E_0)\right)}{\left(1 - \frac{U}{N}\Lambda_{\mathbf{q}\downarrow\uparrow}(E_0)\right)\left(1 - \frac{U}{N}\Lambda_{\mathbf{q}\uparrow\downarrow}(E_0)\right) - \frac{U^2}{N^2}\Lambda_{\mathbf{q}\downarrow\downarrow}(E_0)\Lambda_{\mathbf{q}\uparrow\uparrow}(E_0)}\;}
> $$
> 其余两个元素只需做自旋交换（$\uparrow\leftrightarrow\downarrow$）。Hubbard 模型的有效相互作用因而具有显式的自旋依赖；这里与"裸"相互作用不同，对自旋相同的相互作用伙伴一般也非零（$v_{\text{eff},\sigma\sigma} \neq 0$）。
>
> **（参数情形，顺磁性系统）** 取 $\Lambda_{\mathbf{q}\sigma\sigma'}(E_0) \equiv \Lambda_{\mathbf{q},-\sigma,-\sigma'}(E_0)$，引入缩写
> $$
> \Lambda_{\mathbf{q}}^{(+)}(E_0) = \Lambda_{\mathbf{q}\sigma\sigma}(E_0),\qquad \Lambda_{\mathbf{q}}^{(-)}(E_0) = \Lambda_{\mathbf{q}\sigma,-\sigma}(E_0),
> $$
> 有效相互作用可写成
> $$
> v_{\text{eff},\sigma\sigma}(\mathbf{q},E_0) = \frac{\frac{U^2}{N^2}\Lambda_{\mathbf{q}}^{(+)}(E_0)}{\left(1 - \frac{U}{N}\Lambda_{\mathbf{q}}^{(-)}(E_0)\right)^2 - \frac{U^2}{N^2}\left(\Lambda_{\mathbf{q}}^{(+)}(E_0)\right)^2},
> $$
> $$
> v_{\text{eff},\sigma,-\sigma}(\mathbf{q},E_0) = \frac{\frac{U}{N}\left(1 - \frac{U}{N}\Lambda_{\mathbf{q}}^{(-)}(E_0)\right)}{\left(1 - \frac{U}{N}\Lambda_{\mathbf{q}}^{(-)}(E_0)\right)^2 - \frac{U^2}{N^2}\left(\Lambda_{\mathbf{q}}^{(+)}(E_0)\right)^2}.
> $$
>
> **（零阶极化传播子，RPA）** 特殊情形 $\Lambda_{\mathbf{q}}(E_0) \to \Lambda_{\mathbf{q}}^{(0)}(E_0)$，即 $\Lambda_{\mathbf{q}\sigma\sigma'}(E_0) \to \frac{1}{4}\Lambda_{\mathbf{q}}^{(0)}(E_0)$，给出"无规相近似"：
> $$
> v_{\text{eff},\sigma\sigma}^{RPA}(\mathbf{q},E_0) = \frac{\left(\frac{U}{2N}\right)^2\Lambda_{\mathbf{q}}^{(0)}(E_0)}{1 - \frac{U}{2N}\Lambda_{\mathbf{q}}^{(0)}(E_0)},\qquad
> v_{\text{eff},\sigma,-\sigma}^{RPA}(\mathbf{q},E_0) = \frac{\frac{U}{N}\left(1 - \frac{U}{4N}\Lambda_{\mathbf{q}}^{(0)}(E_0)\right)}{1 - \frac{U}{2N}\Lambda_{\mathbf{q}}^{(0)}(E_0)}.
> $$

## 骨架图、双重重正化与双重计数

有效相互作用为（近似）确定单粒子松原函数或相应自能提供了另一条途径，与 Part II 的方案互为替代：

- 在单粒子松原函数的自能图中，**抑制**一切在至少一条相互作用线中含有极化贡献的图；在剩余的**骨架图**中，把"裸"相互作用替换为有效相互作用。这再次自动地求和了无穷部分级数，给出自能求值的新概念！

- 还可以像 Part II 那样，把剩余骨架图中的自由传播子替换为完整"敷饰"传播子，后者由 Dyson 方程（其图解形式即 Part II 的图 6.9）**自洽**地确定（图 6.42）。

![Fig 6.42–6.43：自能图的"双重"重正化（有效相互作用 + 完整传播子，6.42）与双重计数的例子（6.43）](./fig有限温微扰论3两粒子松原函数/fig642fig643.png)

由于各种重正化，自然存在**重复计数某些图**的危险，这会严重破坏结果，必须不惜一切代价避免。例如，图 6.42 中第一个图的相互作用线**不得**被"敷饰"——否则图 6.43 的图就会出现两次：一次经由完整传播子的一阶自能（Part II 图 6.10 的左图），一次经由一阶有效相互作用（图 6.40 的第二个图）。出于同样原因，Part II 中"直接项"（其图见 Part II 图 6.19）也不以有效相互作用的形式出现在图 6.42 的展开中：可以立即确认，相应的图已经完全包含在图 6.42 的第二求和项之内。

---

# 主题四：顶角函数与梯近似

## 顶角贡献与顶角函数

我们还要讨论另一种通过引入"图块"简化图级数的方案，用极化传播子（图 6.36）为例说明。假设与前面各节相同（如 $v_{\sigma_k\sigma_l}(\mathbf{q})$ 的模型形式）。低阶极化贡献图的例子见图 6.34，与 $T=0$ 的 §5.6.3 图自然拓扑等价。

> **定义（顶角贡献）**：顶角贡献 = （自旋分辨）极化贡献中具有两条粒子线连接与一条相互作用线连接的图部分。

例子见图 6.44。这些图型也出现在 §5.6.3 中，只是那里传播子为 $T=0$ 因果格林函数、相互作用线可特别归诸 Jellium 模型。注意：**零阶顶角贡献仅由一个（顶点）点构成**，但它也必须计入。

![Fig 6.44：顶角贡献图的例子](./fig有限温微扰论3两粒子松原函数/fig644.png)

> **定义（不可约顶角贡献）**：不可约顶角贡献 = 不能通过分离一条传播子而分出一个完整的自能图、也不能通过切断一条相互作用线而分出一个完整极化贡献图的顶角贡献。

图 6.44 中第一、二、四、五个图显然不可约，第三、六个图可约：从第三个图分离一条传播子可得一阶自能图，从第六个图切断一条相互作用线可得零阶极化图。

> **定义（顶角函数）**：顶角函数 $\Gamma_{\sigma\sigma'}(\mathbf{q}E_0;\mathbf{k}E)$ = 所有不可约顶角贡献之和。

顶角函数使用图 6.45 的记号：

![Fig 6.45–6.46：顶角函数的记号（6.45）与借助顶角函数表示极化传播子（6.46）](./fig有限温微扰论3两粒子松原函数/fig645fig646.png)

它应这样解读：由于在右侧顶点处总须接上一条相互作用线，从顶角函数必然伸出两条传播子；$\mathbf{q}$ 与 $E_0$ 是这两条传播子之间转移的波数与能量。回顾本节的一般假设：这两条传播子具有相同的自旋量子数 $\sigma$。记号 $(\mathbf{k}E\sigma')$ 表示左上顶点处的（外）传播子，顶角函数在左侧与之相连；连入左下顶点的传播子则具有量子数 $(\mathbf{k}+\mathbf{q}, E+E_0, \sigma')$。

顶角函数的零阶贡献就是孤立的顶点，以 $\delta_{\sigma\sigma'}$ 进入 $\Gamma_{\sigma\sigma'}$。图 6.44 中第一、二、四、五个图是一、二阶的不可约顶角贡献；其余贡献的结构与 §5.6.3 中定义之后的例子相同。

## 极化传播子的顶角表示

最初的目标——用顶角函数表示极化传播子——可按图 6.46 实现，求值如下：

$$
\boxed{\;-\hbar\Lambda_{\mathbf{q}\sigma\sigma'}(E_0) = \frac{\varepsilon}{\hbar\beta}\sum_{\mathbf{k}E}\left(-G_{\mathbf{k}\sigma}^{M}(E)\right)\left(-G_{\mathbf{k}+\mathbf{q}\sigma}^{M}(E+E_0)\right)\Gamma_{\sigma'\sigma}(\mathbf{q}E_0;\mathbf{k}E)\;}
$$

符号因子 $\varepsilon$ 来自环规则，因子 $1/\hbar\beta$ 来自两条"完整"传播子作为"外接"（与 $-D_{\mathbf{q}\sigma\sigma'}(E_0)$ 中 $1/\hbar\beta$ 因子的来源相同）。注意：这里出现的是**完整**单粒子松原函数 $G^M$。

这个极化传播子的表示原则上仍无法严格求值，我们讨论两个可能的近似：

- 若取顶角函数的零级近似 $\Gamma_{\sigma'\sigma}(\mathbf{q}E_0;\mathbf{k}E) \to \delta_{\sigma'\sigma}$，并把"完整"传播子取为"自由"传播子，则立即回到已经讨论过的简单结果（零阶密度关联 $\Lambda_{\mathbf{q}}^{(0)}$ 的两个等价表达式）。

- 更进一步的近似是**梯近似**（ladder approximation），它可以视为顶角函数的（近似）Dyson 方程，见图 6.47。为与极化传播子精确表示（$\Gamma$ 的表示式）直接衔接，这里出现的传播子理解为完整单粒子松原函数。此外，只限于具有 (Jellium 型) 特殊相互作用矩阵的体系——容易看出，对 Hubbard 模型图 6.47 中所有项都为零，梯近似在这种情形下失去意义（对 Hubbard 模型的替代方案见下一主题）。求值给出

![Fig 6.47：顶角函数的梯近似](./fig有限温微扰论3两粒子松原函数/fig647.png)

$$
\Gamma_{\sigma\sigma'}^{LA}(\mathbf{q}E_0;\mathbf{k}E) = \delta_{\sigma\sigma'} + \delta_{\sigma\sigma'}\left(-\frac{1}{\hbar}\right)\frac{1}{\hbar\beta}\sum_{\mathbf{p}E_1}v(\mathbf{k}-\mathbf{p})\left(-G_{\mathbf{p}\sigma}^{M}(E_1)\right)\left(-G_{\mathbf{p}+\mathbf{q}\sigma}^{M}(E_1+E_0)\right)\Gamma_{\sigma\sigma}^{LA}(\mathbf{q}E_0;\mathbf{p}E_1),
$$

每个顶点处假定的自旋守恒保证了在此特定近似中只有 $\sigma = \sigma'$ 有非零贡献。

## 常数相互作用的极限

若可假设 $v(\mathbf{k}-\mathbf{p})$ 对波数依赖很弱甚至完全不依赖，从而可用

$$
v(\mathbf{k}) \rightarrow v_{0} = \frac{1}{N}\sum_{\mathbf{k}}v(\mathbf{k})
$$

近似，则有

$$
\Gamma_{\sigma\sigma'}^{LA}(\mathbf{q}E_0;\mathbf{k}E) = \delta_{\sigma\sigma'} - \delta_{\sigma\sigma'}\frac{v_{0}}{\hbar^{2}\beta}\sum_{\mathbf{p}E_1}G_{\mathbf{p}\sigma}^{M}(E_1)\,G_{\mathbf{p}+\mathbf{q}\sigma}^{M}(E_1+E_0)\,\Gamma_{\sigma\sigma}^{LA}(\mathbf{q}E_0;\mathbf{p}E_1),
$$

右侧不再依赖 $(\mathbf{k},E)$，顶角函数因此简化成

$$
\Gamma_{\sigma\sigma'}^{LA}(\mathbf{q}E_0;\mathbf{k}E) \to \delta_{\sigma\sigma'}\,\Gamma_{\sigma}^{LA}(\mathbf{q}E_0).
$$

引入缩写

$$
-\hbar\widehat{\Lambda}_{\mathbf{q}\sigma}(E_0) = \frac{\varepsilon}{\hbar\beta}\sum_{\mathbf{p}E_1}G_{\mathbf{p}\sigma}^{M}(E_1)\,G_{\mathbf{p}+\mathbf{q}\sigma}^{M}(E_1+E_0),
$$

是有帮助的。若把"完整"传播子替换为"自由"传播子，则除因子 $1/2$ 外，恰好得到零阶密度关联的 $\hbar\Lambda_{\mathbf{q}}^{(0)}(E_0)$。梯近似于是化简为

$$
\Gamma_{\sigma}^{LA}(\mathbf{q}E_0) = 1 + \varepsilon v_{0}\,\widehat{\Lambda}_{\mathbf{q}\sigma}(E_0)\,\Gamma_{\sigma}^{LA}(\mathbf{q}E_0),
$$

顶角函数的解为

$$
\boxed{\;\Gamma_{\sigma}^{LA}(\mathbf{q}E_0) = \frac{1}{1 - \varepsilon v_{0}\widehat{\Lambda}_{\mathbf{q}\sigma}(E_0)}\;}
$$

把这个结果代入极化传播子的顶角表示，得到相应近似下的自旋分辨极化传播子：

$$
\Lambda_{\mathbf{q}\sigma\sigma'}^{LA}(E_0) = \widehat{\Lambda}_{\mathbf{q}\sigma}(E_0)\,\delta_{\sigma\sigma'}\,\Gamma_{\sigma}^{LA}(\mathbf{q}E_0) = \delta_{\sigma\sigma'}\frac{\widehat{\Lambda}_{\mathbf{q}\sigma}(E_0)}{1 - \varepsilon v_{0}\widehat{\Lambda}_{\mathbf{q}\sigma}(E_0)}.
$$

于是密度关联解中的矩阵 $\widetilde{\Lambda}_{\mathbf{q}}(E_0)$ 现在只有对角元非零。对具有 (Jellium 型) 特殊相互作用矩阵的体系，还可利用总极化传播子 $\Lambda_{\mathbf{q}}(E_0) = \sum_{\sigma\sigma'}\Lambda_{\mathbf{q}\sigma\sigma'}(E_0)$：

$$
\Lambda_{\mathbf{q}}^{LA}(E_0) = \sum_{\sigma\sigma'}\Lambda_{\mathbf{q}\sigma\sigma'}^{LA}(E_0) = \sum_{\sigma}\frac{\widehat{\Lambda}_{\mathbf{q}\sigma}(E_0)}{1 - \varepsilon v_{0}\widehat{\Lambda}_{\mathbf{q}\sigma}(E_0)},
$$

再由密度关联与介电函数的关系式，即近似确定了密度关联与介电函数。

## Λ̂<sub>qσ</sub> 的谱表示与能域求和（推导如下）

还须求值 $\widehat{\Lambda}_{\mathbf{q}\sigma}(E_0)$，至少完成对能量的求和。计算路径与 Part II 中 $I_{nml}(E)$ 的路径相同：用单粒子松原函数的谱表示 + 松原频率求和公式。

> **（推导如下：教材习题 6.3.5）** —— 对
> $$
> -\hbar\widehat{\Lambda}_{\mathbf{q}\sigma}(E_0) = \frac{\varepsilon}{\hbar\beta}\sum_{\mathbf{p}E_1}G_{\mathbf{p}\sigma}^{M}(E_1)\,G_{\mathbf{p}+\mathbf{q}\sigma}^{M}(E_1+E_0)
> $$
> 做能量求和。利用单粒子松原函数的**谱表示**（Part I，内容转写）：
> $$
> G_{\mathbf{p}\sigma}^{M}(E_1) = \int_{-\infty}^{+\infty}dE'\,\frac{S_{\mathbf{p}\sigma}(E')}{\mathrm{i}E_1 - E'},
> $$
> 得
> $$
> \widehat{\Lambda}_{\mathbf{q}\sigma}(E_0) = -\frac{\varepsilon}{\hbar^{2}\beta}\sum_{\mathbf{p}}\int_{-\infty}^{+\infty}dx\int_{-\infty}^{+\infty}dy\,S_{\mathbf{p}\sigma}(x)S_{\mathbf{p}+\mathbf{q}\sigma}(y)\,F_{E_0}(x,y),
> $$
> 其中
> $$
> F_{E_0}(x,y) = \sum_{E_1}\frac{1}{\mathrm{i}E_1 - x}\cdot\frac{1}{\mathrm{i}E_1 + \mathrm{i}E_0 - y} = \sum_{E_1}H_{x,y}(\mathrm{i}E_1).
> $$
> 对松原能量求和可借助（Part I 的）围道积分公式，或等价地使用（Part II 中证明的）残数公式
> $$
> \sum_{E_1}H_{x,y}(\mathrm{i}E_1) = -\varepsilon\beta\sum_{\widehat{E}_i}f_{\varepsilon}(\widehat{E}_i)\,\mathrm{Res}_{\widehat{E}_i}H(E),
> $$
> $H(E)$ 有两个极点：
> $$
> \mathrm{i}\widehat{E}_1 = x,\quad \mathrm{Res}_{\widehat{E}_1} = (x + \mathrm{i}E_0 - y)^{-1};\qquad
> \mathrm{i}\widehat{E}_2 = y - \mathrm{i}E_0,\quad \mathrm{Res}_{\widehat{E}_2} = (y - \mathrm{i}E_0 - x)^{-1}.
> $$
> 于是
> $$
> F_{E_0}(x,y) = -\varepsilon\beta\left(\frac{f_{\varepsilon}(x)}{x + \mathrm{i}E_0 - y} + \frac{f_{\varepsilon}(y - \mathrm{i}E_0)}{y - \mathrm{i}E_0 - x}\right).
> $$
> $E_0$ 是玻色型松原能量，$e^{-\mathrm{i}\beta E_0} = +1$，故 $f_{\varepsilon}(y - \mathrm{i}E_0) = f_{\varepsilon}(y)$，两项合并为
> $$
> F_{E_0}(x,y) = -\varepsilon\beta\frac{f_{\varepsilon}(x) - f_{\varepsilon}(y)}{\mathrm{i}E_0 + x - y},
> $$
> 代回即得
> $$
> \boxed{\;\widehat{\Lambda}_{\mathbf{q}\sigma}(E_0) = \frac{1}{\hbar^{2}}\sum_{\mathbf{p}}\int_{-\infty}^{+\infty}dx\int_{-\infty}^{+\infty}dy\,\frac{S_{\mathbf{p}\sigma}(x)\,S_{\mathbf{p}+\mathbf{q}\sigma}(y)}{\mathrm{i}E_0 + x - y}\left(f_{\varepsilon}(x) - f_{\varepsilon}(y)\right)\;}
> $$
> 其中 $f_{\varepsilon}(x) = 1/(e^{\beta x}-\varepsilon)$ 为广义 Fermi-Dirac/Bose-Einstein 函数。

被积函数中的谱密度最终还必须用单粒子松原函数的 Dyson 方程"以某种方式"确定。替换 $S_{\mathbf{p}\sigma} \to S_{\mathbf{p}}^{(0)}$ 则再次给出零阶结果 $\Lambda_{\mathbf{q}}^{(0)}(E_0)$（即自由粒子-空穴气泡）。

---

# 主题五：横向自旋磁化率

## 定义与运动方程方法回顾

最后讨论顶角函数与梯近似的一个特殊应用：**横向自旋磁化率**。它是确定相互作用电子系统磁性（自旋波、磁振子）的重要量，例如可在 Hubbard 模型框架中描述。自旋磁化率是两粒子松原函数，与前面详细讨论的密度关联有相似之处。它定义为

$$
\chi_{ij}^{+-}(E_0) = -\frac{\gamma}{\hbar^{2}}\left\langle\!\left\langle \sigma_{i}^{+};\sigma_{j}^{-}\right\rangle\!\right\rangle_{E_0}^{M},\qquad \gamma = \frac{\mu_{0}}{V\hbar}g^{2}\mu_{B}^{2},
$$

$\gamma$ 是某个量纲常数，对我们这里的目标不重要。自旋算符作用于巡游电子，可用费米子算符表示（Wannier 算符，$i$ 为格点指标）：

$$
\sigma_{i}^{+} = \hbar a_{i\uparrow}^{+}a_{i\downarrow};\qquad \sigma_{i}^{-} = \hbar a_{i\downarrow}^{+}a_{i\uparrow};\qquad \sigma_{i}^{z} = \frac{\hbar}{2}\left(a_{i\uparrow}^{+}a_{i\uparrow} - a_{i\downarrow}^{+}a_{i\downarrow}\right).
$$

变换到波数：

$$
\chi_{\mathbf{q}}^{+-}(E_0) = \frac{1}{N}\sum_{i,j}\chi_{ij}^{+-}(E_0)\,e^{\mathrm{i}\mathbf{q}\cdot(\mathbf{R}_i-\mathbf{R}_j)} \equiv -\frac{\gamma}{N}\,\widehat{\chi}_{\mathbf{q}}(E_0),
$$

把 $\widehat{\chi}_{\mathbf{q}}(E_0)$ 解释为实际的"自旋磁化率"：

$$
\widehat{\chi}_{\mathbf{q}}(E_0) = \sum_{\mathbf{p},\mathbf{k}}\left\langle\!\left\langle a_{\mathbf{k}\uparrow}^{+}a_{\mathbf{k}+\mathbf{q}\downarrow};\,a_{\mathbf{p}\downarrow}^{+}a_{\mathbf{p}-\mathbf{q}\uparrow}\right\rangle\!\right\rangle_{E_0}^{M}.
$$

![Fig 6.48–6.49：自旋磁化率 $-\widehat{\chi}_{\mathbf{q}}(E_0)$ 的示意图（6.48）与专门用于自旋磁化率的极化传播子（6.49）](./fig有限温微扰论3两粒子松原函数/fig648fig649.png)

除自旋指标外，这个两粒子松原函数与（自旋分辨）密度关联完全对应。差别在于：现在每条外接上耦合的两条传播子**成对携带不同的自旋**。图 6.48 取代图 6.33；自旋指标固定，只对 $\mathbf{k}$ 与 $\mathbf{p}$ 求和。

## 自旋磁化率 = 极化传播子

从自旋分辨密度关联的 Dyson 方程（图 6.37 的内容）可以看出：对自旋磁化率，方程第二行的求和没有非零贡献——由于顶点处自旋守恒与特殊的外接结构，$-\widehat{\chi}_{\mathbf{q}}(E_0)$ 不存在可约极化贡献（图 6.49）：

$$
-\widehat{\chi}_{\mathbf{q}}(E_0) \equiv -\hbar\Lambda_{\mathbf{q}}^{\uparrow\downarrow}(E_0).
$$

**自旋磁化率因此直接对应其极化传播子！** 顶角函数原则上正如上一主题所定义，只需注意自旋的具体枚举（图 6.50）：

![Fig 6.50：专门用于自旋磁化率的顶角函数](./fig有限温微扰论3两粒子松原函数/fig650.png)

$$
-\widehat{\chi}_{\mathbf{q}}(E_0) = \frac{\varepsilon}{\hbar\beta}\sum_{\mathbf{k}E}G_{\mathbf{k}\uparrow}^{M}(E)\,G_{\mathbf{k}+\mathbf{q}\downarrow}^{M}(E+E_0)\,\Gamma_{\uparrow\downarrow}(\mathbf{q}E_0;\mathbf{k}E).
$$

这个表达式仍然是精确的。

## Hubbard 模型的梯近似

在我们将集中讨论的 Hubbard 模型这一特殊情形下，梯近似（上主题中梯近似方程的 Hubbard 版）可以精确地解出：

$$
\Gamma_{\uparrow\downarrow}^{LA}(\mathbf{q}E_0;\mathbf{k}E) = 1 + \left(-\frac{1}{\hbar}\right)\frac{1}{\hbar\beta}\frac{U}{N}\sum_{\mathbf{p},E_1}G_{\mathbf{p}\uparrow}^{M}(E_1)\,G_{\mathbf{p}+\mathbf{q}\downarrow}^{M}(E_1+E_0)\,\Gamma_{\uparrow\downarrow}^{LA}(\mathbf{q}E_0;\mathbf{p}E_1).
$$

右侧不依赖 $(\mathbf{k},E)$，因此可断言

$$
\Gamma_{\uparrow\downarrow}^{LA}(\mathbf{q}E_0;\mathbf{k}E) \equiv \Gamma_{\uparrow\downarrow}^{LA}(\mathbf{q}E_0),
$$

（当然对右侧方程中的顶角函数也同理）。为解出它，与 $\widehat{\Lambda}_{\mathbf{q}\sigma}(E_0)$ 的定义类似地引入（$\varepsilon = -1$ 是 Hubbard 模型中的值）

$$
-\hbar\widehat{\Lambda}_{\mathbf{q}\uparrow\downarrow}(E_0) = \frac{-1}{\hbar\beta}\sum_{\mathbf{p}E_1}G_{\mathbf{p}\uparrow}^{M}(E_1)\,G_{\mathbf{p}+\mathbf{q}\downarrow}^{M}(E_1+E_0),
$$

这个函数在上一个主题末尾已经做过能量求和；只需把谱表示公式中的自旋指标相应调整：

$$
\widehat{\Lambda}_{\mathbf{q}\uparrow\downarrow}(E_0) = \frac{1}{\hbar^{2}}\sum_{\mathbf{p}}\int_{-\infty}^{+\infty}dx\int_{-\infty}^{+\infty}dy\,\frac{S_{\mathbf{p}\uparrow}(x)\,S_{\mathbf{p}+\mathbf{q}\downarrow}(y)}{\mathrm{i}E_0 + x - y}\left(f_{-}(x) - f_{-}(y)\right),
$$

其中 $f_{-}(x) = 1/(e^{\beta x}+1)$ 为费米分布函数。于是梯近似方程化简为

$$
\Gamma_{\uparrow\downarrow}^{LA}(\mathbf{q}E_0) = 1 - \frac{U}{N}\widehat{\Lambda}_{\mathbf{q}\uparrow\downarrow}(E_0)\,\Gamma_{\uparrow\downarrow}^{LA}(\mathbf{q}E_0) = \frac{1}{1 + \frac{U}{N}\widehat{\Lambda}_{\mathbf{q}\uparrow\downarrow}(E_0)}.
$$

梯近似结合自旋磁化率的顶角表示，给出 Hubbard 模型自旋磁化率的表达式：

$$
\boxed{\;\widehat{\chi}_{\mathbf{q}}^{LA}(E_0) = \Gamma_{\uparrow\downarrow}^{LA}(\mathbf{q}E_0)\left(\hbar\widehat{\Lambda}_{\mathbf{q}\uparrow\downarrow}(E_0)\right) = \frac{\hbar\,\widehat{\Lambda}_{\mathbf{q}\uparrow\downarrow}(E_0)}{1 + \frac{U}{N}\widehat{\Lambda}_{\mathbf{q}\uparrow\downarrow}(E_0)}\;}
$$

> **与运动方程方法结果的比较**：第四章用运动方程方法得到的 Hubbard 模型自旋磁化率（(4.183)，内容转写）为
> $$
> \chi_{\mathbf{q}}^{+-}(E) = \left(\chi_{\mathbf{q}}^{+-}(E)\right)^{(S)} + \chi_{\mathbf{q}}^{+-}(E)\left[\frac{U}{\gamma}\left(\chi_{\mathbf{q}}^{+-}(E)\right)^{(S)}\right],
> $$
> $$
> \chi_{\mathbf{q}}^{+-}(E) = \frac{\left(\chi_{\mathbf{q}}^{+-}(E)\right)^{(S)}}{1 - \gamma^{-1}U\left(\chi_{\mathbf{q}}^{+-}(E)\right)^{(S)}},
> $$
> 其中 $(\chi_{\mathbf{q}}^{+-})^{(S)}$ 为"无自旋翻转过程"的部分。注意到 $\chi_{\mathbf{q}}^{+-} = -(\gamma/N)\widehat{\chi}_{\mathbf{q}}$，两种方法的结果在结构上完全一致：分母中都是"单位元 $-$ 相互作用 $\times$ 极化函数"的形式，这正是 RPA 型求和（无穷粒子-空穴链）的共性。

---

# 速查表

## 密度关联与介电函数

- 密度算符：$\rho_{\mathbf{q}} = \sum_{\mathbf{k}\sigma}a_{\mathbf{k}\sigma}^{+}a_{\mathbf{k}+\mathbf{q}\sigma}$；$\rho_{\mathbf{q}}^{+} \equiv \rho_{-\mathbf{q}}$
- 推迟密度关联与介电函数：$1/\varepsilon(\mathbf{q},E) = 1 + \frac{1}{\hbar}v(\mathbf{q})\langle\!\langle\rho_{\mathbf{q}};\rho_{\mathbf{q}}^{+}\rangle\!\rangle_E^{\text{ret}}$；$v(\mathbf{q}) = e^2/\varepsilon_0Vq^2$（来自 Fourier 变换，普遍成立）
- 诱导电荷：$\rho_{\text{ind}} = (1/\varepsilon - 1)\rho_{\text{ext}}$；$\varepsilon\gg 1$ ⟹ 完全屏蔽；$\varepsilon\to 0$ ⟹ 密度涨落发散 ⟹ 等离激元
- 顶点约定：自旋、动量守恒；矩阵元只依赖转移动量：$v(kl;nm) \to v_{\sigma_k\sigma_l}(\mathbf{k}-\mathbf{n})\delta_{\mathbf{k}+\mathbf{l},\mathbf{m}+\mathbf{n}}\delta_{\sigma_k\sigma_n}\delta_{\sigma_l\sigma_m}$
- Jellium：$v_{\sigma_k\sigma_l}(\mathbf{q}) = v(\mathbf{q})$；Hubbard：$v_{\sigma_k\sigma_l}(\mathbf{q}) = \frac{U}{N}\delta_{\sigma_k,-\sigma_l}$
- 两粒子松原函数：$D_{\mathbf{q}}(E_0) = \langle\!\langle\rho_{\mathbf{q}};\rho_{\mathbf{q}}^{+}\rangle\!\rangle_{E_0}^M$；$D_{\mathbf{q}}(\tau-\tau') = -\langle T_\tau\rho_{\mathbf{q}}(\tau)\rho_{\mathbf{q}}^{+}(\tau')\rangle$；自旋分解 $D_{\mathbf{q}} = \sum_{\sigma\sigma'}D_{\mathbf{q}\sigma\sigma'}$
- 连通图定律 ⟹ 分母抵消，只对连通开放图求和；非连通开放图仅 $\mathbf{q}=0$ 可能（Jellium 中 $v(\mathbf{0})=0$ 无贡献）
- 时间均匀性 ⟹ 外线能量满足 $E_2-E_1 = E_3-E_4 \equiv E_0$，$E_0$ 玻色型（$e^{-\mathrm{i}\beta E_0}=+1$）
- 图规则要点：顶点 $\frac{1}{\hbar\beta}v\delta_{\text{能}}\delta_{\text{动量}}\delta_{\text{自旋}}$；内线 $-G^{0,M} = -\hbar/(\mathrm{i}E_n-\varepsilon(\mathbf{n})+\mu)$；非传播线 $e^{\mathrm{i}E_n0^+/\hbar}$；外接两对 $(-G^{0,M})$；只对 $\mathbf{k},\mathbf{p},E,E'$ 求和；因子 $\frac{1}{\hbar\beta}(-1/\hbar)^n\varepsilon^{S}$（$S$ = 环数；$1/\hbar\beta$ 来自四条外接）
- 零阶结果：$\Lambda_{\mathbf{q}}^{(0)}(E_0) = 2\sum_{\mathbf{k}}\frac{\langle n_{\mathbf{k}}\rangle^{(0)}-\langle n_{\mathbf{k}+\mathbf{q}}\rangle^{(0)}}{\mathrm{i}E_0+\varepsilon(\mathbf{k})-\varepsilon(\mathbf{k}+\mathbf{q})}$
- 零阶介电函数：$1/\varepsilon^{(0)} = 1 + 2v\sum_{\mathbf{k}}\frac{\langle n_{\mathbf{k}}\rangle^{(0)}-\langle n_{\mathbf{k}+\mathbf{q}}\rangle^{(0)}}{E+\mathrm{i}0^++\varepsilon(\mathbf{k})-\varepsilon(\mathbf{k}+\mathbf{q})}$；只有粒子-空穴激发，**无等离激元**

## 极化传播子与 RPA

- 极化贡献 = 带两条相互作用线外接（一进一出）的 $-D_{\mathbf{q}\sigma\sigma'}$ 图；不可约 = 不能切一条相互作用线分解成两个低阶极化图
- 极化传播子 $-\hbar\Lambda_{\mathbf{q}\sigma\sigma'}(E_0)$ = 所有不可约极化贡献之和
- 密度关联 Dyson 方程：$D_{\mathbf{q}\sigma\sigma'} = \hbar\Lambda_{\mathbf{q}\sigma\sigma'} + \sum_{\sigma''\sigma'''}\Lambda_{\mathbf{q}\sigma\sigma''}v_{\sigma''\sigma'''}D_{\mathbf{q}\sigma'''\sigma'}$；矩阵解 $\widetilde{D} = \hbar\widetilde{\Lambda}/(\mathbb{1}-\widetilde{\Lambda}\widetilde{V})$
- Jellium：$\widetilde{V} = v(\mathbf{q})\left(\begin{smallmatrix}1&1\\1&1\end{smallmatrix}\right)$ ⟹ $D_{\mathbf{q}} = \frac{\hbar\Lambda_{\mathbf{q}}}{1-v\Lambda_{\mathbf{q}}}$（习题 6.3.1 两种证明）；$\varepsilon(\mathbf{q},E_0) = 1 - v(\mathbf{q})\Lambda_{\mathbf{q}}(E_0)$
- RPA：$\varepsilon_{RPA} = 1 - 2v\sum_{\mathbf{k}}\frac{\langle n_{\mathbf{k}}\rangle^{(0)}-\langle n_{\mathbf{k}+\mathbf{q}}\rangle^{(0)}}{E+\mathrm{i}0^++\varepsilon(\mathbf{k})-\varepsilon(\mathbf{k}+\mathbf{q})}$；与 $\varepsilon^{(0)}$ 不同（$1/\varepsilon = 1+v\Lambda+\cdots$ 只取头两项），RPA 出现等离激元零点；与 EoM 法 (§4.2.2) 结果相同
- Hubbard：$\widetilde{V} = \frac{U}{N}\left(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right)$；$D_{\mathbf{q}\uparrow\uparrow} = \hbar\Lambda_{\mathbf{q}\uparrow\uparrow}/\left(1 - \frac{U}{N}(\Lambda_{\uparrow\downarrow}+\Lambda_{\downarrow\uparrow}) + \frac{U^2}{N^2}(\Lambda_{\downarrow\uparrow}\Lambda_{\uparrow\downarrow}-\Lambda_{\uparrow\uparrow}\Lambda_{\downarrow\downarrow})\right)$，其余元素自旋翻转（习题 6.3.2）
- Hubbard 参数情形：$D_{\mathbf{q}} = \frac{\hbar\Lambda_{\mathbf{q}}}{1 - \frac{U}{2N}\Lambda_{\mathbf{q}}}$；RPA 取 $\Lambda\to\Lambda^{(0)}$；饱和铁磁：$D_{\mathbf{q}} = \hbar\Lambda_{\mathbf{q}}^{(0)}$

## 有效相互作用

- 记号：$_{\sigma}---_{\sigma'} \Leftrightarrow -\frac{1}{\hbar}v_{\sigma\sigma'}(\mathbf{q})$；有效相互作用 $n$ 阶含 $n+1$ 条相互作用线
- Dyson 方程：$\widetilde{v}_{\text{eff}} = \widetilde{V} + \widetilde{V}\widetilde{\Lambda}\widetilde{v}_{\text{eff}} = \widetilde{V}/(\mathbb{1}-\widetilde{V}\widetilde{\Lambda})$
- Jellium：$v_{\text{eff},\sigma\sigma'} = v(\mathbf{q})/\varepsilon(\mathbf{q},E)$——动态屏蔽（习题 6.3.3）；与 $T=0$ 的 (5.189) 对应
- Hubbard：$v_{\text{eff},\uparrow\uparrow} = \frac{U^2}{N^2}\Lambda_{\downarrow\downarrow}/\left((1-\frac{U}{N}\Lambda_{\downarrow\uparrow})(1-\frac{U}{N}\Lambda_{\uparrow\downarrow}) - \frac{U^2}{N^2}\Lambda_{\downarrow\downarrow}\Lambda_{\uparrow\uparrow}\right)$；$v_{\text{eff},\downarrow\uparrow} = \frac{U}{N}(1-\frac{U}{N}\Lambda_{\downarrow\uparrow})/\text{同分母}$；一般 $v_{\text{eff},\sigma\sigma}\neq 0$（习题 6.3.4）；参数情形与 RPA 见正文
- 骨架图方案：抑制含极化贡献的自能图，骨架图中裸相互作用 → 有效相互作用；再以完整传播子自洽敷饰（Dyson 方程）
- 双重计数危险：图 6.42 第一图相互作用线不得敷饰（否则图 6.43 出现两次）；直接项（Part II 图 6.19）不重复出现（已含于第二求和项）

## 顶角函数与梯近似

- 顶角贡献 = 极化贡献中带两条粒子线 + 一条相互作用线的图部分（零阶 = 单个顶点）；不可约 = 不能分出完整自能图或完整极化图
- 顶角函数 $\Gamma_{\sigma\sigma'}(\mathbf{q}E_0;\mathbf{k}E)$ = 所有不可约顶角贡献之和；零阶贡献 $\delta_{\sigma\sigma'}$；记号：右上接相互作用线，$\mathbf{q},E_0$ 为转移量，左上外传播子 $(\mathbf{k}E\sigma')$，左下 $(\mathbf{k}+\mathbf{q},E+E_0,\sigma')$
- 极化传播子的顶角表示：$-\hbar\Lambda_{\mathbf{q}\sigma\sigma'}(E_0) = \frac{\varepsilon}{\hbar\beta}\sum_{\mathbf{k}E}(-G^M)(-G^M)\Gamma_{\sigma'\sigma}$（用完整传播子）
- 零级近似（$\Gamma\to\delta$，$G^M\to G^{0,M}$）⟹ 回到 $\Lambda_{\mathbf{q}}^{(0)}$
- 梯近似（Jellium 型相互作用；Hubbard 全零）：$\Gamma^{LA} = \delta + \delta(-\frac{1}{\hbar})\frac{1}{\hbar\beta}\sum v(-G^M)(-G^M)\Gamma^{LA}$（图 6.47）
- 常数相互作用 $v_0 = \frac{1}{N}\sum_{\mathbf{k}}v(\mathbf{k})$：$\Gamma^{LA}_\sigma = 1/(1-\varepsilon v_0\widehat{\Lambda}_{\mathbf{q}\sigma})$；$\Lambda^{LA}_{\mathbf{q}\sigma\sigma'} = \delta_{\sigma\sigma'}\widehat{\Lambda}_{\mathbf{q}\sigma}/(1-\varepsilon v_0\widehat{\Lambda}_{\mathbf{q}\sigma})$；$\Lambda^{LA}_{\mathbf{q}} = \sum_\sigma\frac{\widehat{\Lambda}_{\mathbf{q}\sigma}}{1-\varepsilon v_0\widehat{\Lambda}_{\mathbf{q}\sigma}}$
- $\widehat{\Lambda}_{\mathbf{q}\sigma}(E_0) = \frac{1}{\hbar^2}\sum_{\mathbf{p}}\int\!\!\int dxdy\,\frac{S_{\mathbf{p}\sigma}(x)S_{\mathbf{p}+\mathbf{q}\sigma}(y)}{\mathrm{i}E_0+x-y}(f_\varepsilon(x)-f_\varepsilon(y))$（习题 6.3.5；谱密度由单粒子 Dyson 方程自洽确定；$S\to S^{(0)}$ 回到 $\Lambda^{(0)}$）

## 横向自旋磁化率

- 定义：$\chi^{+-}_{ij}(E_0) = -\frac{\gamma}{\hbar^2}\langle\!\langle\sigma_i^+;\sigma_j^-\rangle\!\rangle^M_{E_0}$，$\gamma = \mu_0 g^2\mu_B^2/(V\hbar)$；$\sigma_i^+ = \hbar a_{i\uparrow}^+a_{i\downarrow}$ 等（Wannier 表象）
- $\chi^{+-}_{\mathbf{q}} = \frac{1}{N}\sum_{ij}\chi^{+-}_{ij}e^{\mathrm{i}\mathbf{q}\cdot(\mathbf{R}_i-\mathbf{R}_j)} = -\frac{\gamma}{N}\widehat{\chi}_{\mathbf{q}}$；$\widehat{\chi}_{\mathbf{q}} = \sum_{\mathbf{p}\mathbf{k}}\langle\!\langle a_{\mathbf{k}\uparrow}^+a_{\mathbf{k}+\mathbf{q}\downarrow};a_{\mathbf{p}\downarrow}^+a_{\mathbf{p}-\mathbf{q}\uparrow}\rangle\!\rangle^M_{E_0}$
- 外接传播子自旋成对不同 ⟹ Dyson 第二项为零：$-\widehat{\chi}_{\mathbf{q}} \equiv -\hbar\Lambda_{\mathbf{q}}^{\uparrow\downarrow}$（自旋磁化率 = 极化传播子）
- 精确表示：$-\widehat{\chi}_{\mathbf{q}}(E_0) = \frac{\varepsilon}{\hbar\beta}\sum_{\mathbf{k}E}G_{\mathbf{k}\uparrow}^MG_{\mathbf{k}+\mathbf{q}\downarrow}^M\Gamma_{\uparrow\downarrow}$
- Hubbard 梯近似可精确解：$\widehat{\chi}^{LA}_{\mathbf{q}}(E_0) = \frac{\hbar\widehat{\Lambda}_{\mathbf{q}\uparrow\downarrow}(E_0)}{1 + \frac{U}{N}\widehat{\Lambda}_{\mathbf{q}\uparrow\downarrow}(E_0)}$，其中 $\widehat{\Lambda}_{\mathbf{q}\uparrow\downarrow}$ 由谱表示式（自旋指标调整）给出
- 与 EoM 法结果 (4.183)（$\chi^{+-} = (\chi^{+-})^{(S)}/(1-\gamma^{-1}U(\chi^{+-})^{(S)})$）结构一致——RPA 型无穷粒子-空穴链求和
