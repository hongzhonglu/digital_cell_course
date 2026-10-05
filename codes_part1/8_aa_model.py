import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ----------------------------
# 1. Kinetic parameters (with PPP added)
# ----------------------------
params = {
    # Glycolysis
    'Vmax_HK':   10.0,  'Km_HK_Glc': 1.0,  'Km_HK_ATP': 0.8,
    'Vmax_PFK':  8.0,   'Km_PFK_G6P': 0.5, 'Km_PFK_ATP': 0.6,
    'Vmax_GAPD': 12.0,  'Km_GAPD_F16BP': 0.8, 'Km_GAPD_NAD': 0.3,
    'Vmax_PGK':  12.0,  'Km_PGK_3PG': 0.7, 'Km_PGK_ADP': 0.9,
    'Vmax_PYK':  10.0,  'Km_PYK_PEP': 1.2, 'Km_PYK_ADP': 1.0,

    # Ethanol fermentation
    'Vmax_PDC':  4.0,   'Km_PDC_Pyr': 0.6,
    'Vmax_ADH':  5.0,   'Km_ADH_AcAld': 0.1, 'Km_ADH_NADH': 0.2,

    # Cytosolic AcCoA for AA (minor in yeast, but engineered)
    'Vmax_AcCoA_cyt': 0.5, 'Km_AcCoA_cyt_Pyr': 0.4, 'Km_AcCoA_cyt_NAD': 0.5,

    # Mitochondrial PDH (for TCA)
    'Vmax_PDH_mito': 3.0, 'Km_PDH_mito_Pyr': 0.3, 'Km_PDH_mito_NAD': 0.4,

    # TCA cycle (lumped: AcCoA_mito → 3NADH + FADH2 + GTP)
    'Vmax_TCA': 2.5, 'Km_TCA_AcCoA': 0.05,

    # OXPHOS: driven by NADH and FADH2
    'k_OXPHOS_NADH': 1.0,  # rate constant for NADH oxidation
    'k_OXPHOS_FADH2': 0.6, # rate constant for FADH2 oxidation

    # Artemisinic acid pathway
    'Vmax_MVA':  0.8,   'Km_MVA_AcCoA': 0.1,
    'Vmax_ADS':  0.6,   'Km_ADS_FPP': 0.05,
    'Vmax_CYP':  0.5,   'Km_CYP_Amor': 0.02,

    # Pentose Phosphate Pathway (oxidative branch, lumped)
    'Vmax_G6PDH': 2.0,      # moderate flux
    'Km_G6PDH_G6P': 0.2,
    'Km_G6PDH_NADP': 0.05,  # high affinity for NADP+
}

# Conserved pools (mM)
ATP_total = 5.0      # increased due to respiration
NAD_total = 1.5      # higher total redox capacity
NADP_total = 0.6
FAD_total = 0.2      # FAD + FADH2

# Initial concentrations (mM)
y0 = [
    5.0,    # Glc
    0.1,    # G6P
    0.05,   # F16BP
    0.1,    # 3PG
    0.05,   # PEP
    0.01,   # Pyr
    0.0,    # Acetaldehyde
    0.0,    # Ethanol
    0.001,  # Cytosolic AcCoA (for AA)
    0.0,    # FPP
    0.0,    # Amoradiene
    0.0,    # Artemisinic acid (AA)
    0.001,  # Mitochondrial AcCoA (AcCoA_mito)
    3.5,    # ATP
    1.0,    # NAD+
    0.4,    # NADPH
    0.15,   # FAD (oxidized form)
]

names = ['Glc','G6P','F16BP','3PG','PEP','Pyr','AcAld','EtOH',
         'AcCoA_cyt','FPP','Amor','AA','AcCoA_mito',
         'ATP','NAD+','NADPH','FAD']

# ----------------------------
# 2. Reaction rates (with PPP)
# ----------------------------
def reaction_rates(y, params):
    (Glc, G6P, F16BP, _3PG, PEP, Pyr, AcAld, EtOH,
     AcCoA_cyt, FPP, Amor, AA, AcCoA_mito,
     ATP, NAD, NADPH, FAD) = y

    ADP = ATP_total - ATP
    NADH = NAD_total - NAD
    NADP = NADP_total - NADPH
    FADH2 = FAD_total - FAD

    # --- Glycolysis ---
    v_HK = (params['Vmax_HK'] * Glc * ATP) / ((params['Km_HK_Glc'] + Glc) * (params['Km_HK_ATP'] + ATP))
    v_PFK = (params['Vmax_PFK'] * G6P * ATP) / ((params['Km_PFK_G6P'] + G6P) * (params['Km_PFK_ATP'] + ATP))
    v_GAPD = (params['Vmax_GAPD'] * F16BP * NAD) / ((params['Km_GAPD_F16BP'] + F16BP) * (params['Km_GAPD_NAD'] + NAD))
    v_PGK = (params['Vmax_PGK'] * _3PG * ADP) / ((params['Km_PGK_3PG'] + _3PG) * (params['Km_PGK_ADP'] + ADP))
    v_PYK = (params['Vmax_PYK'] * PEP * ADP) / ((params['Km_PYK_PEP'] + PEP) * (params['Km_PYK_ADP'] + ADP))

    # --- Ethanol fermentation ---
    v_PDC = (params['Vmax_PDC'] * Pyr) / (params['Km_PDC_Pyr'] + Pyr)
    v_ADH = (params['Vmax_ADH'] * AcAld * NADH) / (
        (params['Km_ADH_AcAld'] + AcAld) * (params['Km_ADH_NADH'] + NADH)
    )

    # --- Cytosolic AcCoA (for AA pathway) ---
    v_AcCoA_cyt = (params['Vmax_AcCoA_cyt'] * Pyr * NAD) / (
        (params['Km_AcCoA_cyt_Pyr'] + Pyr) * (params['Km_AcCoA_cyt_NAD'] + NAD)
    )

    # --- Mitochondrial PDH ---
    v_PDH_mito = (params['Vmax_PDH_mito'] * Pyr * NAD) / (
        (params['Km_PDH_mito_Pyr'] + Pyr) * (params['Km_PDH_mito_NAD'] + NAD)
    )

    # --- TCA cycle (lumped) ---
    v_TCA = (params['Vmax_TCA'] * AcCoA_mito) / (params['Km_TCA_AcCoA'] + AcCoA_mito)

    # --- OXPHOS ---
    v_OXPHOS_NADH = params['k_OXPHOS_NADH'] * NADH
    v_OXPHOS_FADH2 = params['k_OXPHOS_FADH2'] * FADH2

    # --- AA pathway ---
    v_MVA = (params['Vmax_MVA'] * AcCoA_cyt) / (params['Km_MVA_AcCoA'] + AcCoA_cyt)
    v_ADS = (params['Vmax_ADS'] * FPP) / (params['Km_ADS_FPP'] + FPP)
    v_CYP = (params['Vmax_CYP'] * Amor) / (params['Km_CYP_Amor'] + Amor)

    # --- PPP (oxidative branch) ---
    v_PPP = (params['Vmax_G6PDH'] * G6P * NADP) / (
        (params['Km_G6PDH_G6P'] + G6P) * (params['Km_G6PDH_NADP'] + NADP)
    )

    return (v_HK, v_PFK, v_GAPD, v_PGK, v_PYK, v_PDC, v_ADH,
            v_AcCoA_cyt, v_PDH_mito, v_TCA,
            v_OXPHOS_NADH, v_OXPHOS_FADH2,
            v_MVA, v_ADS, v_CYP, v_PPP,
            ADP, NADH, NADP, FADH2)

# ----------------------------
# 3. ODE system with TCA + OXPHOS + PPP
# ----------------------------
def full_respiratory_model_ode(t, y, params):
    (Glc, G6P, F16BP, _3PG, PEP, Pyr, AcAld, EtOH,
     AcCoA_cyt, FPP, Amor, AA, AcCoA_mito,
     ATP, NAD, NADPH, FAD) = y

    (v_HK, v_PFK, v_GAPD, v_PGK, v_PYK, v_PDC, v_ADH,
     v_AcCoA_cyt, v_PDH_mito, v_TCA,
     v_OXPHOS_NADH, v_OXPHOS_FADH2,
     v_MVA, v_ADS, v_CYP, v_PPP,
     ADP, NADH, NADP, FADH2) = reaction_rates(y, params)

    dGlc_dt       = -v_HK
    dG6P_dt       = v_HK - v_PFK - v_PPP          # <-- PPP分流
    dF16BP_dt     = v_PFK - v_GAPD
    d3PG_dt       = v_GAPD - v_PGK
    dPEP_dt       = v_PGK - v_PYK
    dPyr_dt       = v_PYK - v_PDC - v_AcCoA_cyt - v_PDH_mito
    dAcAld_dt     = v_PDC - v_ADH
    dEtOH_dt      = v_ADH
    dAcCoA_cyt_dt = v_AcCoA_cyt - v_MVA
    dFPP_dt       = v_MVA - v_ADS
    dAmor_dt      = v_ADS - v_CYP
    dAA_dt        = v_CYP
    dAcCoA_mito_dt= v_PDH_mito - v_TCA

    # ATP balance:
    dATP_dt = (-v_HK - v_PFK + v_PGK + v_PYK) \
              + v_TCA \
              + 2.5 * v_OXPHOS_NADH + 1.5 * v_OXPHOS_FADH2 \
              - 3 * v_MVA

    # NAD+ balance:
    dNAD_dt = -v_GAPD + v_ADH - v_AcCoA_cyt - v_PDH_mito + v_OXPHOS_NADH

    # NADPH balance: PPP produces 2 NADPH per G6P
    dNADPH_dt = 2 * v_PPP - 2 * v_MVA - 3 * v_CYP  # <-- key update

    # FAD balance
    dFAD_dt = -v_TCA + v_OXPHOS_FADH2

    return [dGlc_dt, dG6P_dt, dF16BP_dt, d3PG_dt, dPEP_dt, dPyr_dt,
            dAcAld_dt, dEtOH_dt,
            dAcCoA_cyt_dt, dFPP_dt, dAmor_dt, dAA_dt, dAcCoA_mito_dt,
            dATP_dt, dNAD_dt, dNADPH_dt, dFAD_dt]

# ----------------------------
# 4. Simulation
# ----------------------------
t_span = (0, 30)
t_eval = np.linspace(0, 30, 400)

sol = solve_ivp(
    full_respiratory_model_ode,
    t_span,
    y0,
    args=(params,),
    t_eval=t_eval,
    method='LSODA'
)

# Derived variables
ATP_vals = sol.y[13]
ADP_vals = ATP_total - ATP_vals
NAD_vals = sol.y[14]
NADH_vals = NAD_total - NAD_vals
NADPH_vals = sol.y[15]
AA_vals = sol.y[11]
EtOH_vals = sol.y[7]
AcCoA_mito_vals = sol.y[12]
G6P_vals = sol.y[1]
PPP_flux = (params['Vmax_G6PDH'] * G6P_vals * (NADP_total - NADPH_vals)) / (
    (params['Km_G6PDH_G6P'] + G6P_vals) * (params['Km_G6PDH_NADP'] + (NADP_total - NADPH_vals))
)

# ----------------------------
# 5. Plotting
# ----------------------------
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Central carbon
axes[0, 0].plot(sol.t, sol.y[:6].T)
axes[0, 0].set_title('Central Metabolites')
axes[0, 0].legend(names[:6])
axes[0, 0].grid(True, alpha=0.4)

# Fermentation vs respiration
axes[0, 1].plot(sol.t, EtOH_vals, 'brown', label='Ethanol')
axes[0, 1].plot(sol.t, AcCoA_mito_vals, 'g--', label='mito-AcCoA')
axes[0, 1].set_title('Carbon Branching: Fermentation vs Respiration')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.4)

# AA pathway
axes[0, 2].plot(sol.t, sol.y[8], 'c--', label='AcCoA_cyt')
axes[0, 2].plot(sol.t, sol.y[9], 'b-.', label='FPP')
axes[0, 2].plot(sol.t, sol.y[10], 'm:', label='Amor')
axes[0, 2].plot(sol.t, AA_vals, 'r-', linewidth=2, label='Artemisinic Acid')
axes[0, 2].set_title('Artemisinic Acid Synthesis (with PPP)')
axes[0, 2].legend()
axes[0, 2].grid(True, alpha=0.4)

# Energy charge
axes[1, 0].plot(sol.t, ATP_vals, label='ATP')
axes[1, 0].plot(sol.t, ADP_vals, label='ADP')
axes[1, 0].set_title('ATP/ADP (Higher due to Respiration)')
axes[1, 0].set_xlabel('Time (h)')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.4)

# Redox cofactors
axes[1, 1].plot(sol.t, NAD_vals, label='NAD+')
axes[1, 1].plot(sol.t, NADH_vals, label='NADH')
axes[1, 1].plot(sol.t, NADPH_vals, label='NADPH', linewidth=2)
axes[1, 1].set_title('Redox Balance (PPP boosts NADPH)')
axes[1, 1].set_xlabel('Time (h)')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.4)

# Final products
axes[1, 2].plot(sol.t, EtOH_vals, 'brown', label='Ethanol')
axes[1, 2].plot(sol.t, AA_vals, 'r--', label='Artemisinic Acid')
axes[1, 2].set_title('Final Products')
axes[1, 2].set_xlabel('Time (h)')
axes[1, 2].set_ylabel('Concentration (mM)')
axes[1, 2].legend()
axes[1, 2].grid(True, alpha=0.4)

plt.tight_layout()
plt.show()

# ----------------------------
# 6. Summary
# ----------------------------
print(f"Final artemisinic acid: {AA_vals[-1]:.5f} mM")
print(f"Final ethanol:          {EtOH_vals[-1]:.3f} mM")
print(f"Final ATP:              {ATP_vals[-1]:.3f} mM")
print(f"Final NADPH:            {NADPH_vals[-1]:.4f} mM")
print(f"Final NAD+:             {NAD_vals[-1]:.3f} mM")
print(f"mito-AcCoA peak:        {np.max(AcCoA_mito_vals):.4f} mM")
print(f"G6P diverted to PPP (final flux): {PPP_flux[-1]:.4f} mM/h")










### gemini 3 ###
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# --- 1. 定义代谢网络动力学方程 ---
def metabolic_model(y, t, params):
    """
    整合 EMP, PPP, Fermentation, TCA, Artemisinic Acid (AA) 的动力学模型
    """
    # 解包状态变量
    Glc, G6P, Pyr, AcCoA, AA, ATP, NADH, NADPH = y

    # 辅因子守恒计算 (Conservation of Mass)
    # 假设细胞内总库是恒定的
    ADP = params['Atot'] - ATP
    NAD = params['Ntot'] - NADH
    NADP = params['NPtot'] - NADPH

    # 防止数值误差导致的负浓度
    ADP = max(ADP, 1e-4);
    NAD = max(NAD, 1e-4);
    NADP = max(NADP, 1e-4)

    # --- 反应速率 (Flux Rates) ---

    # 1. 葡萄糖摄取 (Uptake)
    v_uptake = params['V_in'] * Glc / (params['Km_glc'] + Glc) * (ADP / (0.5 + ADP))  # ATP低时促进摄取

    # 2. EMP 途径 (Glycolysis): G6P -> Pyr
    # 产物: 1 Pyr, 1 ATP, 1 NADH (集总系数)
    # 调节: 被 ATP 抑制 (巴斯德效应)
    v_emp = params['k_emp'] * G6P * ADP * NAD * (1 / (1 + (ATP / params['Ki_atp']) ** 2))

    # 3. PPP 途径 (Pentose Phosphate Pathway): G6P -> ... -> GAP/F6P
    # 关键功能: 生成 NADPH
    # 简化: G6P -> Pyr (碳流最终回到下划线) + 2 NADPH
    # 调节: 被 NADPH 强烈抑制
    v_ppp = params['k_ppp'] * G6P * NADP * (1 / (1 + (NADPH / params['Ki_nadph'])))

    # 4. 乙醇发酵 (Fermentation): Pyr -> Ethanol
    # 消耗 NADH，再生 NAD+
    v_ferm = params['k_ferm'] * Pyr * NADH

    # 5. PDH (丙酮酸脱氢酶): Pyr -> AcCoA (进入线粒体)
    # 生成 NADH
    v_pdh = params['k_pdh'] * Pyr * NAD

    # 6. TCA + OXPHOS (三羧酸循环 + 氧化磷酸化): AcCoA -> CO2
    # 这是一个巨大的能量工厂: 消耗 NADH, 产生大量 ATP
    # 简化: AcCoA -> 消耗 NADH -> 生成 ATP
    v_resp = params['k_resp'] * AcCoA * ADP * (NADH / (1 + NADH))

    # 7. 青蒿酸合成 (AA Synthesis via MVA pathway): AcCoA -> AA
    # 这是一个极其消耗 NADPH 和 ATP 的过程
    # 假设: 3 AcCoA + 3 ATP + 2 NADPH -> 1 AA (简化化学计量数)
    v_aa = params['k_aa'] * AcCoA * (ATP / (1 + ATP)) * (NADPH / (params['Km_nadph_aa'] + NADPH))

    # 8. 基础代谢消耗 (ATP Maintenance & NADPH Oxidation)
    v_atp_use = params['k_atp_use'] * ATP
    v_nadph_use = params['k_nadph_leak'] * NADPH  # 比如抗氧化压力消耗

    # --- 微分方程 (Mass Balance) ---

    dGlc = -v_uptake

    # G6P: 来源 Uptake, 去路 EMP + PPP
    dG6P = v_uptake - v_emp - v_ppp

    # Pyr: 来源 EMP + PPP(简化归流), 去路 Ferm + PDH
    dPyr = v_emp + v_ppp - v_ferm - v_pdh

    # AcCoA: 来源 PDH, 去路 TCA + AA合成
    dAcCoA = v_pdh - v_resp - 3 * v_aa  # 假设合成需3分子AcCoA

    # AA (青蒿酸): 产物积累
    dAA = v_aa

    # 辅因子动力学
    # dATP: EMP(+), PPP(0), Ferm(0), Resp(+++), AA(-), Use(-)
    dATP = 1 * v_emp + 10 * v_resp - 3 * v_aa - v_atp_use

    # dNADH: EMP(+), PDH(+), Ferm(-), Resp(-)
    dNADH = 1 * v_emp + 1 * v_pdh - 1 * v_ferm - 2 * v_resp  # Resp消耗NADH用于产ATP

    # dNADPH: PPP(++), AA(-), Leak(-)
    dNADPH = 2 * v_ppp - 2 * v_aa - v_nadph_use

    return [dGlc, dG6P, dPyr, dAcCoA, dAA, dATP, dNADH, dNADPH]


# --- 2. 模拟设置与运行 ---

# 基础参数集
base_params = {
    'V_in': 2.0, 'Km_glc': 0.5,
    'k_emp': 0.8, 'Ki_atp': 2.0,  # EMP 强力，受 ATP 抑制
    'k_ppp': 0.1, 'Ki_nadph': 3.0,  # 默认 PPP 很弱 (Wait to change this!)
    'k_ferm': 0.5,  # 发酵能力
    'k_pdh': 0.4,  # 进入线粒体能力
    'k_resp': 0.3,  # 呼吸链能力
    'k_aa': 0.3, 'Km_nadph_aa': 1.0,  # AA 合成酶 (P450 等)
    'k_atp_use': 1.5,
    'k_nadph_leak': 0.05,

    # 总库容量
    'Atot': 10.0,
    'Ntot': 2.0,
    'NPtot': 2.0
}

# 初始条件
# [Glc, G6P, Pyr, AcCoA, AA, ATP, NADH, NADPH]
y0 = [20.0, 0.5, 0.5, 0.1, 0.0, 8.0, 0.2, 0.5]
t = np.linspace(0, 60, 1000)

# --- 场景对比测试 ---

# 场景 A: 野生型 (低 PPP 活性)
params_low_ppp = base_params.copy()
params_low_ppp['k_ppp'] = 0.05  # PPP 弱，NADPH 供应不足
sol_low = odeint(metabolic_model, y0, t, args=(params_low_ppp,))

# 场景 B: 代谢工程改造株 (高 PPP 活性 - 强化 NADPH 供应)
params_high_ppp = base_params.copy()
params_high_ppp['k_ppp'] = 0.8  # 提高 G6P 脱氢酶活性，导流向 PPP
sol_high = odeint(metabolic_model, y0, t, args=(params_high_ppp,))

# --- 3. 结果可视化 ---
fig, axs = plt.subplots(2, 2, figsize=(14, 10))

# 图 1: 青蒿酸产量对比 (核心目标)
axs[0, 0].plot(t, sol_low[:, 4], 'k--', label='Low PPP Flux (Wild Type)', linewidth=2)
axs[0, 0].plot(t, sol_high[:, 4], 'g-', label='High PPP Flux (Optimized)', linewidth=2)
axs[0, 0].set_title('Target Yield: Artemisinic Acid (AA)', fontsize=14, fontweight='bold')
axs[0, 0].set_ylabel('Concentration (mM)')
axs[0, 0].legend()
axs[0, 0].grid(True, alpha=0.3)
axs[0, 0].fill_between(t, sol_low[:, 4], sol_high[:, 4], color='green', alpha=0.1)  # 强调增量

# 图 2: NADPH 供应水平 (机制解释)
axs[0, 1].plot(t, sol_low[:, 7], 'k--', label='NADPH (Low PPP)')
axs[0, 1].plot(t, sol_high[:, 7], 'g-', label='NADPH (High PPP)')
axs[0, 1].axhline(y=base_params['Km_nadph_aa'], color='r', linestyle=':', label='AA Synthase Km_NADPH')
axs[0, 1].set_title('Cofactor Availability: NADPH', fontsize=14)
axs[0, 1].set_ylabel('Concentration (mM)')
axs[0, 1].set_ylim(0, 2.0)
axs[0, 1].legend()
axs[0, 1].grid(True, alpha=0.3)
axs[0, 1].text(5, 0.2, "NADPH limiting\nAA synthesis here", color='black', fontsize=9)

# 图 3: 碳流节点分流 (G6P) - 观察 EMP 和 PPP 的竞争
# 需要反算 Flux
v_emp_low = [base_params['k_emp'] * y[1] * (base_params['Atot'] - y[5]) * (base_params['Ntot'] - y[6]) for y in sol_low]
v_ppp_high = [params_high_ppp['k_ppp'] * y[1] * (base_params['NPtot'] - y[7]) for y in sol_high]
axs[1, 0].plot(t, v_emp_low, 'b-', label='Glycolysis Flux (EMP)', alpha=0.6)
axs[1, 0].plot(t, v_ppp_high, 'orange', label='Pentose Phosphate Flux (PPP)', alpha=0.6)
axs[1, 0].set_title('Carbon Flux Distribution (in Optimized Strain)', fontsize=14)
axs[1, 0].set_ylabel('Flux Rate')
axs[1, 0].legend()
axs[1, 0].grid(True, alpha=0.3)

# 图 4: 能量状态 ATP
axs[1, 1].plot(t, sol_low[:, 5], 'k--', label='ATP (Low PPP)')
axs[1, 1].plot(t, sol_high[:, 5], 'g-', label='ATP (High PPP)')
axs[1, 1].set_title('Energy Status (ATP)', fontsize=14)
axs[1, 1].set_ylabel('Concentration (mM)')
axs[1, 1].set_xlabel('Time (h)')
axs[1, 1].legend()
axs[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()





import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


class BioParameters:
    def __init__(self, ppp_strength=1.0):
        # --- 初始底物与辅因子池 (mM) ---
        self.NADP_pool = 3.0  # NADPH + NADP+
        self.AdN_pool = 6.0  # ATP + ADP + AMP
        self.NAD_pool = 12.0  # NADH + NAD+ (假设胞质与线粒体交换较快，简化为一个池)

        # --- 1. 糖酵解 (Glycolysis) ---
        # HK: Glucose + ATP -> G6P
        self.vmax_HK = 25.0
        self.km_HK_Glc = 0.1
        self.km_HK_ATP = 0.4
        self.ki_HK_G6P = 5.0  # 产物抑制

        # PGI: G6P <-> F6P (可逆，快速)
        self.vmax_PGI = 100.0
        self.Keq_PGI = 0.3  # [F6P]/[G6P] 平衡常数

        # PFK: F6P + ATP -> FBP (关键限速，变构调节)
        self.vmax_PFK = 20.0
        self.km_PFK_F6P = 0.5
        self.ki_PFK_ATP = 1.5  # ATP 变构抑制 (巴斯德效应)

        # ALDO/TPI/GAPDH (Lower Glycolysis Lumped):
        # FBP -> ... -> 2 GAP -> 2 PEP + 2 NADH + 2 ATP
        self.vmax_LowGly = 40.0

        # PYK: PEP + ADP -> Pyruvate + ATP
        self.vmax_PYK = 35.0

        # --- 2. 磷酸戊糖途径 (PPP) ---
        # 氧化段 (G6PDH + 6PGDH): G6P + 2 NADP+ -> Ru5P + 2 NADPH + CO2
        # *ppp_strength 用于模拟工程菌株过表达 G6PDH*
        self.vmax_PPP_Ox = 4.0 * ppp_strength
        self.km_PPP_G6P = 0.8
        self.km_PPP_NADP = 0.1
        self.ki_PPP_NADPH = 0.2  # NADPH 强反馈抑制

        # 非氧化段 (Transketolase/Transaldolase):
        # 3 Ru5P -> 2 F6P + 1 GAP (碳架回补)
        self.k_PPP_NonOx = 10.0  # 一级反应速率常数 (简化)

        # --- 3. 发酵与 PDH Bypass (细胞质乙酰辅酶A来源) ---
        # PDC: Pyruvate -> Acetaldehyde + CO2
        self.vmax_PDC = 20.0

        # ADH: Acetaldehyde + NADH <-> Ethanol + NAD+
        self.vmax_ADH = 25.0
        self.Keq_ADH = 1.0e4  # 倾向于生成乙醇

        # ALD (Aldehyde Dehydrogenase) + ACS (Acetyl-CoA Synthetase):
        # Acetaldehyde -> Acetate -> Acetyl-CoA(cyt) (消耗 ATP, 生成 NADH)
        # 这是青蒿酸合成的前体来源途径 (PDH Bypass)
        self.vmax_ACS = 5.0
        self.km_ACS_Ald = 0.2

        # --- 4. 线粒体呼吸 ---
        # PDH (Mito): Pyruvate -> Acetyl-CoA(mito) + NADH
        self.vmax_PDH = 10.0

        # TCA: Acetyl-CoA(mito) -> ... -> 3 NADH + 1 FADH2 + 1 GTP(ATP)
        self.vmax_TCA = 8.0

        # OXPHOS: NADH + 0.5 O2 -> 2.5 ATP (P/O ratio ~ 2.5)
        # 假设氧气饱和，简化为 NADH 依赖
        self.vmax_OXPHOS = 30.0

        # --- 5. 青蒿酸合成 (MVA Pathway + P450) ---
        # MVA 上游: 3 Acetyl-CoA(cyt) + 3 ATP + 2 NADPH -> FPP
        self.vmax_MVA = 2.5
        self.km_MVA_AcCoA = 0.1
        self.km_MVA_NADPH = 0.05

        # 下游合成: FPP -> Amorphadiene -> Artemisinic Acid
        # 消耗: 3 NADPH, 3 O2
        self.vmax_ADS = 2.0  # FPP -> Amorphadiene
        self.vmax_CYP = 1.5  # Amorphadiene -> AA (限速步骤)
        self.km_CYP_NADPH = 0.05

        # --- 6. 细胞生长 (Biomass) ---
        # 简化方程: 消耗 G6P, GAP, Pyruvate, AcCoA, ATP, NADPH
        self.mu_max = 0.3
        self.Kd = 0.005  # 死亡率
        self.k_maint = 1.0  # 维持能


def kinetic_model(y, t, p):
    # 状态变量解包
    (Glc, G6P, F6P, GAP, PEP, Pyr, Acal, Eth,
     AcCoA_cyt, AcCoA_mit, FPP, Amor, AA,
     Biomass, ATP, NADPH, NADH) = y

    # --- 辅助变量 (守恒约束) ---
    ADP = max(0.01, p.AdN_pool - ATP)
    NADP = max(0.01, p.NADP_pool - NADPH)
    NAD = max(0.01, p.NAD_pool - NADH)

    # --- 1. 反应速率计算 (Fluxes) ---

    # HK: 消耗 ATP
    v_HK = p.vmax_HK * (Glc / (p.km_HK_Glc + Glc)) * (ATP / (p.km_HK_ATP + ATP)) * (
                1 / (1 + G6P / p.ki_HK_G6P)) * Biomass

    # PPP 氧化段: 生成 2 NADPH, 1 CO2 per Flux, 生成 Ru5P (简化直接进入非氧化)
    # 抑制项: (1 + NADPH/Ki)
    v_PPP_Ox = p.vmax_PPP_Ox * (G6P / (p.km_PPP_G6P + G6P)) * (NADP / (p.km_PPP_NADP + NADP)) / (
                1 + NADPH / p.ki_PPP_NADPH) * Biomass

    # PPP 非氧化段 (简化): 假设 Ru5P 瞬时转化为 F6P 和 GAP
    # 化学计量: 3 Ru5P -> 2 F6P + 1 GAP
    # 这里我们不显式模拟 Ru5P，直接将 v_PPP_Ox 的产物分配给 F6P/GAP
    # Flux_recycling = v_PPP_Ox (以 G6P 为基准)

    # PGI: G6P <-> F6P
    v_PGI = (p.vmax_PGI / p.km_HK_Glc) * (G6P - F6P / p.Keq_PGI) * Biomass
    # 注意：上面用了简化的质量作用定律近似

    # PFK: F6P + ATP -> FBP + ADP (消耗 ATP)
    # 强 ATP 抑制 (1 + (ATP/Ki)^n)
    allostery = 1 / (1 + (ATP / p.ki_PFK_ATP) ** 4)
    v_PFK = p.vmax_PFK * (F6P / (p.km_PFK_F6P + F6P)) * (ATP / (0.1 + ATP)) * allostery * Biomass

    # Lower Glycolysis: FBP -> 2 GAP -> 2 PEP
    # 产出: 2 ATP (PGK), 2 NADH (GAPDH) per FBP (即 per PFK flux)
    # 为了动态稳定性，我们将 FBP 视为瞬间稳态，直接写 GAP 动力学
    # 假设 PFK 产物直接变为 GAP: 1 FBP -> 2 GAP
    v_GAP_to_PEP = p.vmax_LowGly * (GAP / (0.1 + GAP)) * (NAD / (0.1 + NAD)) * (ADP / (0.1 + ADP)) * Biomass

    # PYK: PEP -> Pyruvate (产 1 ATP)
    v_PYK = p.vmax_PYK * (PEP / (0.1 + PEP)) * (ADP / (0.1 + ADP)) * Biomass

    # PDC: Pyruvate -> Acetaldehyde + CO2
    v_PDC = p.vmax_PDC * (Pyr / (1.0 + Pyr)) * Biomass

    # ADH: Acal + NADH <-> Ethanol + NAD
    v_ADH = p.vmax_ADH * (Acal * NADH - Eth * NAD / p.Keq_ADH) / (1 + Acal + Eth) * Biomass

    # PDH Bypass (ACS): Acal -> -> Acetyl-CoA(cyt)
    # 简化步骤：Acal + NAD + ATP -> AcCoA_cyt + NADH + AMP + PPi
    # 这是一个耗能步骤！
    v_ACS = p.vmax_ACS * (Acal / (p.km_ACS_Ald + Acal)) * (ATP / (0.2 + ATP)) * (NAD / (0.2 + NAD)) * Biomass

    # PDH (Mito): Pyr -> AcCoA_mit (进入线粒体)
    v_PDH_mit = p.vmax_PDH * (Pyr / (0.5 + Pyr)) * (NAD / (0.2 + NAD)) * Biomass

    # TCA: AcCoA_mit -> Energy
    # 1 AcCoA -> 3 NADH + 1 FADH2(折算为NADH) + 1 ATP
    v_TCA = p.vmax_TCA * (AcCoA_mit / (0.1 + AcCoA_mit)) * (NAD / (0.2 + NAD)) * Biomass

    # OXPHOS: NADH -> ATP
    # 假设 P/O = 2.5
    v_OXPHOS = p.vmax_OXPHOS * (NADH / (0.1 + NADH)) * (ADP / (0.1 + ADP)) * Biomass

    # --- 青蒿酸途径 ---

    # MVA: 3 AcCoA_cyt + 3 ATP + 2 NADPH -> FPP
    v_MVA = p.vmax_MVA * (AcCoA_cyt / (p.km_MVA_AcCoA + AcCoA_cyt)) * (NADPH / (p.km_MVA_NADPH + NADPH)) ** 2 * (
                ATP / (0.2 + ATP)) * Biomass

    # ADS: FPP -> Amorphadiene
    v_ADS = p.vmax_ADS * (FPP / (0.05 + FPP)) * Biomass

    # CYP71AV1: Amor + 3 NADPH + 3 O2 -> AA
    v_CYP = p.vmax_CYP * (Amor / (0.05 + Amor)) * (NADPH / (p.km_CYP_NADPH + NADPH)) * Biomass

    # --- 细胞生长 ---
    # 假设生长需要: G6P, ATP, NADPH, Acetyl-CoA, 氮源等
    # 生长不仅消耗物质，还稀释胞内浓度
    limit_C = G6P / (0.5 + G6P)
    limit_E = ATP / (0.5 + ATP)
    limit_R = NADPH / (0.2 + NADPH)
    mu = p.mu_max * limit_C * limit_E * limit_R
    v_Growth = mu * Biomass

    # --- 2. 微分方程 (Mass Balance ODEs) ---

    # Glc
    dGlc = - v_HK

    # G6P: +HK, -PGI, -PPP_Ox, -Growth
    dG6P = v_HK - v_PGI - v_PPP_Ox - 0.5 * v_Growth

    # F6P: +PGI, -PFK, +PPP_Recycle (2/3 of PPP flux)
    # PPP 产生的 Ru5P 有 2/3 变成 F6P
    dF6P = v_PGI - v_PFK + (2.0 / 3.0) * v_PPP_Ox - 0.1 * v_Growth

    # GAP: +2*PFK(via FBP), -LowGly, +PPP_Recycle (1/3 of PPP flux)
    # PFK 产生 FBP，裂解为 2 GAP
    dGAP = (2.0 * v_PFK) - v_GAP_to_PEP + (1.0 / 3.0) * v_PPP_Ox

    # PEP: +LowGly, -PYK
    dPEP = v_GAP_to_PEP - v_PYK

    # Pyr: +PYK, -PDC, -PDH_mit, -Growth
    dPyr = v_PYK - v_PDC - v_PDH_mit - 0.5 * v_Growth

    # Acal: +PDC, -ADH, -ACS
    dAcal = v_PDC - v_ADH - v_ACS

    # Ethanol: +ADH
    dEth = v_ADH

    # AcCoA_cyt: +ACS, -MVA
    dAcCoA_cyt = v_ACS - (3.0 * v_MVA) - 0.2 * v_Growth

    # AcCoA_mit: +PDH_mit, -TCA
    dAcCoA_mit = v_PDH_mit - v_TCA

    # FPP: +MVA, -ADS
    dFPP = v_MVA - v_ADS

    # Amorphadiene: +ADS, -CYP
    dAmor = v_ADS - v_CYP

    # AA: +CYP
    dAA = v_CYP

    # Biomass
    dBiomass = v_Growth - p.Kd * Biomass

    # --- 3. 辅因子动态 (全局) ---

    # ATP Balance:
    # Prod: +2*v_LowGly(PGK), +1*v_PYK, +1*v_TCA, +2.5*v_OXPHOS
    # Cons: -1*HK, -1*PFK, -1*ACS, -3*MVA, -Growth, -Maint
    dATP = (2 * v_GAP_to_PEP + v_PYK + v_TCA + 2.5 * v_OXPHOS) - (
                v_HK + v_PFK + v_ACS + 3 * v_MVA + 10 * v_Growth + p.k_maint * Biomass)

    # NADPH Balance:
    # Prod: +2*v_PPP_Ox
    # Cons: -2*MVA, -3*CYP, -Growth
    dNADPH = (2.0 * v_PPP_Ox) - (2.0 * v_MVA) - (3.0 * v_CYP) - (2.0 * v_Growth)

    # NADH Balance:
    # Prod: +2*v_LowGly, +1*v_ACS, +1*v_PDH, +4*v_TCA (simplified TCA NADH yield)
    # Cons: -1*ADH (生成乙醇时), -v_OXPHOS
    # *注意 ADH 方向*: v_ADH 正向是消耗 NADH 生成乙醇
    dNADH = (2 * v_GAP_to_PEP + v_ACS + v_PDH_mit + 4 * v_TCA) - v_ADH - v_OXPHOS

    # 增加刚性修正防止浓度为负
    threshold = 1e-4
    if ATP < threshold and dATP < 0: dATP = 0
    if ATP > p.AdN_pool - threshold and dATP > 0: dATP = 0
    if NADPH < threshold and dNADPH < 0: dNADPH = 0
    if NADPH > p.NADP_pool - threshold and dNADPH > 0: dNADPH = 0
    if NADH < threshold and dNADH < 0: dNADH = 0
    if NADH > p.NAD_pool - threshold and dNADH > 0: dNADH = 0

    return [dGlc, dG6P, dF6P, dGAP, dPEP, dPyr, dAcal, dEth,
            dAcCoA_cyt, dAcCoA_mit, dFPP, dAmor, dAA,
            dBiomass, dATP, dNADPH, dNADH]


# ==========================================
# 执行模拟与绘图
# ==========================================

def run_simulation(ppp_factor, label_str):
    # 初始条件 (mM)
    y0 = [
        120.0,  # Glc
        0.5, 0.2, 0.1, 0.1, 0.5,  # G6P ... Pyr
        0.0, 0.0,  # Acal, Eth
        0.1, 0.5,  # AcCoA_cyt, AcCoA_mit
        0.0, 0.0, 0.0,  # FPP, Amor, AA
        0.05,  # Biomass
        5.0, 0.5, 1.0  # ATP, NADPH, NADH
    ]

    t = np.linspace(0, 50, 500)
    params = BioParameters(ppp_strength=ppp_factor)

    sol = odeint(kinetic_model, y0, t, args=(params,))
    return t, sol


# 运行两种场景
t, res_wt = run_simulation(1.0, "Wild Type")  # 正常 PPP
t, res_eng = run_simulation(6.0, "Engineered")  # 增强 PPP (6倍 Vmax)

# 解析结果
aa_wt = res_wt[:, 12]
aa_eng = res_eng[:, 12]
x_wt = res_wt[:, 13]
x_eng = res_eng[:, 13]
nadph_wt = res_wt[:, 15]
nadph_eng = res_eng[:, 15]
eth_wt = res_wt[:, 7]
eth_eng = res_eng[:, 7]

# ==========================================
# 绘制结果
# ==========================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 1. 青蒿酸产量
ax = axes[0, 0]
ax.plot(t, aa_wt, 'k--', label='Wild Type', lw=2)
ax.plot(t, aa_eng, 'r-', label='High PPP Flux', lw=2)
ax.set_title('Artemisinic Acid (AA) Production')
ax.set_ylabel('Concentration (mM)')
ax.legend()
ax.grid(alpha=0.3)

# 2. 细胞生长
ax = axes[0, 1]
ax.plot(t, x_wt, 'k--', label='Wild Type')
ax.plot(t, x_eng, 'g-', label='High PPP Flux')
ax.set_title('Cell Biomass Growth')
ax.set_ylabel('Biomass (AU)')
ax.legend()
ax.grid(alpha=0.3)

# 3. 胞内 NADPH 水平
ax = axes[1, 0]
ax.plot(t, nadph_wt, 'k--', label='Wild Type')
ax.plot(t, nadph_eng, 'b-', label='High PPP Flux')
ax.set_title('Intracellular NADPH Availability')
ax.set_ylabel('[NADPH] (mM)')
ax.legend()
ax.grid(alpha=0.3)

# 4. 乙醇积累 (副产物)
ax = axes[1, 1]
ax.plot(t, eth_wt, 'k--', label='Wild Type')
ax.plot(t, eth_eng, 'orange', label='High PPP Flux')
ax.set_title('Ethanol Accumulation (By-product)')
ax.set_ylabel('[Ethanol] (mM)')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.show()


