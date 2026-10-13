# toy model 1
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 参数设置
k_cat1 = 0.18
K_M1 = 6
K_M2 = 3
k_cat2 = 0.04
K_M3 = 1.5

E1 = 6
E2 = 8

initial_conditions = [6, 4, 0, 0]
t_span = (0, 90)
t_eval = np.linspace(0, 90, 1000)

def kinetic_model(t, y):
    c1, c2, c3, c4 = y
    v1 = (k_cat1 * E1 * c1 * c2) / ((K_M1 + c1) * (K_M2 + c2))
    v2 = (k_cat2 * E2 * c3) / (K_M3 + c3)
    return [-v1, -v1, v1 - v2, v2]

sol = solve_ivp(kinetic_model, t_span, initial_conditions, t_eval=t_eval, rtol=1e-8, atol=1e-8)

time = sol.t
c1, c2, c3, c4 = sol.y

# 绘图（修复版）
plt.figure(figsize=(10, 6))
plt.plot(time, c1, 'b-', label='c1', linewidth=2)
plt.plot(time, c2, '-', color='orange', label='c2', linewidth=2)
plt.plot(time, c3, 'g-', label='c3', linewidth=2)
plt.plot(time, c4, 'r-', label='c4', linewidth=2)

plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Kinetic Simulation of Enzymatic Reaction Pathway')
plt.legend()
plt.grid(True, alpha=0.3)
plt.xlim(0, 90)
plt.ylim(0, 6.5)
plt.tight_layout()
plt.show()







# 糖酵解途径
# ----------------------------
# 1. Define kinetic parameters
# ----------------------------
# Vmax and Km for each enzyme (arbitrary units, for illustration)
params = {
    'Vmax_HK': 10.0,  # Hexokinase
    'Km_HK': 1.0,
    'Vmax_PFK': 8.0,  # Phosphofructokinase
    'Km_PFK': 0.5,
    'Vmax_GAPD': 12.0,  # Glyceraldehyde-3-P dehydrogenase (simplified)
    'Km_GAPD': 0.8,
    'Vmax_PGK': 12.0,  # Phosphoglycerate kinase
    'Km_PGK': 0.7,
    'Vmax_PYK': 10.0,  # Pyruvate kinase
    'Km_PYK': 1.2,
}

# Initial metabolite concentrations (mM)
y0 = [
    5.0,  # Glucose (Glc)
    0.1,  # Glucose-6-P (G6P)
    0.05,  # Fructose-1,6-BP (F16BP)
    0.1,  # 3-Phosphoglycerate (3PG)
    0.05,  # Phosphoenolpyruvate (PEP)
    0.01  # Pyruvate (Pyr)
]

metabolite_names = ['Glc', 'G6P', 'F16BP', '3PG', 'PEP', 'Pyr']


# ----------------------------
# 2. Define reaction rates (fluxes)
# ----------------------------
def reaction_rates(y, params):
    Glc, G6P, F16BP, _3PG, PEP, Pyr = y

    v_HK = params['Vmax_HK'] * Glc / (params['Km_HK'] + Glc)
    v_PFK = params['Vmax_PFK'] * G6P / (params['Km_PFK'] + G6P)
    v_GAPD = params['Vmax_GAPD'] * F16BP / (params['Km_GAPD'] + F16BP)
    v_PGK = params['Vmax_PGK'] * _3PG / (params['Km_PGK'] + _3PG)
    v_PYK = params['Vmax_PYK'] * PEP / (params['Km_PYK'] + PEP)

    return v_HK, v_PFK, v_GAPD, v_PGK, v_PYK


# ----------------------------
# 3. Define ODE system (mass balance)
# ----------------------------
def glycolysis_ode(t, y, params):
    #  t is unused because this is an autonomous system (自洽系统，rates depend only on y)
    v_HK, v_PFK, v_GAPD, v_PGK, v_PYK = reaction_rates(y, params)

    dGlc_dt = -v_HK
    dG6P_dt = v_HK - v_PFK
    dF16BP_dt = v_PFK - v_GAPD
    d3PG_dt = v_GAPD - v_PGK
    dPEP_dt = v_PGK - v_PYK
    dPyr_dt = v_PYK

    return [dGlc_dt, dG6P_dt, dF16BP_dt, d3PG_dt, dPEP_dt, dPyr_dt]


# ----------------------------
# 4. Simulate dynamics
# ----------------------------
t_span = (0, 12)  # Time from 0 to 10 (arbitrary units)
t_eval = np.linspace(0, 2, 100)

sol = solve_ivp(
    glycolysis_ode, # 微分方程的函数
    t_span,         # 积分时间区间
    y0,             # 初始条件
    args=(params,), # 额外参数
    t_eval=t_eval,  # 指定在哪些时间点输出解（否则只返回内部步长点）
    method='LSODA'  # Robust for stiff/non-stiff systems
)

# ----------------------------
# 5. Plot results
# ----------------------------
plt.figure(figsize=(10, 6))
for i, name in enumerate(metabolite_names):
    plt.plot(sol.t, sol.y[i], label=name, linewidth=2)

plt.xlabel('Time (a.u.)')
plt.ylabel('Concentration (mM)')
plt.title('Simplified Glycolysis Dynamics (Michaelis-Menten Kinetics)')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

# Optional: Print final pyruvate concentration
print(f"Final pyruvate concentration: {sol.y[-1][-1]:.3f} mM")







# 上述模型进行拓展 包括辅因子
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ----------------------------
# 1. Define kinetic parameters
# ----------------------------
params = {
    # Enzyme kinetics
    'Vmax_HK':   10.0,
    'Km_HK_Glc': 1.0,
    'Km_HK_ATP': 0.8,

    'Vmax_PFK':  8.0,
    'Km_PFK_G6P': 0.5,
    'Km_PFK_ATP': 0.6,

    'Vmax_GAPD': 12.0,
    'Km_GAPD_F16BP': 0.8,
    'Km_GAPD_NAD': 0.3,

    'Vmax_PGK':  12.0,
    'Km_PGK_3PG': 0.7,
    'Km_PGK_ADP': 0.9,

    'Vmax_PYK':  10.0,
    'Km_PYK_PEP': 1.2,
    'Km_PYK_ADP': 1.0,
}

# Total conserved pools (mM)
ATP_total = 3.0    # [ATP] + [ADP] = constant
NAD_total = 1.0    # [NAD+] + [NADH] = constant

# Initial metabolite concentrations (mM)
y0 = [
    5.0,   # Glc
    0.1,   # G6P
    0.05,  # F16BP
    0.1,   # 3PG
    0.05,  # PEP
    0.01,  # Pyr
    2.5,   # ATP  --> implies ADP = 0.5
    0.8,   # NAD+ --> implies NADH = 0.2
]

metabolite_names = ['Glc', 'G6P', 'F16BP', '3PG', 'PEP', 'Pyr', 'ATP', 'NAD+']

# ----------------------------
# 2. Reaction rates with cofactors
# ----------------------------
def reaction_rates(y, params):
    Glc, G6P, F16BP, _3PG, PEP, Pyr, ATP, NAD = y
    ADP = ATP_total - ATP
    NADH = NAD_total - NAD

    # HK: Glc + ATP → G6P + ADP
    v_HK = (params['Vmax_HK'] * Glc * ATP) / (
        (params['Km_HK_Glc'] + Glc) * (params['Km_HK_ATP'] + ATP)
    )

    # PFK: G6P + ATP → F16BP + ADP
    v_PFK = (params['Vmax_PFK'] * G6P * ATP) / (
        (params['Km_PFK_G6P'] + G6P) * (params['Km_PFK_ATP'] + ATP)
    )

    # GAPD: F16BP + NAD+ + Pi → 3PG + NADH + H+
    v_GAPD = (params['Vmax_GAPD'] * F16BP * NAD) / (
        (params['Km_GAPD_F16BP'] + F16BP) * (params['Km_GAPD_NAD'] + NAD)
    )

    # PGK: 3PG + ADP → PEP + ATP
    v_PGK = (params['Vmax_PGK'] * _3PG * ADP) / (
        (params['Km_PGK_3PG'] + _3PG) * (params['Km_PGK_ADP'] + ADP)
    )

    # PYK: PEP + ADP → Pyr + ATP
    v_PYK = (params['Vmax_PYK'] * PEP * ADP) / (
        (params['Km_PYK_PEP'] + PEP) * (params['Km_PYK_ADP'] + ADP)
    )

    return v_HK, v_PFK, v_GAPD, v_PGK, v_PYK, ADP, NADH

# ----------------------------
# 3. ODE system with cofactor balance
# ----------------------------
def glycolysis_ode(t, y, params):
    Glc, G6P, F16BP, _3PG, PEP, Pyr, ATP, NAD = y
    v_HK, v_PFK, v_GAPD, v_PGK, v_PYK, ADP, NADH = reaction_rates(y, params)

    dGlc_dt     = -v_HK
    dG6P_dt     = v_HK - v_PFK
    dF16BP_dt   = v_PFK - v_GAPD
    d3PG_dt     = v_GAPD - v_PGK
    dPEP_dt     = v_PGK - v_PYK
    dPyr_dt     = v_PYK

    # ATP balance: -v_HK - v_PFK + v_PGK + v_PYK
    dATP_dt     = -v_HK - v_PFK + v_PGK + v_PYK

    # NAD+ balance: -v_GAPD (consumed in GAPD reaction)
    dNAD_dt     = -v_GAPD

    return [dGlc_dt, dG6P_dt, dF16BP_dt, d3PG_dt, dPEP_dt, dPyr_dt, dATP_dt, dNAD_dt]

# ----------------------------
# 4. Simulate dynamics
# ----------------------------
t_span = (0, 15)
t_eval = np.linspace(0, 15, 200)

sol = solve_ivp(
    glycolysis_ode,
    t_span,
    y0,
    args=(params,),
    t_eval=t_eval,
    method='LSODA'
)

# Compute derived variables for plotting
ATP_vals = sol.y[6]
ADP_vals = ATP_total - ATP_vals
NAD_vals = sol.y[7]
NADH_vals = NAD_total - NAD_vals

# ----------------------------
# 5. Plot results
# ----------------------------
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Central metabolites
axes[0, 0].plot(sol.t, sol.y[:6].T)
axes[0, 0].set_title('Central Glycolytic Metabolites')
axes[0, 0].set_ylabel('Concentration (mM)')
axes[0, 0].legend(metabolite_names[:6])
axes[0, 0].grid(True, alpha=0.4)

# ATP/ADP
axes[0, 1].plot(sol.t, ATP_vals, label='ATP', linewidth=2)
axes[0, 1].plot(sol.t, ADP_vals, label='ADP', linewidth=2)
axes[0, 1].set_title('ATP/ADP Dynamics')
axes[0, 1].set_ylabel('Concentration (mM)')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.4)

# NAD+/NADH
axes[1, 0].plot(sol.t, NAD_vals, label='NAD+', linewidth=2)
axes[1, 0].plot(sol.t, NADH_vals, label='NADH', linewidth=2)
axes[1, 0].set_title('NAD+/NADH Dynamics')
axes[1, 0].set_xlabel('Time (a.u.)')
axes[1, 0].set_ylabel('Concentration (mM)')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.4)

# Net ATP production
net_ATP = ATP_vals - y0[6]  # change from initial
axes[1, 1].plot(sol.t, net_ATP, 'g-', linewidth=2)
axes[1, 1].set_title('Net ATP Produced')
axes[1, 1].set_xlabel('Time (a.u.)')
axes[1, 1].set_ylabel('Δ[ATP] (mM)')
axes[1, 1].grid(True, alpha=0.4)

plt.tight_layout()
plt.show()

# Final summary
print(f"Final pyruvate: {sol.y[5][-1]:.3f} mM")
print(f"Final ATP: {ATP_vals[-1]:.3f} mM (ADP: {ADP_vals[-1]:.3f} mM)")
print(f"Final NAD+: {NAD_vals[-1]:.3f} mM (NADH: {NADH_vals[-1]:.3f} mM)")
print(f"Net ATP yield: {net_ATP[-1]:.3f} mM")







# 将上述动力学模型耦合乳酸发酵
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ----------------------------
# 1. Kinetic parameters (including LDH)
# ----------------------------
params = {
    # Glycolysis enzymes
    'Vmax_HK':   10.0,
    'Km_HK_Glc': 1.0,
    'Km_HK_ATP': 0.8,

    'Vmax_PFK':  8.0,
    'Km_PFK_G6P': 0.5,
    'Km_PFK_ATP': 0.6,

    'Vmax_GAPD': 12.0,
    'Km_GAPD_F16BP': 0.8,
    'Km_GAPD_NAD': 0.3,

    'Vmax_PGK':  12.0,
    'Km_PGK_3PG': 0.7,
    'Km_PGK_ADP': 0.9,

    'Vmax_PYK':  10.0,
    'Km_PYK_PEP': 1.2,
    'Km_PYK_ADP': 1.0,

    # Lactate dehydrogenase (LDH)
    'Vmax_LDH':  6.0,      # Adjustable: higher = more lactate
    'Km_LDH_Pyr': 0.5,
    'Km_LDH_NADH': 0.2,
}

# Conserved pools (mM)
ATP_total = 3.0    # ATP + ADP
NAD_total = 1.0    # NAD+ + NADH

# Initial concentrations (mM)
y0 = [
    5.0,   # Glc
    0.1,   # G6P
    0.05,  # F16BP
    0.1,   # 3PG
    0.05,  # PEP
    0.01,  # Pyr
    0.0,   # Lac (lactate)
    2.5,   # ATP
    0.8,   # NAD+
]

metabolite_names = ['Glc', 'G6P', 'F16BP', '3PG', 'PEP', 'Pyr', 'Lac', 'ATP', 'NAD+']

# ----------------------------
# 2. Reaction rates with LDH
# ----------------------------
def reaction_rates(y, params):
    Glc, G6P, F16BP, _3PG, PEP, Pyr, Lac, ATP, NAD = y
    ADP = ATP_total - ATP
    NADH = NAD_total - NAD

    # Glycolysis fluxes
    v_HK = (params['Vmax_HK'] * Glc * ATP) / (
        (params['Km_HK_Glc'] + Glc) * (params['Km_HK_ATP'] + ATP)
    )
    v_PFK = (params['Vmax_PFK'] * G6P * ATP) / (
        (params['Km_PFK_G6P'] + G6P) * (params['Km_PFK_ATP'] + ATP)
    )
    v_GAPD = (params['Vmax_GAPD'] * F16BP * NAD) / (
        (params['Km_GAPD_F16BP'] + F16BP) * (params['Km_GAPD_NAD'] + NAD)
    )
    v_PGK = (params['Vmax_PGK'] * _3PG * ADP) / (
        (params['Km_PGK_3PG'] + _3PG) * (params['Km_PGK_ADP'] + ADP)
    )
    v_PYK = (params['Vmax_PYK'] * PEP * ADP) / (
        (params['Km_PYK_PEP'] + PEP) * (params['Km_PYK_ADP'] + ADP)
    )

    # LDH: Pyr + NADH → Lac + NAD+
    v_LDH = (params['Vmax_LDH'] * Pyr * NADH) / (
        (params['Km_LDH_Pyr'] + Pyr) * (params['Km_LDH_NADH'] + NADH)
    )

    return v_HK, v_PFK, v_GAPD, v_PGK, v_PYK, v_LDH, ADP, NADH

# ----------------------------
# 3. ODE system with lactate branch
# ----------------------------
def glycolysis_lactate_ode(t, y, params):
    Glc, G6P, F16BP, _3PG, PEP, Pyr, Lac, ATP, NAD = y
    v_HK, v_PFK, v_GAPD, v_PGK, v_PYK, v_LDH, ADP, NADH = reaction_rates(y, params)

    dGlc_dt     = -v_HK
    dG6P_dt     = v_HK - v_PFK
    dF16BP_dt   = v_PFK - v_GAPD
    d3PG_dt     = v_GAPD - v_PGK
    dPEP_dt     = v_PGK - v_PYK
    dPyr_dt     = v_PYK - v_LDH          # consumed by LDH
    dLac_dt     = v_LDH                  # produced from Pyr
    dATP_dt     = -v_HK - v_PFK + v_PGK + v_PYK  # unchanged
    dNAD_dt     = -v_GAPD + v_LDH        # GAPD consumes NAD+, LDH regenerates it

    return [dGlc_dt, dG6P_dt, dF16BP_dt, d3PG_dt, dPEP_dt, dPyr_dt, dLac_dt, dATP_dt, dNAD_dt]

# ----------------------------
# 4. Simulation
# ----------------------------
t_span = (0, 20)
t_eval = np.linspace(0, 20, 300)

sol = solve_ivp(
    glycolysis_lactate_ode,
    t_span,
    y0,
    args=(params,),
    t_eval=t_eval,
    method='LSODA'
)

# Derived variables
ATP_vals = sol.y[7]
ADP_vals = ATP_total - ATP_vals
NAD_vals = sol.y[8]
NADH_vals = NAD_total - NAD_vals
Lac_vals = sol.y[6]

# ----------------------------
# 5. Plotting
# ----------------------------
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Central metabolites (excluding Lac)
axes[0, 0].plot(sol.t, sol.y[:6].T)
axes[0, 0].set_title('Glycolytic Intermediates')
axes[0, 0].set_ylabel('Concentration (mM)')
axes[0, 0].legend(metabolite_names[:6])
axes[0, 0].grid(True, alpha=0.4)

# Lactate production
axes[0, 1].plot(sol.t, Lac_vals, 'm-', linewidth=2.5)
axes[0, 1].set_title('Lactate Accumulation')
axes[0, 1].set_ylabel('[Lac] (mM)')
axes[0, 1].grid(True, alpha=0.4)

# Pyruvate dynamics (now consumed)
axes[0, 2].plot(sol.t, sol.y[5], 'orange', linewidth=2)
axes[0, 2].set_title('Pyruvate Dynamics')
axes[0, 2].set_ylabel('[Pyr] (mM)')
axes[0, 2].grid(True, alpha=0.4)

# ATP/ADP
axes[1, 0].plot(sol.t, ATP_vals, label='ATP')
axes[1, 0].plot(sol.t, ADP_vals, label='ADP')
axes[1, 0].set_title('ATP/ADP')
axes[1, 0].set_xlabel('Time (a.u.)')
axes[1, 0].set_ylabel('Concentration (mM)')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.4)

# NAD+/NADH
axes[1, 1].plot(sol.t, NAD_vals, label='NAD+', linewidth=2)
axes[1, 1].plot(sol.t, NADH_vals, label='NADH', linewidth=2)
axes[1, 1].set_title('NAD+/NADH Redox Balance')
axes[1, 1].set_xlabel('Time (a.u.)')
axes[1, 1].set_ylabel('Concentration (mM)')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.4)

# Net ATP and lactate yield
axes[1, 2].plot(sol.t, ATP_vals - y0[7], 'g-', label='Net ATP')
axes[1, 2].plot(sol.t, Lac_vals, 'm--', label='Lactate')
axes[1, 2].set_title('Metabolic Output')
axes[1, 2].set_xlabel('Time (a.u.)')
axes[1, 2].set_ylabel('Concentration (mM)')
axes[1, 2].legend()
axes[1, 2].grid(True, alpha=0.4)

plt.tight_layout()
plt.show()

# ----------------------------
# 6. Summary
# ----------------------------
print(f"Final pyruvate: {sol.y[5][-1]:.4f} mM")
print(f"Final lactate:  {Lac_vals[-1]:.4f} mM")
print(f"Final ATP:      {ATP_vals[-1]:.4f} mM (ADP: {ADP_vals[-1]:.4f})")
print(f"Final NAD+:     {NAD_vals[-1]:.4f} mM (NADH: {NADH_vals[-1]:.4f})")
print(f"Net ATP yield:  {ATP_vals[-1] - y0[7]:.4f} mM")
print(f"Lactate/Pyruvate ratio: {Lac_vals[-1] / max(sol.y[5][-1], 1e-6):.2f}")





## 增加青蒿酸合成模块
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ----------------------------
# 1. Kinetic parameters
# ----------------------------
params = {
    # Glycolysis
    'Vmax_HK':   10.0,  'Km_HK_Glc': 1.0,  'Km_HK_ATP': 0.8,
    'Vmax_PFK':  8.0,   'Km_PFK_G6P': 0.5, 'Km_PFK_ATP': 0.6,
    'Vmax_GAPD': 12.0,  'Km_GAPD_F16BP': 0.8, 'Km_GAPD_NAD': 0.3,
    'Vmax_PGK':  12.0,  'Km_PGK_3PG': 0.7, 'Km_PGK_ADP': 0.9,
    'Vmax_PYK':  10.0,  'Km_PYK_PEP': 1.2, 'Km_PYK_ADP': 1.0,

    # Lactate fermentation
    'Vmax_LDH':  4.0,   'Km_LDH_Pyr': 0.5, 'Km_LDH_NADH': 0.2,

    # Pyruvate to Acetyl-CoA (simplified PDH-like)
    'Vmax_PDH':  2.0,   'Km_PDH_Pyr': 0.4, 'Km_PDH_NAD': 0.5,

    # Artemisinic acid pathway (lumped but step-wise)
    'Vmax_MVA':  0.8,   'Km_MVA_AcCoA': 0.1,  # AcCoA → FPP (consumes ATP/NADPH)
    'Vmax_ADS':  0.6,   'Km_ADS_FPP': 0.05,   # FPP → Amoradiene
    'Vmax_CYP':  0.5,   'Km_CYP_Amor': 0.02,  # Amor → Artemisinic acid (consumes NADPH)
}

# Conserved pools (mM)
ATP_total = 3.0      # ATP + ADP
NAD_total = 1.0      # NAD+ + NADH
NADP_total = 0.6     # NADP+ + NADPH

# Initial concentrations (mM)
y0 = [
    5.0,   # Glc
    0.1,   # G6P
    0.05,  # F16BP
    0.1,   # 3PG
    0.05,  # PEP
    0.01,  # Pyr
    0.0,   # Lac
    0.001, # AcCoA
    0.0,   # FPP
    0.0,   # Amoradiene (Amor)
    0.0,   # Artemisinic acid (AA)
    2.5,   # ATP
    0.8,   # NAD+
    0.4,   # NADPH
]

names = ['Glc','G6P','F16BP','3PG','PEP','Pyr','Lac','AcCoA','FPP','Amor','AA','ATP','NAD+','NADPH']

# ----------------------------
# 2. Reaction rates
# ----------------------------
def reaction_rates(y, params):
    (Glc, G6P, F16BP, _3PG, PEP, Pyr, Lac, AcCoA, FPP, Amor, AA,
     ATP, NAD, NADPH) = y

    ADP = ATP_total - ATP
    NADH = NAD_total - NAD
    NADP = NADP_total - NADPH

    # --- Glycolysis ---
    v_HK = (params['Vmax_HK'] * Glc * ATP) / ((params['Km_HK_Glc'] + Glc) * (params['Km_HK_ATP'] + ATP))
    v_PFK = (params['Vmax_PFK'] * G6P * ATP) / ((params['Km_PFK_G6P'] + G6P) * (params['Km_PFK_ATP'] + ATP))
    v_GAPD = (params['Vmax_GAPD'] * F16BP * NAD) / ((params['Km_GAPD_F16BP'] + F16BP) * (params['Km_GAPD_NAD'] + NAD))
    v_PGK = (params['Vmax_PGK'] * _3PG * ADP) / ((params['Km_PGK_3PG'] + _3PG) * (params['Km_PGK_ADP'] + ADP))
    v_PYK = (params['Vmax_PYK'] * PEP * ADP) / ((params['Km_PYK_PEP'] + PEP) * (params['Km_PYK_ADP'] + ADP))

    # --- Lactate ---
    v_LDH = (params['Vmax_LDH'] * Pyr * NADH) / ((params['Km_LDH_Pyr'] + Pyr) * (params['Km_LDH_NADH'] + NADH))

    # --- Acetyl-CoA synthesis ---
    v_PDH = (params['Vmax_PDH'] * Pyr * NAD) / ((params['Km_PDH_Pyr'] + Pyr) * (params['Km_PDH_NAD'] + NAD))

    # --- Artemisinic acid pathway ---
    v_MVA = (params['Vmax_MVA'] * AcCoA) / (params['Km_MVA_AcCoA'] + AcCoA)
    v_ADS = (params['Vmax_ADS'] * FPP) / (params['Km_ADS_FPP'] + FPP)
    v_CYP = (params['Vmax_CYP'] * Amor) / (params['Km_CYP_Amor'] + Amor)

    return (v_HK, v_PFK, v_GAPD, v_PGK, v_PYK, v_LDH, v_PDH,
            v_MVA, v_ADS, v_CYP, ADP, NADH, NADP)

# ----------------------------
# 3. ODE system
# ----------------------------
def full_model_ode(t, y, params):
    (Glc, G6P, F16BP, _3PG, PEP, Pyr, Lac, AcCoA, FPP, Amor, AA,
     ATP, NAD, NADPH) = y

    (v_HK, v_PFK, v_GAPD, v_PGK, v_PYK, v_LDH, v_PDH,
     v_MVA, v_ADS, v_CYP, ADP, NADH, NADP) = reaction_rates(y, params)

    dGlc_dt     = -v_HK
    dG6P_dt     = v_HK - v_PFK
    dF16BP_dt   = v_PFK - v_GAPD
    d3PG_dt     = v_GAPD - v_PGK
    dPEP_dt     = v_PGK - v_PYK
    dPyr_dt     = v_PYK - v_LDH - v_PDH          # consumed by LDH and PDH
    dLac_dt     = v_LDH
    dAcCoA_dt   = v_PDH - v_MVA
    dFPP_dt     = v_MVA - v_ADS
    dAmor_dt    = v_ADS - v_CYP
    dAA_dt      = v_CYP

    # ATP balance: glycolysis (+2 net) minus MVA cost (3 ATP per v_MVA)
    dATP_dt     = (-v_HK - v_PFK + v_PGK + v_PYK) - 3 * v_MVA

    # NAD+ balance: GAPD consumes, LDH/PDH produce/consume
    dNAD_dt     = -v_GAPD + v_LDH - v_PDH        # PDH consumes NAD+

    # NADPH balance: MVA (2 NADPH) + CYP (3 NADPH) consume
    dNADPH_dt   = -2 * v_MVA - 3 * v_CYP

    return [dGlc_dt, dG6P_dt, dF16BP_dt, d3PG_dt, dPEP_dt, dPyr_dt, dLac_dt,
            dAcCoA_dt, dFPP_dt, dAmor_dt, dAA_dt,
            dATP_dt, dNAD_dt, dNADPH_dt]

# ----------------------------
# 4. Simulation
# ----------------------------
t_span = (0, 30)
t_eval = np.linspace(0, 30, 400)

sol = solve_ivp(
    full_model_ode,
    t_span,
    y0,
    args=(params,),
    t_eval=t_eval,
    method='LSODA'
)

# Derived variables
ATP_vals = sol.y[11]
ADP_vals = ATP_total - ATP_vals
NAD_vals = sol.y[12]
NADH_vals = NAD_total - NAD_vals
NADPH_vals = sol.y[13]
AA_vals = sol.y[10]  # Artemisinic acid

# ----------------------------
# 5. Plotting
# ----------------------------
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Central carbon
axes[0, 0].plot(sol.t, sol.y[:6].T)
axes[0, 0].set_title('Central Metabolites')
axes[0, 0].legend(names[:6])
axes[0, 0].grid(True, alpha=0.4)

# Lactate & Pyr
axes[0, 1].plot(sol.t, sol.y[6], 'm-', label='Lactate')
axes[0, 1].plot(sol.t, sol.y[5], 'orange', label='Pyruvate')
axes[0, 1].set_title('Pyruvate Branching')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.4)

# AcCoA and AA pathway
axes[0, 2].plot(sol.t, sol.y[7], 'g--', label='AcCoA')
axes[0, 2].plot(sol.t, sol.y[8], 'b-.', label='FPP')
axes[0, 2].plot(sol.t, sol.y[9], 'c:', label='Amor')
axes[0, 2].plot(sol.t, AA_vals, 'r-', linewidth=2, label='Artemisinic Acid')
axes[0, 2].set_title('Artemisinic Acid Pathway')
axes[0, 2].legend()
axes[0, 2].grid(True, alpha=0.4)

# Energy charge
axes[1, 0].plot(sol.t, ATP_vals, label='ATP')
axes[1, 0].plot(sol.t, ADP_vals, label='ADP')
axes[1, 0].set_title('ATP/ADP')
axes[1, 0].set_xlabel('Time (h)')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.4)

# Redox cofactors
axes[1, 1].plot(sol.t, NAD_vals, label='NAD+')
axes[1, 1].plot(sol.t, NADH_vals, label='NADH')
axes[1, 1].plot(sol.t, NADPH_vals, label='NADPH', linewidth=2)
axes[1, 1].set_title('Redox Cofactors')
axes[1, 1].set_xlabel('Time (h)')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.4)

# Final product
axes[1, 2].plot(sol.t, AA_vals, 'r-', linewidth=2.5)
axes[1, 2].set_title('Artemisinic Acid Production')
axes[1, 2].set_xlabel('Time (h)')
axes[1, 2].set_ylabel('[AA] (mM)')
axes[1, 2].grid(True, alpha=0.4)

plt.tight_layout()
plt.show()

# ----------------------------
# 6. Summary
# ----------------------------
print(f"Final artemisinic acid: {AA_vals[-1]:.5f} mM")
print(f"Final ATP: {ATP_vals[-1]:.3f} mM")
print(f"Final NADPH: {NADPH_vals[-1]:.4f} mM")
print(f"Total lactate produced: {sol.y[6][-1]:.3f} mM")





# 将乳酸发酵改成为乙醇发酵
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ----------------------------
# 1. Kinetic parameters (replace LDH with PDC/ADH)
# ----------------------------
params = {
    # Glycolysis
    'Vmax_HK':   10.0,  'Km_HK_Glc': 1.0,  'Km_HK_ATP': 0.8,
    'Vmax_PFK':  8.0,   'Km_PFK_G6P': 0.5, 'Km_PFK_ATP': 0.6,
    'Vmax_GAPD': 12.0,  'Km_GAPD_F16BP': 0.8, 'Km_GAPD_NAD': 0.3,
    'Vmax_PGK':  12.0,  'Km_PGK_3PG': 0.7, 'Km_PGK_ADP': 0.9,
    'Vmax_PYK':  10.0,  'Km_PYK_PEP': 1.2, 'Km_PYK_ADP': 1.0,

    # Ethanol fermentation (REPLACES lactate)
    'Vmax_PDC':  6.0,   'Km_PDC_Pyr': 0.6,      # Pyruvate → Acetaldehyde
    'Vmax_ADH':  7.0,   'Km_ADH_AcAld': 0.1,    # Acetaldehyde + NADH → Ethanol
    'Km_ADH_NADH': 0.2,

    # Pyruvate to Acetyl-CoA (for AA pathway)
    'Vmax_PDH':  2.0,   'Km_PDH_Pyr': 0.4, 'Km_PDH_NAD': 0.5,

    # Artemisinic acid pathway
    'Vmax_MVA':  0.8,   'Km_MVA_AcCoA': 0.1,
    'Vmax_ADS':  0.6,   'Km_ADS_FPP': 0.05,
    'Vmax_CYP':  0.5,   'Km_CYP_Amor': 0.02,
}

# Conserved pools (mM)
ATP_total = 3.0
NAD_total = 1.0
NADP_total = 0.6

# Initial concentrations (mM)
y0 = [
    5.0,    # Glc
    0.1,    # G6P
    0.05,   # F16BP
    0.1,    # 3PG
    0.05,   # PEP
    0.01,   # Pyr
    0.0,    # Acetaldehyde (AcAld)
    0.0,    # Ethanol (EtOH)
    0.001,  # AcCoA
    0.0,    # FPP
    0.0,    # Amoradiene (Amor)
    0.0,    # Artemisinic acid (AA)
    2.5,    # ATP
    0.8,    # NAD+
    0.4,    # NADPH
]

names = ['Glc','G6P','F16BP','3PG','PEP','Pyr','AcAld','EtOH','AcCoA','FPP','Amor','AA','ATP','NAD+','NADPH']

# ----------------------------
# 2. Reaction rates (with PDC/ADH)
# ----------------------------
def reaction_rates(y, params):
    (Glc, G6P, F16BP, _3PG, PEP, Pyr, AcAld, EtOH, AcCoA, FPP, Amor, AA,
     ATP, NAD, NADPH) = y

    ADP = ATP_total - ATP
    NADH = NAD_total - NAD
    NADP = NADP_total - NADPH

    # Glycolysis
    v_HK = (params['Vmax_HK'] * Glc * ATP) / ((params['Km_HK_Glc'] + Glc) * (params['Km_HK_ATP'] + ATP))
    v_PFK = (params['Vmax_PFK'] * G6P * ATP) / ((params['Km_PFK_G6P'] + G6P) * (params['Km_PFK_ATP'] + ATP))
    v_GAPD = (params['Vmax_GAPD'] * F16BP * NAD) / ((params['Km_GAPD_F16BP'] + F16BP) * (params['Km_GAPD_NAD'] + NAD))
    v_PGK = (params['Vmax_PGK'] * _3PG * ADP) / ((params['Km_PGK_3PG'] + _3PG) * (params['Km_PGK_ADP'] + ADP))
    v_PYK = (params['Vmax_PYK'] * PEP * ADP) / ((params['Km_PYK_PEP'] + PEP) * (params['Km_PYK_ADP'] + ADP))

    # Ethanol fermentation
    v_PDC = (params['Vmax_PDC'] * Pyr) / (params['Km_PDC_Pyr'] + Pyr)
    v_ADH = (params['Vmax_ADH'] * AcAld * NADH) / (
        (params['Km_ADH_AcAld'] + AcAld) * (params['Km_ADH_NADH'] + NADH)
    )

    # Acetyl-CoA for AA pathway
    v_PDH = (params['Vmax_PDH'] * Pyr * NAD) / ((params['Km_PDH_Pyr'] + Pyr) * (params['Km_PDH_NAD'] + NAD))

    # Artemisinic acid pathway
    v_MVA = (params['Vmax_MVA'] * AcCoA) / (params['Km_MVA_AcCoA'] + AcCoA)
    v_ADS = (params['Vmax_ADS'] * FPP) / (params['Km_ADS_FPP'] + FPP)
    v_CYP = (params['Vmax_CYP'] * Amor) / (params['Km_CYP_Amor'] + Amor)

    return (v_HK, v_PFK, v_GAPD, v_PGK, v_PYK, v_PDC, v_ADH, v_PDH,
            v_MVA, v_ADS, v_CYP, ADP, NADH, NADP)

# ----------------------------
# 3. ODE system
# ----------------------------
def ethanol_aa_model_ode(t, y, params):
    (Glc, G6P, F16BP, _3PG, PEP, Pyr, AcAld, EtOH, AcCoA, FPP, Amor, AA,
     ATP, NAD, NADPH) = y

    (v_HK, v_PFK, v_GAPD, v_PGK, v_PYK, v_PDC, v_ADH, v_PDH,
     v_MVA, v_ADS, v_CYP, ADP, NADH, NADP) = reaction_rates(y, params)

    dGlc_dt     = -v_HK
    dG6P_dt     = v_HK - v_PFK
    dF16BP_dt   = v_PFK - v_GAPD
    d3PG_dt     = v_GAPD - v_PGK
    dPEP_dt     = v_PGK - v_PYK
    dPyr_dt     = v_PYK - v_PDC - v_PDH          # consumed by PDC and PDH
    dAcAld_dt   = v_PDC - v_ADH                  # produced by PDC, consumed by ADH
    dEtOH_dt    = v_ADH                          # final ethanol product

    dAcCoA_dt   = v_PDH - v_MVA
    dFPP_dt     = v_MVA - v_ADS
    dAmor_dt    = v_ADS - v_CYP
    dAA_dt      = v_CYP

    dATP_dt     = (-v_HK - v_PFK + v_PGK + v_PYK) - 3 * v_MVA
    dNAD_dt     = -v_GAPD + v_ADH - v_PDH        # ADH regenerates NAD+, PDH consumes it
    dNADPH_dt   = -2 * v_MVA - 3 * v_CYP

    return [dGlc_dt, dG6P_dt, dF16BP_dt, d3PG_dt, dPEP_dt, dPyr_dt,
            dAcAld_dt, dEtOH_dt,
            dAcCoA_dt, dFPP_dt, dAmor_dt, dAA_dt,
            dATP_dt, dNAD_dt, dNADPH_dt]

# ----------------------------
# 4. Simulation
# ----------------------------
t_span = (0, 30)
t_eval = np.linspace(0, 30, 400)

sol = solve_ivp(
    ethanol_aa_model_ode,
    t_span,
    y0,
    args=(params,),
    t_eval=t_eval,
    method='LSODA'
)

# Derived variables
ATP_vals = sol.y[12]
ADP_vals = ATP_total - ATP_vals
NAD_vals = sol.y[13]
NADH_vals = NAD_total - NAD_vals
NADPH_vals = sol.y[14]
AA_vals = sol.y[11]      # Artemisinic acid
EtOH_vals = sol.y[7]     # Ethanol

# ----------------------------
# 5. Plotting
# ----------------------------
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Central metabolites
axes[0, 0].plot(sol.t, sol.y[:6].T)
axes[0, 0].set_title('Central Glycolytic Metabolites')
axes[0, 0].legend(names[:6])
axes[0, 0].grid(True, alpha=0.4)

# Ethanol pathway
axes[0, 1].plot(sol.t, sol.y[5], 'orange', label='Pyruvate')
axes[0, 1].plot(sol.t, sol.y[6], 'purple', label='Acetaldehyde')
axes[0, 1].plot(sol.t, EtOH_vals, 'brown', linewidth=2, label='Ethanol')
axes[0, 1].set_title('Ethanol Fermentation')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.4)

# AA pathway
axes[0, 2].plot(sol.t, sol.y[8], 'g--', label='AcCoA')
axes[0, 2].plot(sol.t, sol.y[9], 'b-.', label='FPP')
axes[0, 2].plot(sol.t, sol.y[10], 'c:', label='Amor')
axes[0, 2].plot(sol.t, AA_vals, 'r-', linewidth=2, label='Artemisinic Acid')
axes[0, 2].set_title('Artemisinic Acid Synthesis')
axes[0, 2].legend()
axes[0, 2].grid(True, alpha=0.4)

# Energy
axes[1, 0].plot(sol.t, ATP_vals, label='ATP')
axes[1, 0].plot(sol.t, ADP_vals, label='ADP')
axes[1, 0].set_title('ATP/ADP Dynamics')
axes[1, 0].set_xlabel('Time (h)')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.4)

# Redox balance
axes[1, 1].plot(sol.t, NAD_vals, label='NAD+')
axes[1, 1].plot(sol.t, NADH_vals, label='NADH')
axes[1, 1].plot(sol.t, NADPH_vals, label='NADPH', linewidth=2)
axes[1, 1].set_title('Redox Cofactors')
axes[1, 1].set_xlabel('Time (h)')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.4)

# Final products
axes[1, 2].plot(sol.t, EtOH_vals, 'brown', label='Ethanol', linewidth=2)
axes[1, 2].plot(sol.t, AA_vals, 'r--', label='Artemisinic Acid', linewidth=2)
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
print(f"Acetaldehyde peak:      {np.max(sol.y[6]):.4f} mM (transient intermediate)")




# 继续增加氧化磷酸化与TCA途径
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ----------------------------
# 1. Kinetic parameters (add TCA/OXPHOS)
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
}

# Conserved pools (mM)
ATP_total = 5.0      # increased due to respiration
NAD_total = 1.5      # higher total redox capacity
NADP_total = 0.6
FAD_total = 0.2      # FAD + FADH2

# Assume O2 is abundant (not limiting), so not a state variable

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
# 2. Reaction rates
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

    return (v_HK, v_PFK, v_GAPD, v_PGK, v_PYK, v_PDC, v_ADH,
            v_AcCoA_cyt, v_PDH_mito, v_TCA,
            v_OXPHOS_NADH, v_OXPHOS_FADH2,
            v_MVA, v_ADS, v_CYP,
            ADP, NADH, NADP, FADH2)

# ----------------------------
# 3. ODE system with TCA + OXPHOS
# ----------------------------
def full_respiratory_model_ode(t, y, params):
    (Glc, G6P, F16BP, _3PG, PEP, Pyr, AcAld, EtOH,
     AcCoA_cyt, FPP, Amor, AA, AcCoA_mito,
     ATP, NAD, NADPH, FAD) = y

    (v_HK, v_PFK, v_GAPD, v_PGK, v_PYK, v_PDC, v_ADH,
     v_AcCoA_cyt, v_PDH_mito, v_TCA,
     v_OXPHOS_NADH, v_OXPHOS_FADH2,
     v_MVA, v_ADS, v_CYP,
     ADP, NADH, NADP, FADH2) = reaction_rates(y, params)

    dGlc_dt       = -v_HK
    dG6P_dt       = v_HK - v_PFK
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
    #   Glycolysis: -v_HK -v_PFK + v_PGK + v_PYK
    #   TCA: + v_TCA (GTP ≈ ATP)
    #   OXPHOS: + 2.5 * v_OXPHOS_NADH + 1.5 * v_OXPHOS_FADH2
    #   AA pathway: -3 * v_MVA
    dATP_dt = (-v_HK - v_PFK + v_PGK + v_PYK) \
              + v_TCA \
              + 2.5 * v_OXPHOS_NADH + 1.5 * v_OXPHOS_FADH2 \
              - 3 * v_MVA

    # NAD+ balance:
    #   GAPD: -v_GAPD
    #   ADH: +v_ADH (regenerates)
    #   Cytosolic AcCoA & PDH_mito: -v_AcCoA_cyt - v_PDH_mito
    #   OXPHOS: +v_OXPHOS_NADH (oxidizes NADH → NAD+)
    dNAD_dt = -v_GAPD + v_ADH - v_AcCoA_cyt - v_PDH_mito + v_OXPHOS_NADH

    # NADPH balance (unchanged)
    dNADPH_dt = -2 * v_MVA - 3 * v_CYP

    # FAD balance: TCA produces FADH2, OXPHOS consumes it
    dFAD_dt = -v_TCA + v_OXPHOS_FADH2  # FAD consumed in TCA, regenerated in OXPHOS

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
axes[0, 2].set_title('Artemisinic Acid Synthesis')
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
axes[1, 1].set_title('Redox Balance (Respiration oxidizes NADH)')
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
print(f"Final ATP:              {ATP_vals[-1]:.3f} mM (vs ~2.8 without respiration)")
print(f"Final NAD+:             {NAD_vals[-1]:.3f} mM")
print(f"mito-AcCoA peak:        {np.max(AcCoA_mito_vals):.4f} mM")




