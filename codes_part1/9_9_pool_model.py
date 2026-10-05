import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ==============================
# 1. 常量与参数（同前）
# ==============================
gDW_per_CmolX = 24.0
m3_per_kg = 1e-3

# Enzyme parameters
qE1_max = 6.50e-2
delta0 = 5.50e-2
k1_1 = 0.10
kdE1 = 1.46e-2
qE2_max = 3.0e-2
kdE2 = 1.0e-2
qE3_max = 5.0e-2
kdE3 = 1.5e-2
qE4_max = 2.0e-2
kdE4 = 1.47e-2

# Kinetics
KS_1_1 = 9.80e-6
v1_2_max = 0.18
KGlyc_1_2 = 31.38
KAA_1_2 = 870.23
KATP_1_2 = 2.01
v1_3_max = 0.32
KGlyc_1_3 = 38.54
KAA_1_3 = 757.81
KSto_2_1 = 10.0
KSto_2_2 = 5.0
KPAA_3_1 = 0.1
CPAA_in = 0.003
KPAA_4_1 = 0.05

# External
CS_in = 8.30e-2
Fin = 0.0
M_broth0 = 4.0
kd_cell = 5e-3
kd_pen = 1.25e-3

EPS_X = 1e-8


def safe_log(x):
    return np.log(np.clip(x, 1e-12, None))


# ==============================
# 2. ODE System (with extra diagnostics)
# ==============================
def nine_pool_model(t, y):
    X_glyc, X_aa, X_atp, X_sto, X_paa, E1, E2, E3, E4, X_cell, S, P = y

    # 防负
    S = max(S, 0.0)
    X_cell = max(X_cell, EPS_X)

    # 转换为 mmol 单位
    X_glyc_mmol = X_glyc * 1e3
    X_aa_mmol = X_aa * 1e3
    X_sto_mmol = X_sto * 1e3
    X_paa_mmol = X_paa * 1e3

    # 生长速率
    term_glyc = X_glyc_mmol / (KGlyc_1_2 + X_glyc_mmol)
    term_aa = X_aa_mmol / (KAA_1_2 + X_aa_mmol)
    term_atp = X_atp / (KATP_1_2 + X_atp)
    mu = v1_2_max * term_glyc * term_aa * term_atp
    mu = max(mu, 0.0)  # 防负

    # 酶合成
    ratio = (mu + delta0) / k1_1
    qE1_syn = qE1_max * (ratio ** 5) / (1 + ratio ** 5)
    qE2_syn = qE2_max if mu > 1e-6 else 0.0
    qE3_syn = qE3_max * np.exp(-max(X_glyc_mmol, 0) / 50.0)
    qE4_syn = qE4_max if X_paa_mmol > 0.01 else 0.0

    dE1_dt = qE1_syn - kdE1 * E1 - mu * E1
    dE2_dt = qE2_syn - kdE2 * E2 - mu * E2
    dE3_dt = qE3_syn - kdE3 * E3 - mu * E3
    dE4_dt = qE4_syn - kdE4 * E4 - mu * E4

    # 反应速率
    v1_1 = E1 * safe_log(1.0 + S / (KS_1_1 + 1e-15))
    v1_2 = mu * X_cell
    v1_3 = v1_3_max * (X_glyc_mmol / (KGlyc_1_3 + X_glyc_mmol)) * (X_aa_mmol / (KAA_1_3 + X_aa_mmol))
    v2_1 = E2 * safe_log(1.0 + X_glyc_mmol / (KSto_2_1 + 1e-6))
    v2_2 = E2 * safe_log(1.0 + X_sto_mmol / (KSto_2_2 + 1e-6))
    v3_1 = E3 * safe_log(1.0 + X_paa_mmol / (KPAA_3_1 + 1e-6))
    v4_1 = E4 * safe_log(1.0 + CPAA_in / (KPAA_4_1 + 1e-6))
    v4_2 = kdE4 * E4 * X_paa_mmol

    # 内部池
    if X_cell < EPS_X:
        dX_glyc_dt = dX_aa_dt = dX_atp_dt = dX_sto_dt = dX_paa_dt = 0.0
    else:
        dX_glyc_dt = (v1_1 - v1_2 - v1_3 - v2_1 + v2_2) / X_cell
        dX_aa_dt = (v1_3 - v3_1) / X_cell
        dX_atp_dt = 0.0
        dX_sto_dt = (v2_1 - v2_2) / X_cell
        dX_paa_dt = (v4_1 - v4_2 - v3_1) / X_cell

    # 外部池 (batch)
    V = M_broth0 * m3_per_kg
    dX_cell_dt = mu * X_cell - kd_cell * X_cell
    dS_dt = -v1_1 * X_cell
    dP_dt = v3_1 * X_cell - kd_pen * P

    return [dX_glyc_dt, dX_aa_dt, dX_atp_dt, dX_sto_dt, dX_paa_dt,
            dE1_dt, dE2_dt, dE3_dt, dE4_dt,
            dX_cell_dt, dS_dt, dP_dt]


# ==============================
# 3. 事件函数：当系统进入稳态时终止
# ==============================
def steady_state_event(t, y):
    """
    当 biomass 和 penicillin 的变化率都很小时，认为达到稳态
    返回值为 0 时触发事件
    """
    X_cell = y[9]
    P = y[11]

    # 计算当前导数（调用模型）
    dydt = nine_pool_model(t, y)
    dXdt = dydt[9]
    dPdt = dydt[11]

    # 稳态条件：变化率低于阈值
    tol_X = 1e-6  # CmolX/m³/h
    tol_P = 1e-8  # mol/m³/h

    # 返回正值表示未稳态，0 表示稳态
    return min(abs(dXdt) - tol_X, abs(dPdt) - tol_P)


steady_state_event.terminal = True  # 触发时停止积分
steady_state_event.direction = -1  # 仅当从正变负时触发

# ==============================
# 4. 初始条件
# ==============================
y0 = [
    1e-3,  # X_glyc
    5e-2,  # X_aa
    3e-3,  # X_atp
    2e-2,  # X_sto
    1e-4,  # X_paa
    0.5,  # E1
    0.3,  # E2
    0.1,  # E3
    0.2,  # E4
    0.5,  # X_cell
    0.05,  # S
    0.0  # P
]

# ==============================
# 5. 求解（带事件）
# ==============================
t_span = (0, 200)  # 最大 200 小时
t_eval = np.linspace(0, 200, 1000)

sol = solve_ivp(
    nine_pool_model,
    t_span=t_span,
    y0=y0,
    t_eval=t_eval,
    method='LSODA',
    rtol=1e-5,
    atol=1e-8,
    events=steady_state_event,  # 关键！
    max_step=2.0
)

print("Simulation finished.")
if sol.success:
    print(f"Final time: {sol.t[-1]:.2f} h")
    if sol.t_events[0].size > 0:
        print("Terminated by steady-state event.")
else:
    print("Integration failed:", sol.message)

# ==============================
# 6. 绘图
# ==============================
plt.figure(figsize=(12, 10))

plt.subplot(3, 1, 1)
plt.plot(sol.t, sol.y[9], 'b-', linewidth=2, label='Biomass (CmolX/m³)')
plt.ylabel('Biomass')
plt.grid(True, alpha=0.3)
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(sol.t, sol.y[10], 'g-', linewidth=2, label='Glucose (mol/m³)')
plt.ylabel('Glucose')
plt.grid(True, alpha=0.3)
plt.legend()

plt.subplot(3, 1, 3)
plt.plot(sol.t, sol.y[11], 'r-', linewidth=2, label='Penicillin (mol/m³)')
plt.xlabel('Time (h)')
plt.ylabel('Penicillin')
plt.grid(True, alpha=0.3)
plt.legend()

plt.suptitle('9-Pool Model with Steady-State Termination', fontsize=14)
plt.tight_layout()
plt.show()
