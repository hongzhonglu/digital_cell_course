from cobra import Model, Reaction, Metabolite
import pandas as pd
import matplotlib.pyplot as plt


def build_minimal_ecmodel():
    """
    构建一个极简的酶约束代谢模型，只包含3个核心反应：
    1. 糖酵解：葡萄糖 → 丙酮酸
    2. 乳酸合成：丙酮酸 → 乳酸
    3. ATP维持：ATP → ADP

    模型结构模仿GECKO：每个代谢反应消耗特定酶，总酶量有限。
    """
    model = Model("Minimal_EMP_ecModel")

    # ========== 1. 定义代谢物 ==========
    metabolites = {
        'glc': Metabolite('glc_c', name='Glucose', compartment='c'),
        'pyr': Metabolite('pyr_c', name='Pyruvate', compartment='c'),
        'lac': Metabolite('lac_c', name='Lactate', compartment='c'),
        'atp': Metabolite('atp_c', name='ATP', compartment='c'),
        'adp': Metabolite('adp_c', name='ADP', compartment='c'),
        'nad': Metabolite('nad_c', name='NAD+', compartment='c'),
        'nadh': Metabolite('nadh_c', name='NADH', compartment='c'),
        # 酶代谢物（关键！这是GECKO的核心）
        'E_EMP': Metabolite('E_EMP_c', name='Glycolysis_Enzyme', compartment='c'),
        'E_LDH': Metabolite('E_LDH_c', name='LDH_Enzyme', compartment='c'),
        'E_ATPM': Metabolite('E_ATPM_c', name='ATPase_Enzyme', compartment='c'),
        'E_total': Metabolite('E_total_c', name='Total_Enzyme_Pool', compartment='c')
    }
    model.add_metabolites(list(metabolites.values()))

    # ========== 2. 定义代谢反应 ==========
    reactions = []

    # 反应1: 糖酵解 (简化：1葡萄糖 → 2丙酮酸 + 2ATP + 2NADH)
    r_emp = Reaction('R_Glycolysis')
    r_emp.name = 'Glycolysis'
    r_emp.add_metabolites({
        metabolites['glc']: -1.0,
        metabolites['adp']: -2.0,
        metabolites['nad']: -2.0,
        metabolites['pyr']: 2.0,
        metabolites['atp']: 2.0,
        metabolites['nadh']: 2.0
    })
    reactions.append(r_emp)

    # 反应2: 乳酸脱氢酶 (丙酮酸 + NADH → 乳酸 + NAD+)
    r_ldh = Reaction('R_LDH')
    r_ldh.name = 'Lactate Dehydrogenase'
    r_ldh.add_metabolites({
        metabolites['pyr']: -1.0,
        metabolites['nadh']: -1.0,
        metabolites['lac']: 1.0,
        metabolites['nad']: 1.0
    })
    reactions.append(r_ldh)

    # 反应3: ATP维持消耗 (代表基础代谢)
    r_atpm = Reaction('R_ATPM')
    r_atpm.name = 'ATP Maintenance'
    r_atpm.add_metabolites({
        metabolites['atp']: -1.0,
        metabolites['adp']: 1.0
    })
    r_atpm.lower_bound = 0.5  # 最小ATP维持需求
    reactions.append(r_atpm)

    # 交换反应：葡萄糖摄取
    r_ex_glc = Reaction('EX_glc')
    r_ex_glc.name = 'Glucose Uptake'
    r_ex_glc.add_metabolites({metabolites['glc']: -1.0})
    r_ex_glc.lower_bound = -10.0  # 最大摄取10 mmol/gDW/h
    r_ex_glc.upper_bound = 0.0
    reactions.append(r_ex_glc)

    # 交换反应：乳酸分泌
    r_ex_lac = Reaction('EX_lac')
    r_ex_lac.name = 'Lactate Secretion'
    r_ex_lac.add_metabolites({metabolites['lac']: -1.0})
    r_ex_lac.lower_bound = 0.0
    reactions.append(r_ex_lac)

    model.add_reactions(reactions)

    # ========== 3. GECKO风格酶约束 ==========
    # 每个酶代谢物都有独立的合成/降解反应
    enzyme_reactions = []

    # 糖酵解酶的合成/降解
    r_E_EMP = Reaction('SYN_E_EMP')
    r_E_EMP.name = 'Synthesis of Glycolysis Enzymes'
    r_E_EMP.add_metabolites({
        metabolites['E_total']: -1.0,  # 消耗总酶池
        metabolites['E_EMP']: 1.0  # 合成糖酵解酶
    })
    r_E_EMP.lower_bound = 0.0
    r_E_EMP.upper_bound = 0.3  # 最多分配30%总酶给糖酵解
    enzyme_reactions.append(r_E_EMP)

    # LDH酶的合成/降解
    r_E_LDH = Reaction('SYN_E_LDH')
    r_E_LDH.name = 'Synthesis of LDH Enzymes'
    r_E_LDH.add_metabolites({
        metabolites['E_total']: -1.0,
        metabolites['E_LDH']: 1.0
    })
    r_E_LDH.lower_bound = 0.0
    r_E_LDH.upper_bound = 0.3  # 最多分配30%总酶给LDH
    enzyme_reactions.append(r_E_LDH)

    # ATP酶（维持）的合成/降解
    r_E_ATPM = Reaction('SYN_E_ATPM')
    r_E_ATPM.name = 'Synthesis of ATPase Enzymes'
    r_E_ATPM.add_metabolites({
        metabolites['E_total']: -1.0,
        metabolites['E_ATPM']: 1.0
    })
    r_E_ATPM.lower_bound = 0.0
    r_E_ATPM.upper_bound = 0.2  # 最多分配20%总酶给ATP维持
    enzyme_reactions.append(r_E_ATPM)

    # 总酶池的供应（关键约束）
    r_E_total = Reaction('SYN_E_total')
    r_E_total.name = 'Total Enzyme Synthesis'
    r_E_total.add_metabolites({metabolites['E_total']: 1.0})
    r_E_total.lower_bound = 0.0
    r_E_total.upper_bound = 1.0  # 总酶池上限：1 g/gDW
    enzyme_reactions.append(r_E_total)

    model.add_reactions(enzyme_reactions)

    # ========== 4. 连接酶与代谢反应 ==========
    # 关键：代谢反应消耗酶，消耗量 = 反应通量 / kcat

    # 糖酵解反应需要消耗E_EMP酶
    kcat_EMP = 50.0  # 单位：1/h，即每g酶每小时催化50 mmol底物
    r_emp.add_metabolites({metabolites['E_EMP']: -(1.0 / kcat_EMP)})

    # LDH反应需要消耗E_LDH酶
    kcat_LDH = 80.0
    r_ldh.add_metabolites({metabolites['E_LDH']: -(1.0 / kcat_LDH)})

    # ATP维持反应需要消耗E_ATPM酶
    kcat_ATPM = 30.0
    r_atpm.add_metabolites({metabolites['E_ATPM']: -(1.0 / kcat_ATPM)})

    # ========== 5. 设置目标函数 ==========
    # 最大化乳酸分泌
    model.objective = 'EX_lac'

    return model, {'kcat_EMP': kcat_EMP, 'kcat_LDH': kcat_LDH, 'kcat_ATPM': kcat_ATPM}


def predict_lactate_production(model):  # update
    """
    使用酶约束模型预测乳酸合成
    """
    print("=" * 60)
    print("极简酶约束代谢模型 - 乳酸合成预测")
    print("=" * 60)

    # 场景1: 最大乳酸合成
    print("\n1. 最大乳酸合成潜力:")
    model.objective = 'EX_lac'
    solution_max_lac = model.optimize()

    if solution_max_lac.status == 'optimal':
        lac_rate = solution_max_lac.fluxes['EX_lac']
        glc_uptake = -solution_max_lac.fluxes['EX_glc']
        atpm = solution_max_lac.fluxes['R_ATPM']

        print(f"   最大乳酸分泌速率: {lac_rate:.3f} mmol/gDW/h")
        print(f"   所需葡萄糖摄取: {glc_uptake:.3f} mmol/gDW/h")
        print(f"   ATP维持通量: {atpm:.3f} mmol/gDW/h")
        print(f"   乳酸得率: {lac_rate / glc_uptake:.3f} mol乳酸/mol葡萄糖")

        # 酶使用情况
        print(f"\n   酶使用情况:")
        print(f"   糖酵解酶 (E_EMP): {solution_max_lac.fluxes['SYN_E_EMP']:.3f} g/gDW")
        print(f"   LDH酶 (E_LDH): {solution_max_lac.fluxes['SYN_E_LDH']:.3f} g/gDW")
        print(f"   ATP酶 (E_ATPM): {solution_max_lac.fluxes['SYN_E_ATPM']:.3f} g/gDW")
        print(f"   总酶使用量: {solution_max_lac.fluxes['SYN_E_total']:.3f} g/gDW")
    else:
        print(f"   优化失败: {solution_max_lac.status}")

    # 场景2: 不同总酶量下的乳酸合成
    print("\n2. 总酶量对乳酸合成的影响:")

    total_enzyme_levels = [0.1, 0.3, 0.5, 0.8, 1.0, 1.5]
    results = []

    for enzyme_level in total_enzyme_levels:
        # 修改总酶池上限
        model.reactions.get_by_id('SYN_E_total').upper_bound = enzyme_level

        # 优化乳酸合成
        solution = model.optimize()

        if solution.status == 'optimal':
            lac_rate = solution.fluxes['EX_lac']
            glc_uptake = -solution.fluxes['EX_glc']
            results.append({
                'Total_Enzyme': enzyme_level,
                'Lactate_Rate': lac_rate,
                'Glucose_Uptake': glc_uptake,
                'Yield': lac_rate / glc_uptake if glc_uptake > 0 else 0
            })

    # 恢复原始总酶量
    model.reactions.get_by_id('SYN_E_total').upper_bound = 1.0

    # 显示结果表格
    df_results = pd.DataFrame(results)
    print("\n   总酶量(g/gDW) | 乳酸速率(mmol/h) | 葡萄糖摄取(mmol/h) | 得率")
    print("   " + "-" * 55)
    for _, row in df_results.iterrows():
        print(
            f"   {row['Total_Enzyme']:12.1f} | {row['Lactate_Rate']:16.3f} | {row['Glucose_Uptake']:17.3f} | {row['Yield']:.3f}")

    # 场景3: 酶分配策略的影响
    print("\n3. 酶分配策略对乳酸合成的影响:")

    # 固定总酶量为1.0 g/gDW
    model.reactions.get_by_id('SYN_E_total').upper_bound = 1.0

    # 测试不同的酶分配比例
    strategies = [
        ('高EMP', 0.5, 0.3, 0.2),  # 50% EMP, 30% LDH, 20% ATPM
        ('高LDH', 0.3, 0.5, 0.2),  # 30% EMP, 50% LDH, 20% ATPM
        ('平衡', 0.4, 0.4, 0.2),  # 40% EMP, 40% LDH, 20% ATPM
    ]

    print("\n   策略 | EMP上限 | LDH上限 | ATPM上限 | 乳酸速率")
    print("   " + "-" * 50)

    for name, emp_max, ldh_max, atpm_max in strategies:
        # 设置各酶的上限
        model.reactions.get_by_id('SYN_E_EMP').upper_bound = emp_max
        model.reactions.get_by_id('SYN_E_LDH').upper_bound = ldh_max
        model.reactions.get_by_id('SYN_E_ATPM').upper_bound = atpm_max

        solution = model.optimize()

        if solution.status == 'optimal':
            lac_rate = solution.fluxes['EX_lac']
            print(f"   {name:6} | {emp_max:7.1f} | {ldh_max:7.1f} | {atpm_max:8.1f} | {lac_rate:10.3f}")

    # 恢复原始酶分配上限
    model.reactions.get_by_id('SYN_E_EMP').upper_bound = 0.3
    model.reactions.get_by_id('SYN_E_LDH').upper_bound = 0.3
    model.reactions.get_by_id('SYN_E_ATPM').upper_bound = 0.2

    return df_results


def visualize_results(df_results):
    """
    可视化模拟结果
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # 图1: 总酶量 vs 乳酸合成速率
    ax1 = axes[0]
    ax1.plot(df_results['Total_Enzyme'], df_results['Lactate_Rate'], 'bo-', linewidth=2)
    ax1.set_xlabel('Total Enzyme (g/gDW)')
    ax1.set_ylabel('Lactate Production Rate (mmol/gDW/h)')
    ax1.set_title('Effect of Total Enzyme on Lactate Production')
    ax1.grid(True, alpha=0.3)

    # 图2: 总酶量 vs 得率
    ax2 = axes[1]
    ax2.plot(df_results['Total_Enzyme'], df_results['Yield'], 'ro-', linewidth=2)
    ax2.set_xlabel('Total Enzyme (g/gDW)')
    ax2.set_ylabel('Lactate Yield (mol lactate/mol glucose)')
    ax2.set_title('Effect of Total Enzyme on Lactate Yield')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


def analyze_enzyme_bottlenecks(model):
    """
    分析酶瓶颈：如果提高特定酶的kcat，乳酸合成能增加多少？
    """
    print("\n" + "=" * 60)
    print("酶瓶颈分析")
    print("=" * 60)

    # 基础情况
    model.objective = 'EX_lac'
    base_solution = model.optimize()
    base_lac = base_solution.fluxes['EX_lac']

    print(f"\n基础情况乳酸合成: {base_lac:.3f} mmol/gDW/h")

    # 测试提高各酶的kcat
    enzymes = [
        ('糖酵解酶 (E_EMP)', 'R_Glycolysis', 2.0),  # kcat提高2倍
        ('LDH酶 (E_LDH)', 'R_LDH', 2.0),  # kcat提高2倍
        ('ATP酶 (E_ATPM)', 'R_ATPM', 2.0),  # kcat提高2倍
    ]

    print("\n提高kcat对乳酸合成的影响:")
    print("   酶类型          | 提高倍数 | 乳酸速率 | 提升比例")
    print("   " + "-" * 50)

    for enzyme_name, reaction_id, multiplier in enzymes:
        # 获取原始代谢反应和kcat系数
        reaction = model.reactions.get_by_id(reaction_id)

        # 找出该反应中酶代谢物的系数（原始kcat的倒数）
        for met, coeff in reaction.metabolites.items():
            if 'E_' in met.id and reaction_id in ['R_Glycolysis', 'R_LDH', 'R_ATPM']:
                # 计算新的kcat（提高multiplier倍）
                new_kcat = 1.0 / abs(coeff) * multiplier
                new_coeff = -1.0 / new_kcat

                # 临时修改酶系数
                original_coeff = coeff
                reaction.add_metabolites({met: -coeff})  # 移除原系数
                reaction.add_metabolites({met: new_coeff})  # 添加新系数

                # 重新优化
                solution = model.optimize()
                new_lac = solution.fluxes['EX_lac']
                improvement = (new_lac - base_lac) / base_lac * 100

                print(f"   {enzyme_name:16} | {multiplier:8.1f} | {new_lac:8.3f} | {improvement:7.1f}%")

                # 恢复原始系数
                reaction.add_metabolites({met: -new_coeff})
                reaction.add_metabolites({met: original_coeff})

                break


# ===== 主程序 =====
if __name__ == "__main__":
    # 1. 构建极简酶约束模型
    print("构建极简酶约束代谢模型...")
    ecmodel, kcat_values = build_minimal_ecmodel()

    print(f"\n模型包含:")
    print(f"  代谢物: {len(ecmodel.metabolites)} 个")
    print(f"  反应: {len(ecmodel.reactions)} 个")
    print(f"  kcat值: EMP={kcat_values['kcat_EMP']}, LDH={kcat_values['kcat_LDH']}, ATPM={kcat_values['kcat_ATPM']}")

    # 2. 预测乳酸合成
    results_df = predict_lactate_production(ecmodel) # update

    # 3. 可视化结果
    visualize_results(results_df)

    # 4. 分析酶瓶颈
    analyze_enzyme_bottlenecks(ecmodel)

    print("\n" + "=" * 60)
    print("模型构建与预测完成！")
    print("=" * 60)