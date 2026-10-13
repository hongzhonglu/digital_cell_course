from cobra import Model, Reaction, Metabolite


def build_corrected_model():
    # 初始化模型
    model = Model('Glycolysis_Lactate_Corrected')

    # --- 1. 定义代谢物 ---
    # 核心碳架
    glc_c = Metabolite('glc_c', name='Glucose', formula='C6H12O6', compartment='c')
    g6p_c = Metabolite('g6p_c', name='Glucose-6-Phosphate', compartment='c')
    f6p_c = Metabolite('f6p_c', name='Fructose-6-Phosphate', compartment='c')
    fdp_c = Metabolite('fdp_c', name='Fructose-1,6-Bisphosphate', compartment='c')
    dhap_c = Metabolite('dhap_c', name='Dihydroxyacetone phosphate', compartment='c')
    g3p_c = Metabolite('g3p_c', name='Glyceraldehyde-3-Phosphate', compartment='c')
    dpg13_c = Metabolite('13dpg_c', name='1,3-Bisphosphoglycerate', compartment='c')
    pg3_c = Metabolite('3pg_c', name='3-Phosphoglycerate', compartment='c')
    pg2_c = Metabolite('2pg_c', name='2-Phosphoglycerate', compartment='c')
    pep_c = Metabolite('pep_c', name='Phosphoenolpyruvate', compartment='c')
    pyr_c = Metabolite('pyr_c', name='Pyruvate', compartment='c')
    lac_c = Metabolite('lac_c', name='Lactate', formula='C3H5O3', compartment='c')  # 乳酸根

    # 辅因子
    atp_c = Metabolite('atp_c', name='ATP', compartment='c')
    adp_c = Metabolite('adp_c', name='ADP', compartment='c')
    nadh_c = Metabolite('nadh_c', name='NADH', compartment='c')
    nad_c = Metabolite('nad_c', name='NAD+', compartment='c')
    pi_c = Metabolite('pi_c', name='Phosphate', compartment='c')
    h_c = Metabolite('h_c', name='H+', compartment='c')
    h2o_c = Metabolite('h2o_c', name='Water', formula='H2O', compartment='c')

    # --- 2. 定义内部反应 (糖酵解 + 乳酸发酵) ---

    # HEX1: Glc + ATP -> G6P + ADP + H
    hex1 = Reaction('HEX1')
    hex1.add_metabolites({glc_c: -1, atp_c: -1, g6p_c: 1, adp_c: 1, h_c: 1})

    # PGI: G6P <-> F6P
    pgi = Reaction('PGI')
    pgi.add_metabolites({g6p_c: -1, f6p_c: 1})
    pgi.lower_bound = -1000

    # PFK: F6P + ATP -> FDP + ADP + H
    pfk = Reaction('PFK')
    pfk.add_metabolites({f6p_c: -1, atp_c: -1, fdp_c: 1, adp_c: 1, h_c: 1})

    # FBA: FDP <-> G3P + DHAP
    fba = Reaction('FBA')
    fba.add_metabolites({fdp_c: -1, g3p_c: 1, dhap_c: 1})
    fba.lower_bound = -1000

    # TPI: DHAP <-> G3P
    tpi = Reaction('TPI')
    tpi.add_metabolites({dhap_c: -1, g3p_c: 1})
    tpi.lower_bound = -1000

    # GAPD: G3P + NAD + Pi <-> 13DPG + NADH + H
    gapd = Reaction('GAPD')
    gapd.add_metabolites({g3p_c: -1, nad_c: -1, pi_c: -1, dpg13_c: 1, nadh_c: 1, h_c: 1})
    gapd.lower_bound = -1000

    # PGK: 13DPG + ADP <-> 3PG + ATP
    pgk = Reaction('PGK')
    pgk.add_metabolites({dpg13_c: -1, adp_c: -1, pg3_c: 1, atp_c: 1})
    pgk.lower_bound = -1000

    # PGM: 3PG <-> 2PG
    pgm = Reaction('PGM')
    pgm.add_metabolites({pg3_c: -1, pg2_c: 1})
    pgm.lower_bound = -1000

    # ENO: 2PG <-> PEP + H2O
    eno = Reaction('ENO')
    eno.add_metabolites({pg2_c: -1, pep_c: 1, h2o_c: 1})
    eno.lower_bound = -1000

    # PYK: PEP + ADP + H -> Pyr + ATP (注意这里通常消耗一个质子)
    pyk = Reaction('PYK')
    pyk.add_metabolites({pep_c: -1, adp_c: -1, h_c: -1, pyr_c: 1, atp_c: 1})

    # LDH: Pyr + NADH + H <-> Lac + NAD
    ldh = Reaction('LDH_L')
    ldh.add_metabolites({pyr_c: -1, nadh_c: -1, h_c: -1, lac_c: 1, nad_c: 1})
    ldh.lower_bound = -1000

    # ATPM: ATP + H2O -> ADP + Pi + H (模拟ATP消耗/维持)
    atpm = Reaction('ATPM')
    atpm.add_metabolites({atp_c: -1, h2o_c: -1, adp_c: 1, pi_c: 1, h_c: 1})
    atpm.lower_bound = 0
    atpm.upper_bound = 1000

    # --- 3. 定义交换反应 (System Boundaries) ---
    # 关键修复：必须允许 H+ 和 H2O 进出，否则质量无法平衡

    # EX_glc: 葡萄糖进入
    ex_glc = Reaction('EX_glc_e')
    ex_glc.add_metabolites({glc_c: -1})
    ex_glc.lower_bound = -10.0
    ex_glc.upper_bound = 0

    # EX_lac: 乳酸排出
    ex_lac = Reaction('EX_lac_e')
    ex_lac.add_metabolites({lac_c: -1})
    ex_lac.lower_bound = 0
    ex_lac.upper_bound = 1000

    # EX_h: 质子交换 (关键修复！！！)
    ex_h = Reaction('EX_h_e')
    ex_h.add_metabolites({h_c: -1})
    ex_h.lower_bound = -1000  # 允许进
    ex_h.upper_bound = 1000  # 允许出

    # EX_h2o: 水交换 (关键修复！！！)
    ex_h2o = Reaction('EX_h2o_e')
    ex_h2o.add_metabolites({h2o_c: -1})
    ex_h2o.lower_bound = -1000
    ex_h2o.upper_bound = 1000

    # EX_pi: 磷酸盐交换 (为了保险起见，防止Pi死锁，虽然理论上ATP循环可闭合)
    ex_pi = Reaction('EX_pi_e')
    ex_pi.add_metabolites({pi_c: -1})
    ex_pi.lower_bound = -1000
    ex_pi.upper_bound = 1000

    # 添加所有反应
    model.add_reactions([hex1, pgi, pfk, fba, tpi, gapd, pgk, pgm, eno, pyk, ldh, atpm,
                         ex_glc, ex_lac, ex_h, ex_h2o, ex_pi])

    return model


if __name__ == "__main__":
    model = build_corrected_model()

    # 设定目标函数
    model.objective = "EX_lac_e"

    solution = model.optimize()

    print(f"模型求解状态: {solution.status}")
    print(f"葡萄糖摄入速率: {-model.reactions.EX_glc_e.flux:.2f}")
    print(f"乳酸合成速率 (预测): {solution.objective_value:.2f}")
    print("-" * 30)
    print("关键内部通量详情:")
    print(f"  HEX1 (Input):   {solution.fluxes['HEX1']:.2f}")
    print(f"  GAPD (NADH Gen):{solution.fluxes['GAPD']:.2f}")
    print(f"  LDH  (Lac Gen): {solution.fluxes['LDH_L']:.2f}")
    print(f"  ATPM (ATP Use): {solution.fluxes['ATPM']:.2f}")
    print(f"  EX_h_e (H+ Out):{solution.fluxes['EX_h_e']:.2f}")


