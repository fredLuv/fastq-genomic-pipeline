import os
import csv
import json
import time
import urllib.request
import io
from datetime import datetime

# Setup directories relative to the script location
SCRATCH_DIR = os.path.dirname(os.path.abspath(__file__))
WORKDIR = os.path.dirname(SCRATCH_DIR)
RESULTS_DIR = os.path.join(WORKDIR, "results")
READS_DIR = os.path.join(WORKDIR, "reads")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(READS_DIR, exist_ok=True)

# 27 targeted SNPs from generic_gene_feature_matrix.py
snp_database = {
    # Methylation & Brain Health
    "rs1801133": {"Gene": "MTHFR", "Category": "Methylation", "Trait": "Folate metabolism (C677T)", "Risk_Allele": "T", "Notes": "TT: ~16% folate conversion; CT: ~65%; CC: Typical"},
    "rs1801131": {"Gene": "MTHFR", "Category": "Methylation", "Trait": "Folate metabolism (A1298C)", "Risk_Allele": "C", "Notes": "CC: Reduced folate conversion; AC: Mildly reduced"},
    "rs4680": {"Gene": "COMT", "Category": "Methylation/Brain", "Trait": "Dopamine clearance (Val158Met)", "Risk_Allele": "A", "Notes": "AA (Met/Met): Slow dopamine breakdown ('Worrier'); GG (Val/Val): Fast breakdown ('Warrior'); AG: Intermediate"},
    "rs6265": {"Gene": "BDNF", "Category": "Brain Health", "Trait": "Neuroplasticity (Val66Met)", "Risk_Allele": "A", "Notes": "AA/AG: Met allele associated with lower BDNF secretion and memory/stress resilience variation"},
    "rs1800497": {"Gene": "DRD2", "Category": "Brain Health", "Trait": "Dopamine D2 Receptor (Taq1A)", "Risk_Allele": "A", "Notes": "AA/AG: Reduced receptor density; associated with reward-seeking behavior and addiction susceptibility"},
    "rs6318": {"Gene": "MAOA", "Category": "Brain Health", "Trait": "Monoamine oxidase A activity", "Risk_Allele": "G", "Notes": "G allele associated with altered serotonin/dopamine degradation rates"},

    # Cardiovascular, Lipids & Longevity
    "rs429358": {"Gene": "APOE", "Category": "Cardiovascular/Lipids", "Trait": "Alzheimer's & lipid risk (C112R)", "Risk_Allele": "C", "Notes": "C allele determines APOE-epsilon4 (Increased Alzheimer's & lipid risk)"},
    "rs7412": {"Gene": "APOE", "Category": "Cardiovascular/Lipids", "Trait": "Alzheimer's & lipid risk (R158C)", "Risk_Allele": "T", "Notes": "T allele determines APOE-epsilon2 (Reduced Alzheimer's risk)"},
    "rs1799983": {"Gene": "NOS3", "Category": "Cardiovascular", "Trait": "Nitric Oxide Synthase (vasodilation)", "Risk_Allele": "T", "Notes": "TT/GT: Reduced nitric oxide production, associated with increased hypertension risk"},
    "rs2802292": {"Gene": "FOXO3", "Category": "Longevity", "Trait": "Healthy aging & longevity", "Risk_Allele": "C", "Notes": "CC/CG: Associated with increased lifespan, protection against cardiovascular disease"},

    # Diet, Lactose & Metabolism
    "rs4988235": {"Gene": "LCT", "Category": "Dietary Tolerance", "Trait": "Lactase persistence (Lactose tolerance)", "Risk_Allele": "G", "Notes": "GG: Lactose intolerant; AA/AG: Lactose tolerant"},
    "rs9939609": {"Gene": "FTO", "Category": "Metabolism", "Trait": "Weight and obesity risk", "Risk_Allele": "A", "Notes": "AA: Increased obesity risk and fat mass; TA: Intermediate; TT: Typical risk"},
    "rs1801282": {"Gene": "PPARG", "Category": "Metabolism", "Trait": "Insulin sensitivity & fat storage", "Risk_Allele": "G", "Notes": "GG/CG: Typical sensitivity; CC (Pro/Pro): Associated with higher BMI but lower diabetes risk"},
    "rs1801260": {"Gene": "CLOCK", "Category": "Circadian Rhythm", "Trait": "Sleep duration & night owl tendency", "Risk_Allele": "G", "Notes": "GG/AG: Associated with evening chronotype ('night owl') and shorter sleep duration"},

    # Vitamins & Iron
    "rs1800562": {"Gene": "HFE", "Category": "Iron Metabolism", "Trait": "Hemochromatosis (C282Y)", "Risk_Allele": "A", "Notes": "AA: High risk of iron overload; GA: Carrier; GG: Typical"},
    "rs1799945": {"Gene": "HFE", "Category": "Iron Metabolism", "Trait": "Hemochromatosis (H63D)", "Risk_Allele": "G", "Notes": "GG: Moderate risk of iron overload; CG: Carrier; CC: Typical"},
    "rs731236": {"Gene": "VDR", "Category": "Vitamin Absorption", "Trait": "Vitamin D Receptor (TaqI)", "Risk_Allele": "G", "Notes": "GG (tt): Associated with reduced vitamin D absorption/bone density; AA: Normal"},
    "rs7501331": {"Gene": "BCMO1", "Category": "Vitamin Absorption", "Trait": "Beta-carotene to Vitamin A conversion", "Risk_Allele": "C", "Notes": "CC/TC: Up to 57% reduction in conversion of plant-based beta-carotene to active Vitamin A"},
    "rs602662": {"Gene": "FUT2", "Category": "Vitamin Absorption", "Trait": "Vitamin B12 absorption", "Risk_Allele": "A", "Notes": "AA: Associated with lower B12 absorption and plasma levels; GG: Normal"},
    "rs4654748": {"Gene": "NBPF3", "Category": "Vitamin Absorption", "Trait": "Vitamin B6 clearance rate", "Risk_Allele": "T", "Notes": "TT/CT: Associated with faster B6 clearance and lower plasma levels"},

    # Detoxification & Inflammation
    "rs1695": {"Gene": "GSTP1", "Category": "Detoxification", "Trait": "Glutathione S-transferase (detox speed)", "Risk_Allele": "G", "Notes": "GG: Reduced detoxification capacity for toxins/heavy metals; AG: Intermediate; AA: Normal"},
    "rs4880": {"Gene": "SOD2", "Category": "Antioxidant", "Trait": "Superoxide dismutase (antioxidant defense)", "Risk_Allele": "T", "Notes": "TT: Reduced transport of SOD2 to mitochondria, associated with higher oxidative stress; CC: Normal"},
    "rs671": {"Gene": "ALDH2", "Category": "Detoxification", "Trait": "Alcohol flush reaction (Glu504Lys)", "Risk_Allele": "A", "Notes": "AA: Severe flush/hangover; GA: Mild flush; GG: Normal tolerance"},
    "rs1800795": {"Gene": "IL6", "Category": "Inflammation", "Trait": "Interleukin-6 systemic inflammation", "Risk_Allele": "G", "Notes": "GG/CG: Associated with higher baseline IL-6 levels and systemic inflammation"},
    "rs1800629": {"Gene": "TNF", "Category": "Inflammation", "Trait": "Tumor Necrosis Factor (inflammatory response)", "Risk_Allele": "A", "Notes": "AA/GA: Higher baseline TNF-alpha, associated with increased inflammatory/autoimmune risk"},

    # Traits & Physical Characteristics
    "rs1815739": {"Gene": "ACTN3", "Category": "Physical Traits", "Trait": "Muscle performance (R577X)", "Risk_Allele": "T", "Notes": "CC: Sprinters/Power athletes; CT: Intermediate; TT: Endurance athletes"},
    "rs1805007": {"Gene": "MC1R", "Category": "Physical Traits", "Trait": "Red hair & fair skin", "Risk_Allele": "T", "Notes": "TT: High probability of red hair; CT: Carrier; CC: Typical"}
}

# Fallback cache to guarantee offline functionality and rapid rendering
FALLBACK_METADATA = {
  "rs1801133": {"Chromosome": "1", "Position": 11796321, "Consequence": "missense_variant", "Official_Gene": "MTHFR"},
  "rs1801131": {"Chromosome": "1", "Position": 11794419, "Consequence": "missense_variant", "Official_Gene": "C1orf167"},
  "rs4680": {"Chromosome": "22", "Position": 19963748, "Consequence": "missense_variant", "Official_Gene": "COMT"},
  "rs6265": {"Chromosome": "11", "Position": 27658369, "Consequence": "missense_variant", "Official_Gene": "BDNF"},
  "rs1800497": {"Chromosome": "11", "Position": 113400106, "Consequence": "missense_variant", "Official_Gene": "ANKK1"},
  "rs6318": {"Chromosome": "X", "Position": 114731326, "Consequence": "missense_variant", "Official_Gene": "HTR2C"},
  "rs429358": {"Chromosome": "19", "Position": 44908684, "Consequence": "missense_variant", "Official_Gene": "APOE"},
  "rs7412": {"Chromosome": "19", "Position": 44908822, "Consequence": "missense_variant", "Official_Gene": "APOE"},
  "rs1799983": {"Chromosome": "7", "Position": 150999023, "Consequence": "missense_variant", "Official_Gene": "NOS3"},
  "rs2802292": {"Chromosome": "6", "Position": 108587315, "Consequence": "intron_variant", "Official_Gene": "FOXO3"},
  "rs4988235": {"Chromosome": "2", "Position": 135851076, "Consequence": "intron_variant", "Official_Gene": "MCM6"},
  "rs9939609": {"Chromosome": "16", "Position": 53786615, "Consequence": "intron_variant", "Official_Gene": "FTO"},
  "rs1801282": {"Chromosome": "3", "Position": 12351626, "Consequence": "missense_variant", "Official_Gene": "PPARG"},
  "rs1801260": {"Chromosome": "4", "Position": 55435202, "Consequence": "3_prime_UTR_variant", "Official_Gene": "CLOCK"},
  "rs1800562": {"Chromosome": "6", "Position": 26092913, "Consequence": "missense_variant", "Official_Gene": "HFE"},
  "rs1799945": {"Chromosome": "6", "Position": 26090951, "Consequence": "missense_variant", "Official_Gene": "HFE"},
  "rs731236": {"Chromosome": "12", "Position": 47844974, "Consequence": "synonymous_variant", "Official_Gene": "VDR"},
  "rs7501331": {"Chromosome": "16", "Position": 81280891, "Consequence": "missense_variant", "Official_Gene": "BCO1"},
  "rs602662": {"Chromosome": "19", "Position": 48703728, "Consequence": "missense_variant", "Official_Gene": "FUT2"},
  "rs4654748": {"Chromosome": "1", "Position": 21459575, "Consequence": "intron_variant", "Official_Gene": "NBPF3"},
  "rs1695": {"Chromosome": "11", "Position": 67585218, "Consequence": "missense_variant", "Official_Gene": "GSTP1"},
  "rs4880": {"Chromosome": "6", "Position": 159692840, "Consequence": "missense_variant", "Official_Gene": "SOD2"},
  "rs671": {"Chromosome": "12", "Position": 111803962, "Consequence": "missense_variant", "Official_Gene": "ALDH2"},
  "rs1800795": {"Chromosome": "7", "Position": 22727026, "Consequence": "5_prime_UTR_variant", "Official_Gene": "IL6"},
  "rs1800629": {"Chromosome": "HSCHR6_MHC_SSTO_CTG1", "Position": 2874534, "Consequence": "upstream_gene_variant", "Official_Gene": "LTA"},
  "rs1815739": {"Chromosome": "11", "Position": 66560624, "Consequence": "stop_gained", "Official_Gene": "CTSF"},
  "rs1805007": {"Chromosome": "16", "Position": 89919709, "Consequence": "missense_variant", "Official_Gene": "TUBB3"}
}

def get_patient_interpretation(rsid, risk_count, genotype):
    """Applying Chinese-English bilingual interpretations from generate_html_report.py."""
    genotype = (genotype or "--").upper().strip()
    
    if rsid == "rs1801133":
        if risk_count == 2:
            return "<b>MTHFR C677T 突变 (TT)</b>：MTHFR enzyme activity 较常人降低约 84%，导致 homocysteine (同型半胱氨酸) 升高风险显著增加。建议补充 active methylfolate (5-MTHF)。"
        elif risk_count == 1:
            return "<b>MTHFR C677T 杂合 (CT)</b>：MTHFR enzyme activity 降低约 35%。日常建议增加 leafy green vegetables，适当补充 active folate。"
        else:
            return "<b>MTHFR C677T 野生型 (CC)</b>：未发现 C677T mutation，具有典型的 folate conversion 和 methylation 代谢速度。"
            
    elif rsid == "rs1801131":
        if risk_count == 2:
            return "<b>MTHFR A1298C 突变 (CC)</b>：MTHFR enzyme activity 降低约 40%，伴随中度的 methyl donor 不足倾向。"
        elif risk_count == 1:
            return "<b>MTHFR A1298C 杂合 (GT/AC)</b>：典型的 MTHFR enzyme activity，体内 methylation 处于平衡水平。"
        else:
            return "<b>MTHFR A1298C 野生型 (AA)</b>：未检出 A1298C mutation，具有典型的 folate utilization 能力。"
            
    elif rsid == "rs4680":
        if risk_count == 0:
            return "<b>Warrior Genotype (GG)</b>：dopamine breakdown (多巴胺降解) 速率极快。在 high-pressure environments 中冷静沉着、focus 极高，但在 low-stress 日常环境中容易觉得 bored、缺乏 motivation。"
        elif risk_count == 1:
            return "<b>Balanced Genotype (AG)</b>：dopamine breakdown 速率中等。在日常 focus 与应对 sudden stress 之间具有健康的 balance。"
        else:
            return "<b>Worrier Genotype (AA)</b>：dopamine breakdown 极其缓慢。在 calm environment 中拥有极强的 cognitive depth、focus 与 memory，但在 high-pressure 压力下极易感到 anxiety 与 stress sensitivity。"
            
    elif rsid == "rs6265":
        if risk_count >= 1:
            return "<b>BDNF Val66Met 突变 (AG/AA)</b>：Met allele 会导致 activity-dependent BDNF secretion 轻度下调，略微影响 neuroplasticity (神经可塑性)。建议通过 regular cardiovascular exercise 有效刺激 BDNF 分泌。"
        else:
            return "<b>Typical BDNF Secretion (CC)</b>：典型的 BDNF 分泌活性与 synaptic plasticity，认知保护与 brain repair 机制正常。"
            
    elif rsid == "rs1800497":
        if risk_count == 1:
            return "<b>Dopamine D2 Receptor 密度降低 (AG)</b>：Taq1A carrier。大脑 dopamine D2 receptor 数量较典型值减少约 30%，容易表现出更高的 reward threshold (较易沉迷甜食、刺激或 novelty-seeking 行为)。"
        elif risk_count == 2:
            return "<b>Dopamine D2 Receptor 密度显著降低 (AA)</b>：dopamine receptor 数量减少超 30%。天然 reward 满足感较低，容易陷入冲动行为或寻求外界强刺激。"
        else:
            return "<b>Dopamine Receptor 正常 (GG)</b>：具有典型的 dopamine D2 receptor 密度与 reward pathways 响应机制。"
            
    elif rsid == "rs6318":
        return "<b>Typical MAO-A Activity</b>：典型的 monoamine oxidase A 活性，血清素和多巴胺清除速度正常，情绪耐受力处于 baseline 范围。"
        
    elif rsid == "rs429358":
        if "C" in genotype:
            return "<b>APOE-e4 Risk Allele 携带者</b>：Alzheimer's disease 及 LDL-cholesterol (低密度脂蛋白胆固醇) 升高的遗传风险增加。日常建议采用 Mediterranean diet 并定期复查 lipid panel。"
        else:
            return "<b>APOE-e4 野生型背景 (TT)</b>：未检测到 APOE-e4 风险基因，该位点表现为健康的 typical 单倍型背景。"
            
    elif rsid == "rs7412":
        if "T" in genotype:
            return "<b>APOE-e2 Protective Allele 携带者</b>：与降低 Alzheimer's disease 发病风险以及改善 lipid clearance 脂质清除相关。"
        else:
            return "<b>APOE-e3 正常背景 (CC)</b>：未发现该位点突变，脂质代谢基准正常。"
            
    elif rsid == "rs1799983":
        if risk_count >= 1:
            return "<b>Nitric Oxide 生成受限 (TT/GT)</b>：eNOS (血管内皮型一氧化氮合成酶) 活性略微降低，可能影响 vasodilation (血管舒张) 和微循环。日常可多补充 beetroot 等 nitric oxide 食物来源。"
        else:
            return "<b>Nitric Oxide 生成正常 (GG)</b>：血管收缩与舒张平衡正常，具有健康的 microcirculation (微循环) 机能。"
            
    elif rsid == "rs2802292":
        if "G" in genotype or "C" in genotype:
            return "<b>FOXO3 Longevity Variant 携带者 (GT)</b>：能有效激活 cell autophagy (细胞自噬) 与 antioxidant enzymes 的转录，增强 cardiovascular 血管内皮的健康老化保护。"
        else:
            return "<b>Typical Aging Profile (TT)</b>：细胞抗氧化酶系统基底功能正常，老化速度处于平均状态。"
            
    elif rsid == "rs4988235":
        if risk_count == 2:
            return "<b>Lactose Intolerance (GG)</b>：小肠粘膜 lactase enzyme 发生 age-related 完全停产。饮用普通 milk 或 dairy 后极易发生 bloating (腹胀)、gas 或 IBS (肠易激)。"
        else:
            return "<b>Lactose Tolerant (AA/AG)</b>：lactase enzyme 在成年期持续产生，能正常水解和消化乳糖。"
            
    elif rsid == "rs9939609":
        if risk_count == 2:
            return "<b>FTO Obesity & Satiety 风险 (AA)</b>：ghrelin (饥饿素) 抑制迟缓，satiety (饱腹感) 信号延迟，极易无意识摄入额外卡路里。建议执行 high-satiety、high-protein 膳食策略。"
        elif risk_count == 1:
            return "<b>FTO 中度 Obesity 风险 (TA)</b>：食欲控制及脂肪囤积倾向中等，建议适当控制 refined carbs 摄入。"
        else:
            return "<b>Typical FTO Metabolic Profile (TT)</b>：正常的 satiety control 反应，未表现出明显的易胖体质遗传倾向。"
            
    elif rsid == "rs1801282":
        if risk_count >= 1:
            return "<b>PPARG Insulin Sensitivity 优化 (CG)</b>：虽然与略微升高的 BMI 敏感度相关，但在防范 type 2 diabetes 及保护血管内皮 insulin sensitivity 方面具备遗传优势。"
        else:
            return "<b>Typical PPARG Activity (CC)</b>：处于人群均值水平的脂肪分化与 insulin 结合度。"
            
    elif rsid == "rs1801260":
        if risk_count >= 1:
            return "<b>Night Owl Chronotype (AG/GG)</b>：CLOCK 基因轻微位移，表现为典型的夜间型 chronotype，入睡时间偏晚且清晨精力恢复较慢。"
        else:
            return "<b>Morning Chronotype (AA)</b>：自然健康的 circadian rhythm (昼夜节律)，晨间精力充沛。"
            
    elif rsid == "rs1800562":
        if risk_count == 2:
            return "<b>Hemochromatosis 铁过载高危 (AA)</b>：肠道铁吸收失控风险高。需定期监测 Ferritin (铁蛋白) 与 transferrin saturation (转铁蛋白饱和度)。"
        elif risk_count == 1:
            return "<b>Hemochromatosis Carrier (GA)</b>：通常不会引发任何铁沉积症状，体内铁负荷完全正常。"
        else:
            return "<b>Typical Iron Absorption (GG)</b>：体内铁吸收与代谢负荷受健康的 Hepcidin 系统严格调控，无积铁危险。"
            
    elif rsid == "rs1799945":
        if risk_count >= 1:
            return "<b>Hemochromatosis H63D Carrier (CG)</b>：肠道铁吸收无显著异常。除非与 C282Y 双重杂合，否则 ferritin 水平正常。"
        else:
            return "<b>Typical Iron Metabolism (CC)</b>：正常的血色素与铁离子吸收机制。"
            
    elif rsid == "rs731236":
        if risk_count >= 1:
            return "<b>Vitamin D Receptor 效率降低 (tt/Gt)</b>：肠道对 calcium 和 active Vitamin D3 的吸收与结合效率略低。建议定期监测 Vit D3 浓度，必要时补充 Vit D3 + K2。"
        else:
            return "<b>Typical Vitamin D Receptor (AA)</b>：VDR 结合性能极佳，拥有健康的 bone mineralization (骨矿化) 保护力。"
            
    elif rsid == "rs7501331":
        if risk_count == 2:
            return "<b>Vitamin A Conversion 降低 57% (CC)</b>：BCMO1 酶活性降低达 57%。从 plant-based beta-carotene 转化为 active Vitamin A (retinol) 效率极差。日常必需直接摄入 retinol (蛋黄、动物肝脏)。"
        elif risk_count == 1:
            return "<b>Vitamin A Conversion 降低 32% (TC)</b>：$\beta$-carotene 转化 active Vit A 效率轻度降低。"
        else:
            return "<b>Typical Vitamin A Conversion (TT)</b>：能将食物中的 beta-carotene 高效转化为 active Vitamin A。"
            
    elif rsid == "rs602662":
        if risk_count == 2:
            return "<b>Vitamin B12 Absorption 降低 (AA)</b>：胃肠道对 B12 的特异性结合吸收偏慢，易发生 B12 不足。建议补充 active methylcobalamin。"
        elif risk_count == 1:
            return "<b>Vitamin B12 Absorption 轻度受限 (AG)</b>：体内 active B12 蓄积率处于轻微偏低范围。"
        else:
            return "<b>Typical B12 Absorption (GG)</b>：具有高水准的 B12 生物利用度与胃粘膜吸收速度。"
            
    elif rsid == "rs4654748":
        if risk_count >= 1:
            return "<b>Fast Vitamin B6 Clearance</b>：由于转运体表达差异，体内的 B6 在血液中清除降解偏快。建议适当高频获取富含 B6 的膳食。"
        else:
            return "<b>Typical Vitamin B6 Levels</b>：正常且稳定的 B6 留存率。"
            
    elif rsid == "rs1695":
        if risk_count == 2:
            return "<b>Glutathione Detoxification 降低 (GG)</b>：GSTP1 纯合突变导致细胞清除重金属、环境毒素的速度降低。建议多吃 cruciferous vegetables (十字花科蔬菜) 提升 GST 活性。"
        elif risk_count == 1:
            return "<b>Glutathione Detoxification 轻度降低 (AG)</b>：细胞抗毒物与解毒酶活性中等。"
        else:
            return "<b>Typical Glutathione Detoxification (AA)</b>：高效且典型的细胞第一/二阶段解毒排毒循环系统。"
            
    elif rsid == "rs4880":
        if risk_count == 2:
            return "<b>Mitochondrial Antioxidant 偏低 (SOD2 TT)</b>：锰型超氧化物歧化酶向线粒体转移受阻。mitochondrial oxidative stress 升高，需增加饮食中强抗氧化剂(浆果、绿茶)。"
        else:
            return "<b>Typical Mitochondrial Defense (GG/GT)</b>：充足的线粒体内部自由基清除能力，能量工厂保护力强。"
            
    elif rsid == "rs671":
        if risk_count == 2:
            return "<b>Severe Alcohol Flush Reaction (AA)</b>：ALDH2 酶活性接近为 0。乙醛 (acetaldehyde) 在体内迅速蓄积并引起 facial flushing、头痛、心悸等，严禁饮酒。"
        elif risk_count == 1:
            return "<b>Moderate Alcohol Flush Reaction (GA)</b>：ALDH2 活性降低约 60-80%，酒后极易脸红，宿醉感强烈。建议限制酒精摄入。"
        else:
            return "<b>Typical Alcohol Tolerance (GG)</b>：体内 ALDH2 酶活性充沛，酒精与乙醛转化排泄速度正常，无红脸反应。"
            
    elif rsid == "rs1800795":
        if risk_count == 2:
            return "<b>Elevated Baseline Systemic Inflammation (GG)</b>：促炎因子 IL-6 基因处于活跃表达状态，机体易进入 low-grade chronic inflammation (慢性低度炎症) 状态。建议遵循 anti-inflammatory diet。"
        else:
            return "<b>Typical Inflammatory Baseline</b>：典型的促炎及抑炎反应通路，机体炎性基准处于平均水平。"
            
    elif rsid == "rs1800629":
        if risk_count >= 1:
            return "<b>Elevated Inflammatory Response (AG/AA)</b>：TNF-alpha 反应偏高。受到细菌或毒素刺激后免疫炎性反应启动更快。推荐补充 Omega-3 脂肪酸。"
        else:
            return "<b>Typical TNF-alpha Response (GG)</b>：TNF-alpha 的基线生成完全正常，自身免疫反应平稳。"
            
    elif rsid == "rs1815739":
        if risk_count == 2:
            return "<b>Endurance Adapted (TT)</b>：快肌纤维中 complete absence of actinin-3 protein。天然缺乏 explosive power 项目优势，但在 endurance sports (马拉松、长跑) 中表现出卓越的 oxygen efficiency。"
        elif risk_count == 1:
            return "<b>Balanced Muscle Composition (CT)</b>：肌纤维兼具爆发力蛋白质与有氧代谢的折中配置，适合各类复合型运动。"
        else:
            return "<b>Sprint & Power Adapted (CC)</b>：快肌纤维富含 actinin-3 protein，收缩速度优异，适合 explosive power 项目（冲刺、举重）。"
            
    elif rsid == "rs1805007":
        if risk_count == 2:
            return "<b>Red Hair / Fair Skin (TT)</b>：MC1R 黑色素受体功能完全失活。皮肤极易 sunburnt (晒伤) 而无法 tanned (晒黑)，有较高的 UV sensitivity。"
        elif risk_count == 1:
            return "<b>Red Hair Carrier (CT)</b>：正常的黑色素分泌及典型防晒伤皮层状态。"
        else:
            return "<b>Typical Pigmentation (CC)</b>：真黑色素分泌正常，皮肤能健康的适应紫外线并生成防晒古铜色。"
            
    return "检出正常野生型或中性基因型。"

def get_snp_metadata(rsid):
    """Dynamic resolution of missing SNPs metadata using local cache or API backup."""
    cache_path = os.path.join(RESULTS_DIR, "ensembl_metadata_cache.json")
    cache = {}
    if os.path.exists(cache_path):
        try:
            with open(cache_path, 'r') as f:
                cache = json.load(f)
        except Exception:
            pass
            
    if rsid in cache:
        return cache[rsid]
        
    if rsid in FALLBACK_METADATA:
        cache[rsid] = FALLBACK_METADATA[rsid]
        try:
            with open(cache_path, 'w') as f:
                json.dump(cache, f, indent=2)
        except Exception:
            pass
        return FALLBACK_METADATA[rsid]
        
    # Ensembl lookup fallback
    url = f"https://rest.ensembl.org/vep/human/id/{rsid}?content-type=application/json"
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req) as r:
            res = json.loads(r.read().decode('utf-8'))
            if res and len(res) > 0:
                consequence = res[0].get("most_severe_consequence", "unknown")
                chrom = res[0].get("seq_region_name", "unknown")
                start = res[0].get("start", "unknown")
                
                gene_symbol = "unknown"
                tc = res[0].get("transcript_consequences", [])
                for t in tc:
                    if "gene_symbol" in t:
                        gene_symbol = t["gene_symbol"]
                        break
                        
                meta = {
                    "Chromosome": chrom,
                    "Position": start,
                    "Consequence": consequence,
                    "Official_Gene": gene_symbol
                }
                cache[rsid] = meta
                with open(cache_path, 'w') as f:
                    json.dump(cache, f, indent=2)
                return meta
    except Exception as e:
        print(f"Ensembl API lookup failed for {rsid}: {e}")
        
    return {
        "Chromosome": "unknown",
        "Position": "unknown",
        "Consequence": "unknown",
        "Official_Gene": "unknown"
    }

def run_mapping_logic(file_content, filename):
    """Processes raw genotyping file in memory and compiles genomic data payload."""
    genotype_map = {}
    total_parsed_lines = 0
    
    # Parse lines
    for line in file_content.splitlines():
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            genotype_map[parts[0]] = parts[3]
            total_parsed_lines += 1
            
    print(f"Parsed {len(genotype_map)} genotypes from {filename}")
    
    # Process SNPs
    snp_data = {}
    for rsid, target in snp_database.items():
        meta = get_snp_metadata(rsid)
        genotype = genotype_map.get(rsid, "--")
        
        # Calculate Risk Allele Counts
        risk_allele = target["Risk_Allele"]
        if genotype == "--" or genotype == "null" or not genotype:
            risk_count = 0
            status = "Not Genotyped"
        else:
            risk_count = sum(1 for base in genotype if base == risk_allele)
            if risk_count == 2:
                status = "Homozygous Risk"
            elif risk_count == 1:
                status = "Heterozygous Risk"
            else:
                status = "Homozygous Wildtype"
                
        patient_implication = get_patient_interpretation(rsid, risk_count, genotype)
        
        gene_symbol = target["Gene"]
        if gene_symbol not in snp_data:
            snp_data[gene_symbol] = []
            
        snp_data[gene_symbol].append({
            "rsID": rsid,
            "Chromosome": meta["Chromosome"],
            "Position": meta["Position"],
            "Consequence": meta["Consequence"],
            "Genotype": genotype,
            "Risk_Allele": risk_allele,
            "Risk_Count": risk_count,
            "Status": status,
            "Trait": target["Trait"],
            "Notes": target["Notes"],
            "interpretation": patient_implication
        })
        
    # Aggregate by Gene
    genes = []
    for gene_name, snps in snp_data.items():
        category = snps[0]["Category"] if snps else "General"
        chromosome = snps[0]["Chromosome"] if snps else "unknown"
        
        total_snps = len(snps)
        genotyped_snps = sum(1 for s in snps if s["Genotype"] != "--")
        total_risk_alleles = sum(s["Risk_Count"] for s in snps)
        
        total_alleles = genotyped_snps * 2
        avg_risk_score = round(total_risk_alleles / total_alleles, 3) if total_alleles > 0 else 0.0
        
        # Consequences set
        consequences = set(s["Consequence"] for s in snps if s["Consequence"] != "unknown")
        consequence_string = "; ".join(consequences) if consequences else "unknown"
        
        if len(snps) == 1:
            implication = snps[0]["interpretation"]
            trait = snps[0]["Trait"]
        elif len(snps) > 1:
            implications = []
            for s in snps:
                if s["Genotype"] != "--":
                    clean_interp = s["interpretation"].replace("<b>", "").replace("</b>", "")
                    implications.append(f"【{s['rsID']}】{clean_interp}")
            implication = "; <br>".join(implications) if implications else "未检出有效测序位点。"
            trait = snps[0]["Trait"].split("(")[0].strip() + " (多位点组合)"
        else:
            implication = "无可用突变解释。"
            trait = "未知突变性状"
            
        genes.append({
            "name": gene_name,
            "category": category,
            "chromosome": chromosome,
            "total_snps": total_snps,
            "genotyped_snps": genotyped_snps,
            "risk_alleles": total_risk_alleles,
            "risk_score": avg_risk_score,
            "consequences": consequence_string,
            "snps": snps,
            "trait": trait,
            "implication": implication
        })
        
    return genes, total_parsed_lines

# Flask server declaration
try:
    from flask import Flask, request, render_template, make_response, redirect, url_for
    FLASK_AVAILABLE = True
except ImportError:
    FLASK_AVAILABLE = False

if FLASK_AVAILABLE:
    app = Flask(__name__, template_folder='templates')
    
    @app.route('/')
    def index():
        return render_template('upload.html')
        
    @app.route('/upload', methods=['POST'])
    def upload_file():
        use_demo = request.form.get('use_demo')
        
        if use_demo == 'true':
            # Load reads/23andme.txt
            demo_path = os.path.join(READS_DIR, "23andme.txt")
            if not os.path.exists(demo_path):
                # Try downloading or fallback to generating toy data
                # Since 23andme.txt exists under reads, let's load it
                demo_path = "/Users/fred/Library/CloudStorage/GoogleDrive-fredlu0521@gmail.com/My Drive/Ethan-Resume/fastq_pipeline/reads/23andme.txt"
                
            if os.path.exists(demo_path):
                with open(demo_path, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                filename = "Lorena_Sandoval_v5_Full.txt"
            else:
                return "Sample 23andMe file not found on system. Please upload a file.", 404
        else:
            # Handle uploaded file
            if 'file' not in request.files:
                return redirect(url_for('index'))
            file = request.files['file']
            if file.filename == '':
                return redirect(url_for('index'))
            
            content = file.read().decode('utf-8', errors='ignore')
            filename = file.filename
            
        # Parse and process
        genes, total_genotypes = run_mapping_logic(content, filename)
        
        # Render the template and inject data
        genes_json = json.dumps(genes, indent=4)
        upload_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        return render_template(
            'report.html',
            filename=filename,
            upload_time=upload_time,
            genotype_count=total_genotypes,
            genes_json=genes_json
        )
        
    @app.route('/download', methods=['POST'])
    def download_csv():
        download_type = request.form.get('type')
        payload_json = request.form.get('payload_json')
        
        if not payload_json:
            return "Missing payload", 400
            
        try:
            genes_data = json.loads(payload_json)
        except Exception as e:
            return f"Invalid JSON payload: {e}", 400
            
        output = io.StringIO()
        
        if download_type == 'snps':
            fieldnames = [
                "rsID", "Gene_Symbol", "Category", "Chromosome", "Position_GRCh38", "Consequence",
                "Genotype", "Risk_Allele", "Risk_Allele_Count", "Genotype_Status", "Trait_Association", "Clinical_Implications"
            ]
            writer = csv.DictWriter(output, fieldnames=fieldnames)
            writer.writeheader()
            
            seen_snps = set()
            for gene in genes_data:
                for s in gene.get('snps', []):
                    rsid = s.get('rsID')
                    if rsid and rsid not in seen_snps:
                        seen_snps.add(rsid)
                        writer.writerow({
                            "rsID": rsid,
                            "Gene_Symbol": gene.get('name'),
                            "Category": gene.get('category'),
                            "Chromosome": s.get('Chromosome'),
                            "Position_GRCh38": s.get('Position'),
                            "Consequence": s.get('Consequence'),
                            "Genotype": s.get('Genotype'),
                            "Risk_Allele": s.get('Risk_Allele'),
                            "Risk_Allele_Count": s.get('Risk_Count'),
                            "Genotype_Status": s.get('Status'),
                            "Trait_Association": s.get('Trait'),
                            "Clinical_Implications": s.get('Notes')
                        })
            filename = "generic_snp_features.csv"
            
        elif download_type == 'genes':
            fieldnames = [
                "Gene_Symbol", "Category", "Chromosome", "Total_SNPs_Measured", 
                "Genotyped_SNPs", "Total_Risk_Alleles", "Average_Risk_Score", 
                "SNP_Genotypes", "Most_Severe_Consequences"
            ]
            writer = csv.DictWriter(output, fieldnames=fieldnames)
            writer.writeheader()
            
            for gene in genes_data:
                snp_genotypes = []
                for s in gene.get('snps', []):
                    if s.get('Genotype') != '--':
                        snp_genotypes.append(f"{s.get('rsID')}:{s.get('Genotype')}")
                snp_string = "; ".join(snp_genotypes) if snp_genotypes else "None"
                
                writer.writerow({
                    "Gene_Symbol": gene.get('name'),
                    "Category": gene.get('category'),
                    "Chromosome": gene.get('chromosome'),
                    "Total_SNPs_Measured": gene.get('total_snps'),
                    "Genotyped_SNPs": gene.get('genotyped_snps'),
                    "Total_Risk_Alleles": gene.get('risk_alleles'),
                    "Average_Risk_Score": gene.get('risk_score'),
                    "SNP_Genotypes": snp_string,
                    "Most_Severe_Consequences": gene.get('consequences')
                })
            filename = "generic_gene_features.csv"
        else:
            return "Invalid download type", 400
            
        resp = make_response(output.getvalue())
        resp.headers["Content-Disposition"] = f"attachment; filename={filename}"
        resp.headers["Content-type"] = "text/csv"
        return resp

if __name__ == '__main__':
    # Verify imports
    print("=========================================")
    print(" Checking Local Dependencies...")
    print("=========================================")
    import sys
    try:
        import flask
        print("[-] Flask imported successfully.")
    except ImportError:
        print("[!] Flask is NOT installed.")
        print("[!] Please run: pip install flask")
        sys.exit(1)
        
    print("\nStarting local Genomic SaaS Web Prototype server...")
    print("Please open: http://127.0.0.1:5050/")
    print("=========================================")
    app.run(host='127.0.0.1', port=5050, debug=True)
