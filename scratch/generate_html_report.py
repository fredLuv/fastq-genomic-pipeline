import os
import csv
import json

# Setup directories
WORKDIR = "/Users/fred/Library/CloudStorage/GoogleDrive-fredlu0521@gmail.com/My Drive/Ethan-Resume/fastq_pipeline"
RESULTS_DIR = os.path.join(WORKDIR, "results")

gene_csv_path = os.path.join(RESULTS_DIR, "generic_gene_features.csv")
snp_csv_path = os.path.join(RESULTS_DIR, "generic_snp_features.csv")
html_output_path = os.path.join(RESULTS_DIR, "genomic_report.html")

# Define Patient-Specific Clinical Interpretations (Bilingual / Chinglish Style)
def get_patient_interpretation(rsid, risk_count, genotype):
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
            return "<b>APOE-e4 野生型背景 (TT)</b>：未检测到 APOE-e4 风险基因，该位点表现为健康的典型单倍型背景。"
            
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

# Read SNP data
snp_data = {}
with open(snp_csv_path, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        gene = row["Gene_Symbol"]
        rsid = row["rsID"]
        risk_count = int(row["Risk_Allele_Count"])
        genotype = row["Genotype"]
        
        patient_implication = get_patient_interpretation(rsid, risk_count, genotype)
        
        if gene not in snp_data:
            snp_data[gene] = []
        snp_data[gene].append({
            "rsID": rsid,
            "Chromosome": row["Chromosome"],
            "Position": row["Position_GRCh38"],
            "Consequence": row["Consequence"],
            "Genotype": genotype,
            "Risk_Allele": row["Risk_Allele"],
            "Risk_Count": risk_count,
            "Status": row["Genotype_Status"],
            "Trait": row["Trait_Association"],
            "Notes": row["Clinical_Implications"],
            "interpretation": patient_implication
        })

# Read Gene data and merge
genes = []
with open(gene_csv_path, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        gene_name = row["Gene_Symbol"]
        snps = snp_data.get(gene_name, [])
        
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
            "category": row["Category"],
            "chromosome": row["Chromosome"],
            "total_snps": int(row["Total_SNPs_Measured"]),
            "genotyped_snps": int(row["Genotyped_SNPs"]),
            "risk_alleles": int(row["Total_Risk_Alleles"]),
            "risk_score": float(row["Average_Risk_Score"]),
            "consequences": row["Most_Severe_Consequences"],
            "snps": snps,
            "trait": trait,
            "implication": implication
        })

# HTML/CSS/JS template content
html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Anonymous Patient - Personal Genomic Health Report</title>
    <meta name="description" content="Interactive visual report of personal genetic traits, methylation, cardiovascular risk, and dietary tolerances.">
    <!-- Google Fonts -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=Outfit:wght@400;500;600;700;800&display=swap" rel="stylesheet">
    <!-- Chart.js -->
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        :root {
            --bg-primary: #0b0f19;
            --bg-secondary: #111827;
            --bg-card: rgba(30, 41, 59, 0.45);
            --bg-card-hover: rgba(30, 41, 59, 0.7);
            --border-color: rgba(255, 255, 255, 0.08);
            
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --text-muted: #64748b;
            
            --color-wildtype: #10b981;
            --color-carrier: #f59e0b;
            --color-variant: #f43f5e;
            
            --glow-wildtype: rgba(16, 185, 129, 0.15);
            --glow-carrier: rgba(245, 158, 11, 0.15);
            --glow-variant: rgba(244, 63, 94, 0.15);
        }

        * {
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }

        body {
            font-family: 'Inter', sans-serif;
            background-color: var(--bg-primary);
            color: var(--text-primary);
            line-height: 1.6;
            padding: 2rem 1.5rem;
            max-width: 1400px;
            margin: 0 auto;
        }

        h1, h2, h3, h4 {
            font-family: 'Outfit', sans-serif;
            font-weight: 700;
        }

        /* Top Header Card */
        header {
            background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
            border: 1px solid var(--border-color);
            border-radius: 1.25rem;
            padding: 2.5rem;
            margin-bottom: 2rem;
            box-shadow: 0 10px 30px rgba(0, 0, 0, 0.3);
            position: relative;
            overflow: hidden;
        }

        header::after {
            content: '';
            position: absolute;
            top: -50%;
            right: -20%;
            width: 400px;
            height: 400px;
            background: radial-gradient(circle, rgba(99, 102, 241, 0.15) 0%, transparent 70%);
            pointer-events: none;
        }

        .header-title-container {
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
            gap: 1.5rem;
        }

        .header-text h1 {
            font-size: 2.25rem;
            background: linear-gradient(to right, #6366f1, #a5b4fc);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin-bottom: 0.5rem;
        }

        .header-text p {
            color: var(--text-secondary);
            font-size: 1rem;
        }

        .header-meta {
            display: flex;
            gap: 2rem;
        }

        .meta-item {
            background: rgba(0, 0, 0, 0.2);
            border: 1px solid var(--border-color);
            border-radius: 0.75rem;
            padding: 0.75rem 1.25rem;
            text-align: center;
        }

        .meta-val {
            font-family: 'Outfit', sans-serif;
            font-size: 1.25rem;
            font-weight: 700;
            color: #818cf8;
        }

        .meta-lbl {
            font-size: 0.75rem;
            color: var(--text-muted);
            text-transform: uppercase;
            letter-spacing: 0.05em;
            margin-top: 0.25rem;
        }

        /* Dashboard Overview Grid */
        .dashboard-grid {
            display: grid;
            grid-template-columns: 2fr 1fr;
            gap: 2rem;
            margin-bottom: 2.5rem;
        }

        @media (max-width: 1024px) {
            .dashboard-grid {
                grid-template-columns: 1fr;
            }
        }

        .card {
            background: var(--bg-card);
            border: 1px solid var(--border-color);
            backdrop-filter: blur(12px);
            border-radius: 1.25rem;
            padding: 2rem;
            box-shadow: 0 10px 25px rgba(0, 0, 0, 0.2);
        }

        .chart-card {
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
        }

        .chart-container {
            position: relative;
            width: 100%;
            height: 350px;
        }

        .stats-summary {
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 1.5rem;
            height: 100%;
        }

        .stat-box {
            background: rgba(15, 23, 42, 0.5);
            border: 1px solid var(--border-color);
            border-radius: 1rem;
            padding: 1.5rem;
            display: flex;
            flex-direction: column;
            justify-content: center;
            position: relative;
            overflow: hidden;
            transition: transform 0.3s ease;
        }

        .stat-box:hover {
            transform: translateY(-3px);
        }

        .stat-val {
            font-size: 2.5rem;
            font-weight: 800;
            font-family: 'Outfit', sans-serif;
            margin-bottom: 0.25rem;
        }

        .stat-lbl {
            font-size: 0.85rem;
            color: var(--text-secondary);
            font-weight: 500;
        }

        .stat-box.wildtype { color: var(--color-wildtype); border-left: 4px solid var(--color-wildtype); }
        .stat-box.carrier { color: var(--color-carrier); border-left: 4px solid var(--color-carrier); }
        .stat-box.variant { color: var(--color-variant); border-left: 4px solid var(--color-variant); }
        .stat-box.total { color: #818cf8; border-left: 4px solid #818cf8; }

        /* Filter Controls */
        .controls-section {
            background: rgba(30, 41, 59, 0.2);
            border: 1px solid var(--border-color);
            border-radius: 1rem;
            padding: 1.25rem;
            margin-bottom: 2rem;
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 1.5rem;
            flex-wrap: wrap;
        }

        .tabs {
            display: flex;
            gap: 0.5rem;
            flex-wrap: wrap;
        }

        .tab-btn {
            background: transparent;
            border: 1px solid transparent;
            color: var(--text-secondary);
            padding: 0.5rem 1rem;
            border-radius: 0.5rem;
            cursor: pointer;
            font-size: 0.875rem;
            font-weight: 500;
            transition: all 0.2s ease;
        }

        .tab-btn:hover {
            color: var(--text-primary);
            background: rgba(255, 255, 255, 0.05);
        }

        .tab-btn.active {
            color: var(--text-primary);
            background: #4f46e5;
            border-color: #6366f1;
        }

        .search-box {
            position: relative;
            max-width: 320px;
            width: 100%;
        }

        .search-input {
            width: 100%;
            background: rgba(15, 23, 42, 0.7);
            border: 1px solid var(--border-color);
            border-radius: 0.5rem;
            padding: 0.6rem 1rem;
            color: var(--text-primary);
            font-size: 0.875rem;
            outline: none;
            transition: border-color 0.2s ease;
        }

        .search-input:focus {
            border-color: #6366f1;
        }

        /* Matrix Gene Grid */
        .gene-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(320px, 1fr));
            gap: 1.5rem;
        }

        .gene-card {
            background: var(--bg-card);
            border: 1px solid var(--border-color);
            border-radius: 1rem;
            padding: 1.5rem;
            cursor: pointer;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            position: relative;
            overflow: hidden;
            display: flex;
            flex-direction: column;
            justify-content: space-between;
        }

        .gene-card::after {
            content: '';
            position: absolute;
            bottom: 0;
            left: 0;
            width: 100%;
            height: 3px;
            background: transparent;
        }

        .gene-card:hover {
            transform: translateY(-5px);
            background: var(--bg-card-hover);
            box-shadow: 0 15px 30px rgba(0, 0, 0, 0.4);
        }

        .gene-card.wildtype:hover {
            box-shadow: 0 15px 30px var(--glow-wildtype);
        }
        .gene-card.wildtype::after { background: var(--color-wildtype); }

        .gene-card.carrier:hover {
            box-shadow: 0 15px 30px var(--glow-carrier);
        }
        .gene-card.carrier::after { background: var(--color-carrier); }

        .gene-card.variant:hover {
            box-shadow: 0 15px 30px var(--glow-variant);
        }
        .gene-card.variant::after { background: var(--color-variant); }

        .gene-card-header {
            display: flex;
            justify-content: space-between;
            align-items: flex-start;
            margin-bottom: 1rem;
        }

        .gene-name {
            font-size: 1.5rem;
            font-weight: 700;
        }

        .gene-cat {
            font-size: 0.75rem;
            color: var(--text-muted);
            text-transform: uppercase;
            letter-spacing: 0.05em;
            margin-top: 0.15rem;
        }

        /* Risk Gauge */
        .gauge-wrapper {
            position: relative;
            width: 50px;
            height: 50px;
        }

        .gauge-svg {
            transform: rotate(-90deg);
        }

        .gauge-bg {
            fill: none;
            stroke: rgba(255, 255, 255, 0.05);
            stroke-width: 4;
        }

        .gauge-val {
            fill: none;
            stroke-width: 4;
            stroke-linecap: round;
            transition: stroke-dasharray 0.3s ease;
        }

        .gene-card.wildtype .gauge-val { stroke: var(--color-wildtype); }
        .gene-card.carrier .gauge-val { stroke: var(--color-carrier); }
        .gene-card.variant .gauge-val { stroke: var(--color-variant); }

        .gauge-text {
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            font-size: 0.8rem;
            font-weight: 700;
            font-family: 'Outfit', sans-serif;
        }

        .gene-card.wildtype .gauge-text { color: var(--color-wildtype); }
        .gene-card.carrier .gauge-text { color: var(--color-carrier); }
        .gene-card.variant .gauge-text { color: var(--color-variant); }

        .gene-card-body {
            flex-grow: 1;
            margin-bottom: 1rem;
        }

        .gene-trait {
            font-size: 0.95rem;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 0.5rem;
        }

        .gene-desc {
            font-size: 0.825rem;
            color: var(--text-secondary);
            display: -webkit-box;
            -webkit-line-clamp: 3;
            -webkit-box-orient: vertical;
            overflow: hidden;
            text-overflow: ellipsis;
        }

        .gene-card-footer {
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-top: 1px solid rgba(255, 255, 255, 0.05);
            padding-top: 0.75rem;
            margin-top: auto;
        }

        .snps-badge {
            background: rgba(255, 255, 255, 0.05);
            border: 1px solid var(--border-color);
            border-radius: 0.35rem;
            padding: 0.2rem 0.5rem;
            font-size: 0.75rem;
            color: var(--text-secondary);
            font-family: monospace;
        }

        .status-badge {
            font-size: 0.75rem;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }

        .gene-card.wildtype .status-badge { color: var(--color-wildtype); }
        .gene-card.carrier .status-badge { color: var(--color-carrier); }
        .gene-card.variant .status-badge { color: var(--color-variant); }

        /* Modal Structure */
        .modal-overlay {
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(15, 23, 42, 0.85);
            backdrop-filter: blur(8px);
            z-index: 100;
            display: none;
            align-items: center;
            justify-content: center;
            padding: 1.5rem;
        }

        .modal-content {
            background: #111827;
            border: 1px solid var(--border-color);
            border-radius: 1.5rem;
            max-width: 680px;
            width: 100%;
            max-height: 90vh;
            overflow-y: auto;
            box-shadow: 0 25px 50px rgba(0, 0, 0, 0.5);
            animation: modalFadeIn 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        }

        @keyframes modalFadeIn {
            from { opacity: 0; transform: scale(0.95); }
            to { opacity: 1; transform: scale(1); }
        }

        .modal-header {
            padding: 1.75rem 2rem;
            border-bottom: 1px solid rgba(255, 255, 255, 0.05);
            display: flex;
            justify-content: space-between;
            align-items: center;
        }

        .modal-title-left {
            display: flex;
            align-items: center;
            gap: 1.25rem;
        }

        .modal-title-left h2 {
            font-size: 1.75rem;
        }

        .close-btn {
            background: transparent;
            border: none;
            color: var(--text-secondary);
            font-size: 1.5rem;
            cursor: pointer;
            transition: color 0.2s ease;
        }

        .close-btn:hover {
            color: var(--text-primary);
        }

        .modal-body {
            padding: 2rem;
        }

        .section-title {
            font-size: 0.9rem;
            color: var(--text-muted);
            text-transform: uppercase;
            letter-spacing: 0.05em;
            margin-bottom: 0.5rem;
        }

        .modal-trait {
            font-size: 1.1rem;
            font-weight: 700;
            color: var(--text-primary);
            margin-bottom: 1.5rem;
        }

        .implications-block {
            background: rgba(30, 41, 59, 0.35);
            border: 1px solid var(--border-color);
            border-radius: 0.75rem;
            padding: 1.25rem 1.5rem;
            margin-bottom: 2rem;
        }

        .implications-block p {
            font-size: 0.95rem;
            color: var(--text-primary);
        }

        .snp-list {
            display: flex;
            flex-direction: column;
            gap: 1.25rem;
        }

        .snp-item {
            border: 1px solid var(--border-color);
            border-radius: 0.75rem;
            padding: 1.25rem;
            background: rgba(15, 23, 42, 0.4);
        }

        .snp-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 0.75rem;
        }

        .snp-id {
            font-weight: 700;
            font-size: 1rem;
            color: #818cf8;
        }

        .snp-status {
            font-size: 0.75rem;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            padding: 0.2rem 0.5rem;
            border-radius: 0.25rem;
            background: rgba(255, 255, 255, 0.03);
            border: 1px solid var(--border-color);
        }

        .snp-status.wildtype { color: var(--color-wildtype); }
        .snp-status.carrier { color: var(--color-carrier); }
        .snp-status.variant { color: var(--color-variant); }

        .snp-meta-grid {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 1rem;
            margin-bottom: 0.75rem;
        }

        .snp-meta-item {
            font-size: 0.8rem;
            color: var(--text-secondary);
        }

        .snp-meta-lbl {
            color: var(--text-muted);
            font-size: 0.7rem;
            text-transform: uppercase;
            margin-bottom: 0.15rem;
        }

        .snp-meta-val {
            font-weight: 500;
        }

        .snp-desc {
            font-size: 0.825rem;
            color: var(--text-secondary);
            border-top: 1px dashed rgba(255, 255, 255, 0.05);
            padding-top: 0.5rem;
            margin-top: 0.5rem;
        }
    </style>
</head>
<body>

    <!-- Header -->
    <header>
        <div class="header-title-container">
            <div class="header-text">
                <h1>Anonymous Patient</h1>
                <p>Personal Genomic Health & Variant-Calling Analysis</p>
                <div style="margin-top: 1rem; display: flex; gap: 0.75rem; align-items: center;">
                    <a href="https://github.com/fredLuv/fastq-genomic-pipeline" target="_blank" style="display: inline-flex; align-items: center; gap: 0.5rem; background: rgba(255,255,255,0.08); border: 1px solid var(--border-color); color: var(--text-primary); text-decoration: none; padding: 0.4rem 0.8rem; border-radius: 0.5rem; font-size: 0.8rem; font-weight: 500; transition: background 0.2s ease;">
                        <svg height="16" width="16" viewBox="0 0 16 16" fill="currentColor" style="vertical-align: middle;">
                            <path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path>
                        </svg>
                        View on GitHub
                    </a>
                </div>
            </div>
            <div class="header-meta">
                <div class="meta-item">
                    <div class="meta-val">631,455</div>
                    <div class="meta-lbl">SNPs Genotyped</div>
                </div>
                <div class="meta-item">
                    <div class="meta-val">GRCh38</div>
                    <div class="meta-lbl">Genome Build</div>
                </div>
            </div>
        </div>
    </header>

    <!-- Dashboard -->
    <div class="dashboard-grid">
        <!-- Statistics Summary -->
        <div class="card stats-card">
            <h3 style="margin-bottom: 1.25rem; font-size: 1.25rem; color: var(--text-secondary);">Genetic Risk Profiles</h3>
            <div class="stats-summary">
                <div class="stat-box total">
                    <div class="stat-val" id="stat-total">0</div>
                    <div class="stat-lbl">Genes Analyzed</div>
                </div>
                <div class="stat-box variant">
                    <div class="stat-val" id="stat-high">0</div>
                    <div class="stat-lbl">High Risk Genes (Risk = 1.0)</div>
                </div>
                <div class="stat-box carrier">
                    <div class="stat-val" id="stat-carrier">0</div>
                    <div class="stat-lbl">Carrier Genes (Risk = 0.5)</div>
                </div>
                <div class="stat-box wildtype">
                    <div class="stat-val" id="stat-wildtype">0</div>
                    <div class="stat-lbl">Wildtype Genes (Risk = 0.0)</div>
                </div>
            </div>
        </div>
        
        <!-- Radar Category Chart -->
        <div class="card chart-card">
            <h3 style="margin-bottom: 1.25rem; font-size: 1.25rem; color: var(--text-secondary); width: 100%; text-align: left;">Category Risk Assessment</h3>
            <div class="chart-container">
                <canvas id="categoryChart"></canvas>
            </div>
        </div>
    </div>

    <!-- Filtering & Controls -->
    <div class="controls-section">
        <div class="tabs" id="categoryTabs">
            <button class="tab-btn active" data-category="All">All Categories</button>
            <button class="tab-btn" data-category="Methylation">Methylation</button>
            <button class="tab-btn" data-category="Brain Health">Brain Health</button>
            <button class="tab-btn" data-category="Cardiovascular">Cardiovascular</button>
            <button class="tab-btn" data-category="Dietary">Dietary Tolerance</button>
            <button class="tab-btn" data-category="Vitamins">Vitamins</button>
            <button class="tab-btn" data-category="Detox">Detox & AntiOx</button>
            <button class="tab-btn" data-category="Inflammation">Inflammation</button>
            <button class="tab-btn" data-category="Physical Traits">Physical Traits</button>
        </div>
        <div class="search-box">
            <input type="text" class="search-input" id="searchInput" placeholder="Search genes (e.g. MTHFR)...">
        </div>
    </div>

    <!-- Gene Matrix Grid -->
    <div class="gene-grid" id="geneGrid"></div>

    <!-- Gene Details Modal -->
    <div class="modal-overlay" id="detailsModal">
        <div class="modal-content">
            <div class="modal-header">
                <div class="modal-title-left">
                    <h2 id="modalGeneName">MTHFR</h2>
                    <div class="snps-badge" id="modalGeneChr">Chr 1</div>
                </div>
                <button class="close-btn" id="modalClose">&times;</button>
            </div>
            <div class="modal-body">
                <div class="section-title">Genetic Trait</div>
                <div class="modal-trait" id="modalGeneTrait">Folate Metabolism</div>
                
                <div class="section-title">Clinical Implications</div>
                <div class="implications-block">
                    <p id="modalGeneImplication">Reduced folate conversion...</p>
                </div>
                
                <div class="section-title" style="margin-bottom: 0.75rem;">Measured Genotypes (SNPs)</div>
                <div class="snp-list" id="modalSnpList"></div>
            </div>
        </div>
    </div>

    <script>
        // Inject compiled genomic data directly
        const genomicData = %s;

        let activeCategory = "All";
        let searchQuery = "";

        // Elements
        const geneGrid = document.getElementById("geneGrid");
        const categoryTabs = document.getElementById("categoryTabs");
        const searchInput = document.getElementById("searchInput");
        
        // Modal elements
        const detailsModal = document.getElementById("detailsModal");
        const modalGeneName = document.getElementById("modalGeneName");
        const modalGeneChr = document.getElementById("modalGeneChr");
        const modalGeneTrait = document.getElementById("modalGeneTrait");
        const modalGeneImplication = document.getElementById("modalGeneImplication");
        const modalSnpList = document.getElementById("modalSnpList");
        const modalClose = document.getElementById("modalClose");

        // Initialization
        document.addEventListener("DOMContentLoaded", () => {
            renderOverviewStats();
            renderCategoryChart();
            renderGrid();
            
            // Tab click handling
            categoryTabs.addEventListener("click", (e) => {
                if (e.target.classList.contains("tab-btn")) {
                    document.querySelectorAll(".tab-btn").forEach(btn => btn.classList.remove("active"));
                    e.target.classList.add("active");
                    activeCategory = e.target.getAttribute("data-category");
                    renderGrid();
                }
            });

            // Search typing handling
            searchInput.addEventListener("input", (e) => {
                searchQuery = e.target.value.toLowerCase().trim();
                renderGrid();
            });

            // Modal Close
            modalClose.addEventListener("click", () => { detailsModal.style.display = "none"; });
            detailsModal.addEventListener("click", (e) => {
                if (e.target === detailsModal) detailsModal.style.display = "none";
            });
        });

        // Compute and Display Overview Stats
        function renderOverviewStats() {
            let total = genomicData.length;
            let high = 0;
            let carrier = 0;
            let wildtype = 0;

            genomicData.forEach(g => {
                if (g.risk_score === 1.0) high++;
                else if (g.risk_score > 0.0) carrier++;
                else wildtype++;
            });

            document.getElementById("stat-total").innerText = total;
            document.getElementById("stat-high").innerText = high;
            document.getElementById("stat-carrier").innerText = carrier;
            document.getElementById("stat-wildtype").innerText = wildtype;
        }

        // Draw Category Risk Chart
        function renderCategoryChart() {
            // Group and compute average risk score per category
            const categories = {};
            genomicData.forEach(g => {
                let cat = g.category;
                if (cat.includes("Brain")) cat = "Brain Health";
                if (cat.includes("Vitamins") || cat.includes("Vitamin")) cat = "Vitamins";
                if (cat.includes("Detox") || cat.includes("Antioxidant")) cat = "Detox";
                
                if (!categories[cat]) categories[cat] = { totalScore: 0, count: 0 };
                categories[cat].totalScore += g.risk_score;
                categories[cat].count++;
            });

            const labels = [];
            const dataVals = [];

            for (const [cat, val] of Object.entries(categories)) {
                labels.push(cat);
                dataVals.push(parseFloat((val.totalScore / val.count).toFixed(2)));
            }

            const ctx = document.getElementById('categoryChart').getContext('2d');
            new Chart(ctx, {
                type: 'radar',
                data: {
                    labels: labels,
                    datasets: [{
                        label: 'Average Genetic Load',
                        data: dataVals,
                        backgroundColor: 'rgba(99, 102, 241, 0.2)',
                        borderColor: '#6366f1',
                        borderWidth: 2,
                        pointBackgroundColor: '#818cf8',
                        pointBorderColor: '#fff',
                        pointRadius: 4
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {
                        r: {
                            angleLines: { color: 'rgba(255, 255, 255, 0.08)' },
                            grid: { color: 'rgba(255, 255, 255, 0.08)' },
                            pointLabels: { color: '#94a3b8', font: { family: 'Outfit', size: 11 } },
                            ticks: { display: false },
                            min: 0,
                            max: 1.0
                        }
                    },
                    plugins: {
                        legend: { display: false }
                    }
                }
            });
        }

        // Render Gene Cards Grid
        function renderGrid() {
            geneGrid.innerHTML = "";

            const filteredGenes = genomicData.filter(g => {
                let catMatch = false;
                if (activeCategory === "All") {
                    catMatch = true;
                } else {
                    const searchCat = activeCategory.toLowerCase();
                    const geneCat = g.category.toLowerCase();
                    if (searchCat === "detox" && (geneCat.includes("detox") || geneCat.includes("antioxidant"))) {
                        catMatch = true;
                    } else if (searchCat === "dietary" && geneCat.includes("diet")) {
                        catMatch = true;
                    } else {
                        catMatch = geneCat.includes(searchCat);
                    }
                }

                const nameMatch = g.name.toLowerCase().includes(searchQuery) || 
                                  g.trait.toLowerCase().includes(searchQuery);

                return catMatch && nameMatch;
            });

            if (filteredGenes.length === 0) {
                geneGrid.innerHTML = `<div style="grid-column: 1/-1; text-align: center; padding: 3rem; color: var(--text-muted);">No genetic markers match the current filter criteria.</div>`;
                return;
            }

            filteredGenes.forEach(g => {
                let statusClass = "wildtype";
                let statusText = "Homozygous Wildtype";
                if (g.risk_score === 1.0) {
                    statusClass = "variant";
                    statusText = "Homozygous Risk";
                } else if (g.risk_score > 0.0) {
                    statusClass = "carrier";
                    statusText = "Heterozygous Risk";
                }

                const radius = 22;
                const circumference = 2 * Math.PI * radius;
                const offset = circumference - (g.risk_score * circumference);

                const card = document.createElement("div");
                card.className = `gene-card ${statusClass}`;
                card.innerHTML = `
                    <div class="gene-card-header">
                        <div>
                            <div class="gene-name">${g.name}</div>
                            <div class="gene-cat">${g.category}</div>
                        </div>
                        <div class="gauge-wrapper">
                            <svg width="50" height="50" class="gauge-svg">
                                <circle class="gauge-bg" cx="25" cy="25" r="${radius}"></circle>
                                <circle class="gauge-val" cx="25" cy="25" r="${radius}" 
                                        stroke-dasharray="${circumference}" 
                                        stroke-dashoffset="${offset}"></circle>
                            </svg>
                            <div class="gauge-text">${g.risk_score.toFixed(1)}</div>
                        </div>
                    </div>
                    <div class="gene-card-body">
                        <div class="gene-trait">${g.trait}</div>
                        <div class="gene-desc">${g.implication}</div>
                    </div>
                    <div class="gene-card-footer">
                        <span class="snps-badge">${g.snps.length} SNP${g.snps.length > 1 ? 's' : ''}</span>
                        <span class="status-badge">${statusText}</span>
                    </div>
                `;

                card.addEventListener("click", () => showGeneDetails(g));
                geneGrid.appendChild(card);
            });
        }

        // Show Modal Dialog on Gene Card Click
        function showGeneDetails(g) {
            modalGeneName.innerText = g.name;
            modalGeneChr.innerText = `Chr ${g.chromosome}`;
            modalGeneTrait.innerText = g.trait;
            modalGeneImplication.innerHTML = g.implication; // Allow HTML tags in modal

            modalSnpList.innerHTML = "";
            g.snps.forEach(s => {
                let statusClass = "wildtype";
                if (s.Risk_Count === 2) statusClass = "variant";
                else if (s.Risk_Count === 1) statusClass = "carrier";

                const snpItem = document.createElement("div");
                snpItem.className = "snp-item";
                snpItem.innerHTML = `
                    <div class="snp-header">
                        <span class="snp-id">${s.rsID}</span>
                        <span class="snp-status ${statusClass}">${s.Status}</span>
                    </div>
                    <div class="snp-meta-grid">
                        <div class="snp-meta-item">
                            <div class="snp-meta-lbl">Genotype</div>
                            <div class="snp-meta-val">${s.Genotype || '--'}</div>
                        </div>
                        <div class="snp-meta-item">
                            <div class="snp-meta-lbl">Risk Allele</div>
                            <div class="snp-meta-val">${s.Risk_Allele}</div>
                        </div>
                        <div class="snp-meta-item">
                            <div class="snp-meta-lbl">Risk Count</div>
                            <div class="snp-meta-val">${s.Risk_Count}</div>
                        </div>
                    </div>
                    <div class="snp-meta-grid" style="grid-template-columns: 2fr 1fr; margin-bottom: 0;">
                        <div class="snp-meta-item">
                            <div class="snp-meta-lbl">Consequence</div>
                            <div class="snp-meta-val" style="font-family: monospace; font-size: 0.75rem;">${s.Consequence}</div>
                        </div>
                        <div class="snp-meta-item">
                            <div class="snp-meta-lbl">Position (GRCh38)</div>
                            <div class="snp-meta-val" style="font-family: monospace; font-size: 0.75rem;">${s.Position}</div>
                        </div>
                    </div>
                    <div class="snp-desc">
                        <strong>Associated Trait</strong>: ${s.Trait}<br>
                        <strong>Implications</strong>: ${s.interpretation}
                    </div>
                `;
                modalSnpList.appendChild(snpItem);
            });

            detailsModal.style.display = "flex";
        }
    </script>
</body>
</html>
"""

# Format HTML using simple string replace
json_data = json.dumps(genes, indent=4)
final_html_content = html_template.replace("%s", json_data)

# Write self-contained HTML file
with open(html_output_path, 'w', encoding='utf-8') as f:
    f.write(final_html_content)

print(f"Interactive genomic HTML report successfully written to {html_output_path}")
