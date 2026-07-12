import os
import csv
import json

# Setup directories
WORKDIR = "/Users/fred/Library/CloudStorage/GoogleDrive-fredlu0521@gmail.com/My Drive/Ethan-Resume/fastq_pipeline"
RESULTS_DIR = os.path.join(WORKDIR, "results")

gene_csv_path = os.path.join(RESULTS_DIR, "generic_gene_features.csv")
snp_csv_path = os.path.join(RESULTS_DIR, "generic_snp_features.csv")
html_output_path = os.path.join(RESULTS_DIR, "genomic_report.html")

# Define Patient-Specific Clinical Interpretations based on rsID, Risk Count, and Genotype
def get_patient_interpretation(rsid, risk_count, genotype):
    # Normalize genotype
    genotype = (genotype or "--").upper().strip()
    
    if rsid == "rs1801133":
        if risk_count == 2:
            return "<b>叶酸代谢能力重度降低 (TT)</b>：MTHFR 酶活性较常人降低约 84%，导致同型半胱氨酸升高风险显著增加。建议补充活性甲基叶酸 (5-MTHF)。"
        elif risk_count == 1:
            return "<b>叶酸代谢能力轻度降低 (CT)</b>：MTHFR 酶活性降低约 35%。建议日常膳食中增加绿叶蔬菜，适当补充活性叶酸。"
        else:
            return "<b>叶酸代谢正常 (CC)</b>：未发现 MTHFR C677T 突变，具有典型的叶酸转化及甲基化代谢速度。"
            
    elif rsid == "rs1801131":
        if risk_count == 2:
            return "<b>叶酸代谢中度降低 (CC)</b>：MTHFR 酶活性降低约 40%，伴随中度的甲基供体不足倾向。"
        elif risk_count == 1:
            return "<b>叶酸代谢正常/杂合携带者 (GT/AC)</b>：典型的 MTHFR 酶活性，体内甲基化处于平衡水平。"
        else:
            return "<b>叶酸代谢正常 (AA)</b>：未检出 A1298C 突变，具备典型且高效的叶酸利用能力。"
            
    elif rsid == "rs4680":
        if risk_count == 0:
            return "<b>战士基因型 (Warrior - GG)</b>：多巴胺降解速率极快。在高压环境中冷静沉着、专注度高，但在低压日常环境中容易觉得无聊、缺乏动力。"
        elif risk_count == 1:
            return "<b>平衡基因型 (AG)</b>：多巴胺降解速率中等。在日常学习工作与应对突发压力之间具有健康的平衡力。"
        else:
            return "<b>谋士基因型 (Worrier - AA)</b>：多巴胺降解极为缓慢。在舒适冷静的环境中拥有极强的认知深度、专注力与记忆力，但在高压环境下极易感到焦虑和敏感。"
            
    elif rsid == "rs6265":
        if risk_count >= 1:
            return "<b>神经营养因子分泌降低 (AG/AA)</b>：Met 突变型会导致活动依赖型 BDNF 分泌轻度下调，略微影响神经突触可塑性。建议通过规律有氧运动有效代偿并刺激 BDNF 分泌。"
        else:
            return "<b>神经营养因子正常 (CC)</b>：典型的 BDNF 分泌活性与大脑神经突触形成能力，认知保护与神经细胞修护机制正常。"
            
    elif rsid == "rs1800497":
        if risk_count == 1:
            return "<b>多巴胺 D2 受体密度降低 (AG)</b>：Taq1A 突变携带者。大脑多巴胺 D2 受体数量较典型值减少约 30%，容易表现出更高的奖赏阈值（较易沉迷甜食、刺激或探索行为）。"
        elif risk_count == 2:
            return "<b>多巴胺 D2 受体密度显著降低 (AA)</b>：多巴胺受体数量减少超 30%。天然满足感较低，容易陷入冲动行为或寻求外界强刺激。建议建立良性生活自律习惯。"
        else:
            return "<b>多巴胺受体正常 (GG)</b>：具有典型的多巴胺 D2 受体密度与多巴胺通路响应机制，内在奖赏系统处于平衡。"
            
    elif rsid == "rs6318":
        return "<b>单胺氧化酶活性正常</b>：典型的神经递质（血清素、多巴胺）清除降解速率，情绪耐受力正常。"
        
    elif rsid == "rs429358":
        if "C" in genotype:
            return "<b>检测到 APOE-e4 风险等位基因</b>：阿尔茨海默病以及低密度脂蛋白胆固醇升高的遗传风险增加。日常建议采用地中海饮食并定期复查血脂。"
        else:
            return "<b>未检测到 APOE-e4 风险基因 (TT)</b>：该位点表现为健康的典型单倍型背景。"
            
    elif rsid == "rs7412":
        if "T" in genotype:
            return "<b>检测到 APOE-e2 保护性等位基因</b>：与降低阿尔茨海默病发病风险以及改善脂质清除相关。"
        else:
            return "<b>典型的 APOE-e3 正常背景 (CC)</b>：未发现该位点突变，脂质代谢基准正常。"
            
    elif rsid == "rs1799983":
        if risk_count >= 1:
            return "<b>一氧化氮生成受限 (TT/GT)</b>：血管内皮型一氧化氮合成酶活性略微降低，可能稍微影响微循环血管舒张。日常可多补充甜菜根等一氧化氮食物来源。"
        else:
            return "<b>一氧化氮生成典型 (GG)</b>：心血管内皮血管收缩与舒张平衡正常，具有健康的微循环机能。"
            
    elif rsid == "rs2802292":
        if "G" in genotype or "C" in genotype:
            return "<b>长寿相关 FOXO3 保护基因携带者 (GT)</b>：能有效激活细胞自噬与抗氧化酶基因的转录，增强心血管血管内皮的健康老化保护能力。"
        else:
            return "<b>典型抗衰老能力 (TT)</b>：细胞抗氧化酶系统基底功能正常，老化速度处于平均状态。"
            
    elif rsid == "rs4988235":
        if risk_count == 2:
            return "<b>重度乳糖不耐受 (GG)</b>：小肠粘膜乳糖酶发生年龄相关性完全停产。饮用普通牛奶或乳制品后极易发生腹胀、腹泻或肠易激，日常建议选用无乳糖奶或植物奶。"
        else:
            return "<b>乳糖耐受 (AA/AG)</b>：乳糖酶在成年期持续产生，能毫无障碍地消化和分解动物奶中的乳糖。"
            
    elif rsid == "rs9939609":
        if risk_count == 2:
            return "<b>肥胖与食欲失控高风险 (AA)</b>：饥饿素（Ghrelin）抑制迟缓。饱腹感信号延迟，极易无意识摄入额外卡路里。建议执行高饱腹感、高蛋白质膳食策略。"
        elif risk_count == 1:
            return "<b>中度肥胖与能量代谢风险 (TA)</b>：食欲控制及脂肪囤积倾向中等，建议适当控制精制碳水摄入。"
        else:
            return "<b>正常能量代谢 (TT)</b>：正常的饱腹感控制反应，未表现出明显的易胖体质遗传倾向。"
            
    elif rsid == "rs1801282":
        if risk_count >= 1:
            return "<b>胰岛素敏感度优化 variant (CG)</b>：虽然与略微升高的 BMI 敏感度相关，但在防范糖尿病及保护血管内皮胰岛素敏感性方面具备遗传优势。"
        else:
            return "<b>典型 PPARG 胰岛素敏感性 (CC)</b>：处于人群均值水平的脂肪分化与胰岛素结合度。"
            
    elif rsid == "rs1801260":
        if risk_count >= 1:
            return "<b>夜猫子倾向 (AG/GG)</b>：昼夜节律 CLOCK 基因轻微位移，表现为典型的夜间型 chronotype，入睡时间偏晚且清晨精力恢复较慢。"
        else:
            return "<b>晨型人（云雀）基因 (AA)</b>：拥有规律且高度自洽的白昼节律系统，晨间精力充沛。"
            
    elif rsid == "rs1800562":
        if risk_count == 2:
            return "<b>遗传性铁过载高危 (AA)</b>：血色病高风险（Hemochromatosis），肠道铁吸收失控。需定期监测铁蛋白 (Ferritin) 与转铁蛋白饱和度。"
        elif risk_count == 1:
            return "<b>血色病基因携带者 (GA)</b>：通常不会引发任何铁沉积症状，体内铁负荷完全正常。"
        else:
            return "<b>典型铁吸收速度 (GG)</b>：体内铁吸收与代谢负荷受健康的 Hepcidin 系统严格调控，无积铁危险。"
            
    elif rsid == "rs1799945":
        if risk_count >= 1:
            return "<b>血色病协同携带 (H63D - CG)</b>：肠道铁吸收无显著异常。除非与 C282Y 双重杂合，否则铁蛋白水平正常。"
        else:
            return "<b>典型铁通道功能 (CC)</b>：正常的血色素与铁离子吸收机制。"
            
    elif rsid == "rs731236":
        if risk_count >= 1:
            return "<b>维生素 D 受体效率降低 (tt/Gt)</b>：肠道对钙和活性维生素 D3 的吸收与骨骼结合效率略低。建议定期监测维 D 浓度，必要时补充维 D3 和 K2。"
        else:
            return "<b>维生素 D 结合正常 (AA)</b>：VDR 结合性能极佳，拥有健康的骨矿化与骨骼密度保护力。"
            
    elif rsid == "rs7501331":
        if risk_count == 2:
            return "<b>胡萝卜素转化率重度降低 (CC)</b>：体内 BCMO1 酶活性降低达 57%。从植物（胡萝卜、红薯）中转化活性维 A (Retinol) 的能力极差。日常必需直接摄入蛋黄、动物肝脏以保证维 A 充足。"
        elif risk_count == 1:
            return "<b>胡萝卜素转化率轻度降低 (TC)</b>：维 A 转化能力降低约 32%。"
        else:
            return "<b>典型的维 A 转化速度 (TT)</b>：能将食物中的 $\beta$-胡萝卜素高效裂解转化为活性维生素 A。"
            
    elif rsid == "rs602662":
        if risk_count == 2:
            return "<b>维生素 B12 吸收能力重度降低 (AA)</b>：胃肠道对 B12 的特异性结合吸收偏慢。极易发生 B12 不足，建议补充活性甲钴胺。"
        elif risk_count == 1:
            return "<b>维生素 B12 吸收轻度受限 (AG)</b>：体内活性 B12 蓄积率处于轻微偏低范围。"
        else:
            return "<b>典型 B12 吸收率 (GG)</b>：具有高水准的 B12 生物利用度与胃粘膜吸收转化速度。"
            
    elif rsid == "rs4654748":
        if risk_count >= 1:
            return "<b>维生素 B6 清除速率加快 (CT/TT)</b>：由于运输载体表达差异，体内的 B6 在血液中降解清除偏快。建议适当高频获取富含 B6 的膳食。"
        else:
            return "<b>典型 B6 留存率 (CC)</b>：血浆中维生素 B6 的利用和存储时间稳定。"
            
    elif rsid == "rs1695":
        if risk_count == 2:
            return "<b>谷胱甘肽解毒速率重度降低 (GG)</b>：GSTP1 纯合突变导致细胞清除重金属、空气污染微粒和环境雌激素毒素的速度显著降低。建议多吃十字花科蔬菜提升解毒储备。"
        elif risk_count == 1:
            return "<b>谷胱甘肽解毒速率轻度降低 (AG)</b>：细胞抗毒物与解毒酶活性中等。"
        else:
            return "<b>谷胱甘肽解毒正常 (AA)</b>：高效且典型的细胞第一/二阶段解毒排毒循环系统。"
            
    elif rsid == "rs4880":
        if risk_count == 2:
            return "<b>线粒体抗氧化力偏低 (SOD2 TT)</b>：锰型超氧化物歧化酶向线粒体的转移受阻。线粒体抗自由基损伤力降低，需要增加饮食中的强抗氧化剂（茶多酚、浆果）。"
        else:
            return "<b>线粒体防线正常 (GG/GT)</b>：充足的线粒体内部自由基清除能力，细胞能量工厂衰老保护力强。"
            
    elif rsid == "rs671":
        if risk_count == 2:
            return "<b>重度酒精红脸/乙醛毒性 (AA)</b>：乙醛脱氢酶 (ALDH2) 活性接近为 0。即使摄入极微量酒精也会导致乙醛迅速蓄积并引起面红、剧烈头痛、心悸等，严禁饮酒。"
        elif risk_count == 1:
            return "<b>轻中度酒精脸红反应 (GA)</b>：解酒酶速度减缓约 60-80%，酒后易快速脸红，宿醉感强烈。建议控制酒精摄入。"
        else:
            return "<b>解酒能力完全正常 (GG)</b>：体内 ALDH2 酶活性充沛，乙醇和乙醛转化排泄速度典型，无红脸反应。"
            
    elif rsid == "rs1800795":
        if risk_count == 2:
            return "<b>慢性炎症因子基线升高 (GG)</b>：促炎因子白介素-6 基因处于活跃表达状态。机体极易进入慢性低度炎性负担状态（Low-grade chronic inflammation），日常宜遵循抗炎饮食模式。"
        else:
            return "<b>炎性细胞因子典型 (CC/CG)</b>：典型的促炎及抑炎反应通路，机体炎性基准处于平均水平。"
            
    elif rsid == "rs1800629":
        if risk_count >= 1:
            return "<b>促炎反应敏感性升高 (AG/AA)</b>：肿瘤坏死因子-alpha（TNF-$\alpha$）反应偏高。受到细菌或毒素刺激后免疫炎性风暴启动更快。推荐补充 Omega-3 脂肪酸。"
        else:
            return "<b>炎症细胞因子反应稳定 (GG)</b>：TNF-alpha 的基线生成完全正常，自身免疫反应平稳。"
            
    elif rsid == "rs1815739":
        if risk_count == 2:
            return "<b>纯耐力型体质 (TT)</b>：快肌纤维中彻底缺失 $\alpha$-actinin-3 爆发力蛋白。天然缺乏短跑、瞬间爆发力项目的优势，但**在马拉松、越野、长距离耐力运动中表现出卓越的肌肉代谢和低耗氧效率**。"
        elif risk_count == 1:
            return "<b>全能平衡型肌纤维 (CT)</b>：肌纤维兼具爆发力蛋白质组与有氧代谢的折中配置，适合各类复合型体育活动。"
        else:
            return "<b>爆发/短跑型体质 (CC)</b>：快肌纤维富含 actinin-3 爆发蛋白。拥有优异的快肌收缩速度，适合极瞬爆发力项目（冲刺、举重）。"
            
    elif rsid == "rs1805007":
        if risk_count == 2:
            return "<b>红发/极度白肤基因 (TT)</b>：黑色素受体 MC1R 功能完全失活，导致真黑色素极其缺乏。皮肤极易晒伤而几乎无法晒黑，罹患皮肤光损伤的概率极高。"
        elif risk_count == 1:
            return "<b>红发基因携带者 (CT)</b>：正常的黑色素分泌及典型防晒伤皮层状态。"
        else:
            return "<b>典型黑色素生成 (CC)</b>：真黑色素分泌正常，皮肤能健康的适应紫外线并生成防晒古铜色。"
            
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
        
        # Calculate patient-specific interpretation
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
        
        # Build consolidated clinical implication for the Gene
        if len(snps) == 1:
            implication = snps[0]["interpretation"]
            trait = snps[0]["Trait"]
        elif len(snps) > 1:
            # Filter out non-genotyped and combine descriptions of genotyped SNPs
            implications = []
            for s in snps:
                if s["Genotype"] != "--":
                    # Format as: "【rsID】Interpretation"
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
    <title>Lorena Sandoval - Personal Genomic Health Report</title>
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
                <h1>Lorena Sandoval</h1>
                <p>Personal Genomic Health & Variant-Calling Analysis</p>
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
