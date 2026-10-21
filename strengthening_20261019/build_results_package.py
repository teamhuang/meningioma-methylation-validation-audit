from pathlib import Path
import json,hashlib,zipfile
import numpy as np,pandas as pd
from run_strengthening import OUT,DEST,BASE,CORE

NAMES={
'baylor_published':'Baylor 已发表 MenG 分组',
'baylor_transferred':'Baylor 研究用质心模型',
'tubingen_transferred':'Tübingen 研究用质心模型',
'tubingen_published':'Tübingen 作者 Toronto 分组 两指示变量',
'ucsf_published':'UCSF 已发表分组',
'tubingen_four_groups':'Tübingen 完整四分组',
'tubingen_categorical_grade':'Tübingen WHO分类编码',
'ucsf_categorical_grade':'UCSF WHO分类编码',
'ucsf_mib1':'UCSF 临床加MIB1',
'ucsf_groups_beyond_mib1':'UCSF 临床加MIB1再加分组',
'ucsf_positive_followup':'UCSF 排除零时点事件',
'baylor_published_positive_followup':'Baylor 排除零随访',
'baylor_author_knn_beta_adapted':'Baylor 作者KNN 已处理β值适配',
'baylor_author_knn_beta_highconfidence':'Baylor 作者KNN 高置信子集'}
s=pd.read_csv(OUT/'validation_summary.csv').set_index('comparison')
lr=pd.read_csv(OUT/'unpenalized_association_tests.csv').set_index('comparison')
risk=pd.read_csv(OUT/'oof_risk_summary.csv')

def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(x) for x in row)+' |' for row in rows])
def f(x):return f'{x:.3f}'
def sign(x):return f'{x:+.3f}'
allrows=[]
for key in CORE+[k for k in s.index if k not in CORE]:
    z=s.loc[key]
    allrows.append([NAMES[key],f'{int(z.n)}/{int(z.events)}',f(z.C_app_clinical)+' → '+f(z.C_app_molecular),
      f(z.cv_C_clinical)+' → '+f(z.cv_C_molecular),sign(z.cv_delta_mean),
      sign(z.cv_delta_split_p025)+' 至 '+sign(z.cv_delta_split_p975),
      sign(z.delta_corrected) if pd.notna(z.delta_corrected) else '未计算',f'{int(z.cv_positive_repeats)}/20'])
alltable=table(['比较','记录/事件','表观C','重复交叉验证C','平均折外ΔC','折分2.5至97.5百分位','乐观校正ΔC','正向重复'],allrows)
lrtrows=[]
for key in CORE:
    z=lr.loc[key]
    lrtrows.append([NAMES[key],f(z.LR),int(z['df']),f'{z.P:.6f}',f'{z.P_Holm_core_family:.6f}',f'{z.P_Holm_all_tested:.6f}'])
lrtable=table(['比较','全局LR统计量','自由度','原始P','五核心比较Holm P','全部14比较Holm P'],lrtrows)
rr=[]
for key in CORE+['ucsf_groups_beyond_mib1','baylor_author_knn_beta_adapted','baylor_author_knn_beta_highconfidence']:
    z=risk[(risk.comparison==key)&risk.horizon.isin([5.,60.])].iloc[0]
    rr.append([NAMES[key],f(z.OP_molecular_mean),f(z.brier_clinical_mean),f(z.brier_molecular_mean),sign(z.delta_brier_mean)])
risktable=table(['比较','折外O/P平均值','基线Brier','加分子层Brier','ΔBrier 越小越好'],rr)
refs='''1. Bayley, J. C. et al. Multiple approaches converge on three biological subtypes of meningioma and extract new insights from published studies. *Sci. Adv.* **8**, eabm6247 (2022). https://doi.org/10.1126/sciadv.abm6247
2. Duckett, D. et al. Development, validation, and utility of a clinically applicable methylation classifier for recurrence risk prediction in meningiomas. *Acta Neuropathol. Commun.* **13**, 259 (2025). https://doi.org/10.1186/s40478-025-02184-9
3. Choudhury, A. et al. Meningioma DNA methylation groups identify biological drivers and therapeutic vulnerabilities. *Nat. Genet.* **54**, 649–659 (2022). https://doi.org/10.1038/s41588-022-01061-8
4. Wang, J. Z. et al. Molecular classification to refine surgical and radiotherapeutic decision-making in meningioma. *Nat. Med.* **30**, 3173–3183 (2024). https://doi.org/10.1038/s41591-024-03167-4
5. Merk, D. J. et al. DNA methylation profiling identifies long-range epigenetic silencing of clustered protocadherins as a key determinant of meningioma progression. *Nat. Commun.* **17**, 9236 (2026). https://doi.org/10.1038/s41467-026-77170-3
6. Collins, G. S. et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ* **385**, e078378 (2024). https://doi.org/10.1136/bmj-2023-078378
7. Riley, R. D. et al. Evaluation of clinical prediction models (part 2): how to undertake an external validation study. *BMJ* **384**, e074820 (2024). https://doi.org/10.1136/bmj-2023-074820
'''
report=f'''# npj Digital Medicine 数据与统计补强报告

2026年10月19日。本轮未制作图，未覆盖原稿、原始数据或既有结果。完成14个比较、每个比较20次重复五折验证，6个比较各1000次配对bootstrap；原始数值、折分、逐例折外预测、拟合失败记录均保留。

## 当前结论

最稳健的证据是UCSF已发表分组在其临床记录中的增量：折外C从0.710升至0.741，ΔC +0.032；乐观校正ΔC +0.033；全局关联检验经全部14比较Holm校正仍支持关联。加入MIB1后的额外ΔC +0.012，但此增补检验经全部比较校正后不达到0.05。Baylor和Tübingen存在正向增量，精度、患者身份与模型泛化仍限制推断。不能将本研究描述为全部阳性、完整复现原分类器或已经证明临床效用。

## 一 核心结果和全部敏感性结果

{alltable}

两模型在每个比较中使用同一完整病例集与同一折分。C为同一测试折内可比较患者对累加后的值，再对20次重复取平均。折分百分位仅描述随机划分变异，**不是总体95%置信区间，也不等于显著性检验**。6个乐观校正比较全部完成1000次有效重拟合，1000次并不构成校正后ΔC的抽样置信区间。各规格完整表见 `results/validation_summary.csv`。

## 二 已修复的实质性问题

### 原论文分类器与本项目质心模型不同

作者公开模型为141例平衡训练参考、403个CpG、K=10的KNN；本项目旧模型使用186例release、20,000个CpG和最近质心规则，两者不是同一模型。原稿关于冻结原分类器的说法必须撤回，质心结果保留为研究用替代实现。已从作者公开镜像提取真实训练参考与脚本并核验镜像层SHA256，在R中用作者相同的caret::knn3执行，另用欧氏距离手算核对概率，已检查的样本概率差异为0。

Baylor公开β值覆盖401/403探针，按训练均值填补2个探针；3例最高票并列者不强制分类，临床模型为106条记录/18事件。真实作者KNN的β值适配分支折外ΔC +0.018，表观bootstrap区间跨零；5年Brier由0.268增至0.287，说明排序的正向增量不等于绝对风险误差改善。按论文0.75阈值保留的临床完整子集为95条记录/16事件，ΔC +0.031；这是被选择的不同人群，不能用来声称拒判提高整体临床预测能力。

Tübingen的公开处理后矩阵只保留164/403探针（40.7%），未通过50%最低覆盖要求，因此没有为该队列生成作者KNN适配预测。UCSF的565个处理后甲基化样本已生成分布级结果，但未与未经证明的临床位置桥接强行关联。以上均不是作者原始IDAT预处理流程的完整复现。

### 原37例release含作者训练样本

403探针数值指纹逐维最大误差<1e-8，证明原186例release含121个作者训练样本，原37例release含20个；合并正好覆盖141个不同训练指纹。141个指纹锚点与GEO样本顺序中的barcode均一致。原37例release不能称原模型独立留出集。剔除训练重叠、无标签及并列者后，剩11例可作依赖SAMPLE编号/series顺序桥接的条件性标签比较（10/11一致）；该小规模、尚未全量ID证实的结果不进入主要独立验证结论。不能将旧稿0.774/0.806当成作者原模型独立验证性能。

### Bootstrap尾部比例不是经零假设校准的P值

原有5,000次结果和分布保留，并在 `legacy_bootstrap_reinterpreted.csv` 将尾部字段重命名。原0.034/0.037不能再解释为确认性P值。本轮以效应、区间、折外性能和乐观校正为主；正式关联检验另用未惩罚嵌套Cox模型，不能把两种证据混合。

### 惩罚模型的似然差不能直接沿用常规LR检验

已重新拟合未惩罚模型并报告两侧全局组效应，模型均成功拟合。检验为审阅后探索性分析，且不检验ΔC是否大于0。Baylor少事件，Wald及LR的渐近近似本身仍受限。

{lrtable}

UCSF的分子组别关联最稳健；其余核心比较不能因单个HR的名义P小于0.05就称整体模型已证实增量。PH诊断也全部保存：Tübingen的性别项出现名义PH偏离，不能据多重校正后不显著就宣布PH假设成立。

### 校准与净获益计算修正

旧校准脚本把标准化LP的exp(coef)当成校准斜率。已改为原始LP的未指数化Cox系数，并另外计算按折分层的折外斜率。旧IPCW漏掉超过时点非事件者的1/G(t)权重，并用不恰当的加权分母形成观察风险；这影响校准、O/P和决策曲线，旧图与旧临床效用结果不可直接沿用。

校正后的表观5年O/P：Baylor已发表分组1.105，Tübingen质心1.039，UCSF已发表分组0.962；原2.19/1.33并非可信的低估程度。相应表观校准斜率为1.208、1.133、1.149，而非原指数化数值。折外校准斜率仍会明显偏离1，尤其少事件队列，表观O/P接近1也不证明个体校准良好。

{risktable}

Brier为固定5年（UCSF以年记录，其余以月记录）的IPCW平均平方误差。ΔBrier<0表示改善；Tübingen质心分支与Baylor作者KNN全适配分支出现误差增加，必须保留。决策曲线数值已按正确权重重算，但目前阈值未绑定具体治疗决策，不据此声称临床获益。Baylor 99例比较在60月只有11例仍在随访，长时点估计精度很差，完整风险集、G(t)、KM风险及区间已保存。

### P3探针排序错误已纠正

原脚本将MAD降序名次与缺失率升序名次相加后，却按总分降序选取2000探针，选的是该评分定义下较差的一端。已改为升序，正确面板与原面板0探针重叠，并按各目标共享面板在训练数据重算P99阈值。正确面板在Baylor仅共享679探针、UCSF共享1002、994例队列共享1999，进一步表明训练内变异与缺失排名不保证外部覆盖。

正确P3：Baylor ARI 0.808、best-map 0.910、拒判3/110；UCSF有参考分区的200例ARI 0.485、best-map 0.795、拒判13/200；994例类别434/225/335、拒判31/994（3.12%）。旧P3的UCSF best-map 0.86不能继续称正确执行稳健面板后的主要性能。P1/P2未重定义；994例旧P1=439/237/318、旧P2=437/238/319，均须注明路径。分布相似不等于准确或预后有效。

### 样本身份和纵向基线

Baylor原文是110个肿瘤、108位患者；公开临床表缺患者簇标识，本轮不能完成患者簇级重采样/分组交叉验证。现有Baylor结果为肿瘤记录级探索性分析。UCSF 201条临床记录到205条记录的八字段匹配中，199条只有1候选、2条各有2候选，旧脚本字典覆盖掩盖了非唯一匹配；保存全部候选，临床行到fingerprint列的位置映射仍未证明。临床标签模型直接使用MOESM19，不依赖该位置桥接。

纵向15/17是观察到的“全class3簇数”，不是常量分类器的基线。常量输出class3将产生17/17簇内一致性，高于观察16/17；患者配对身份又未验证，故不将其作为稳定性的支持证据。未检出CpG变化也不等于无变化。

### 标签编码与终点

Tübingen旧规格是两个独立指示变量：proliferative和hypermetabolic，共同参照为NF2-wildtype与immunogenic的合并人群，**不是四水平分类变量，也不是把两个不良组合成一个二分类变量**。完整四组敏感性以NF2-wildtype为参照、设3个哑变量，折外ΔC约+0.007，须与旧规格+0.023并列。质心class3 HR不等于Toronto proliferative HR。

UCSF LFFR应展开为local freedom from recurrence，不能写late failure-free survival。时间以年记录。UCSF旧列名male实际编码I(Sex=2)，本轮按原始性别编码的二分类指示变量解释；在来源未明确2的性别含义前，不把该列HR标为男性效应。反转二分类编码不改变模型预测。

## 三 可用于稿件的阳性发现

1. 三套队列作者已发表的分组在统一的记录级内部验证中，平均折外ΔC方向均为正：+0.073、+0.023、+0.032。这是不同标签体系的证据，不能包装成同一冻结分类器的三队列验证。
2. UCSF分组在临床变量之外的增量同时见于折外C、乐观校正C及5年Brier；其关联在多重比较校正后仍有支持。
3. UCSF分组在MIB1之外仍有小幅平均预测增量，但属于事后增补分析，全部14检验Holm校正后P约0.113，应以效应而非显著性表述。
4. 质心替代模型在Baylor/Tübingen有正向平均折外ΔC，但不能声称普遍提高绝对风险预测或原模型的临床效度。

不建议保留不同分子对比的主要pooled HR；队列HR逐个报告，来源、参照与终点列明。I²=0不能消除参数不可比性。本轮未为异质对比再计算合并“阳性”。

## 四 尚需外部证据的问题

- Baylor肿瘤与患者对应关系：取得108位患者的簇ID后，重做簇级验证和重采样；当前记录级结果不能称患者独立验证。
- UCSF临床与甲基化的唯一ID桥接：不以聚类/生物学梯度或统计显著性证明身份。
- 作者KNN完整IDAT流程：已核验模型来源并完成处理后β值适配，但Tübingen处理后矩阵覆盖不足；完整一致预处理与外部性能仍待补充。
- 跨来源患者重叠：已证明本研究训练release之间的重叠；其他队列间完整患者级排重仍受公开标识限制。
- 外部校准及临床决策验证、公开代码访问安排、作者/伦理/资助/贡献声明仍需要真实材料。不能以本轮统计补强代替这些材料，也不保证期刊接收。

## 五 文件与验证

`results/validation_summary.csv` 为14规格汇总；`oof_predictions_*` 与 `fold_assignments_*` 支持逐例核查；`bootstrap_optimism_*` 保存每次重采样；`unpenalized_association_tests.csv` 为正式但探索性的全局关联检验；`oof_risk_summary.csv` 为正确删失权重下的折外误差；`corrected_P3_transfer_summary.csv` 与 `official_knn_audit_summary.csv` 记录部署及训练重叠。

67项包内一致性核查及8项独立手算统计检查通过。原稿与两份原关键结果的SHA256未改变。无新增图文件。结果表、英文可替换段落、代码与哈希清单均随包提供。本报告不将“已完成本轮计算”写成“已解决所有数据身份与独立外部验证问题”。

## 已核验来源

{refs}

期刊结构核对来源：[npj Digital Medicine content types](https://www.nature.com/npjdigitalmed/content-types)。本轮先完成数据，不据此声称投稿材料已齐备。原作者模型来源：[NWM-MGP Public Resources](https://github.com/NWM-MGP/Public_Resources/tree/develop/meningioma_risk_classifier)，镜像v0.1.0的实际层摘要保存在包内。
'''
(DEST/'数据统计补强报告_20261019.md').write_text(report,encoding='utf-8')

engrows=[]
for key in CORE:
    z=s.loc[key];engrows.append([key,f'{int(z.n)}/{int(z.events)}',f(z.C_app_clinical),f(z.C_app_molecular),f(z.cv_C_clinical),f(z.cv_C_molecular),sign(z.cv_delta_mean),sign(z.delta_corrected)])
engtable=table(['Comparison','Records/events','Apparent C clinical','Apparent C molecular','CV C clinical','CV C molecular','CV delta C','Optimism corrected delta C'],engrows)
english=f'''# Revised data results and methods for the meningioma manuscript

These sections replace the corresponding statistical claims in manuscript_npjDM.md. This is a retrospective audit revision dated 19 October 2026. Figures and submission declarations are outside this data revision.

## Results

### Cohort roles and analysis units

The Baylor series comprised 110 primary meningioma specimens from 108 patients [1], rather than 110 independent patients. The public clinical table does not establish the patient clusters, so analyses of this series were performed at the tumour-record level and remain exploratory. Complete-case comparisons included 99 records and 17 recurrence events for the published MenG labels and 109 records and 19 events for the centroid surrogate. In Tübingen, 158 primary tumours without prior radiotherapy contributed 149 complete clinical records and 55 progression events. The UCSF clinical table contained 205 labelled records, of which 201 with 76 local-recurrence events had complete model covariates. Its endpoint was local freedom from recurrence, recorded in years [3]. The 200 methylation profiles with discovery fingerprints were analysed separately from these clinical models. The 994-sample multicentre methylation release contributed distribution and quality-control information only.

### Identity audits change the interpretation of classifier validation

The existing research implementation is a 20,000-CpG nearest-centroid surrogate learned from a public methylation release. It is not the published Duckett classifier, which uses a 141-case training reference, 403 CpGs and K-nearest-neighbour classification with K=10 [2]. Using all 403 CpGs and a maximum absolute numerical difference below 1e-8, we matched 121 reference training profiles to the 186-profile release and 20 to the 37-profile release. These matches covered all 141 unique reference profiles and agreed with their barcode annotations under the series-order bridge. Thus, the 37-profile release is not an independent held-out sample set for the author's model. Its prior label agreement estimates must not be presented as independent validation of that model.

For UCSF, 199 of the 201 clinical records had a unique match on the shared eight-field key and two had two candidate source records. A dictionary-based linkage had hidden this non-uniqueness. More fundamentally, clinical-row order has not been proven to correspond to fingerprint-column order. Biological gradients cannot establish individual identity. UCSF prognostic analyses therefore used the native clinical table and its published group labels, without importing methylation predictions through the positional bridge.

### Molecular labels retain cohort-dependent prognostic increments after internal validation

With identical complete cases and paired event-stratified folds, mean repeated cross-validated discrimination increased for all five core comparisons. Published MenG labels increased C from 0.634 to 0.707 in Baylor; the published Tübingen Toronto labels under the existing two-indicator coding increased C from 0.541 to 0.564; and the UCSF published groups increased C from 0.710 to 0.741. The centroid surrogate increased C from 0.609 to 0.638 in Baylor and from 0.541 to 0.552 in Tübingen. These comparisons evaluate different molecular variables; they do not establish a single classifier's prognostic performance across all three cohorts.

{engtable}

The split-to-split 2.5th–97.5th percentiles of delta C were −0.0003 to 0.161 for Baylor published labels, −0.015 to 0.060 for Tübingen published labels, and 0.006 to 0.045 for UCSF published groups. These ranges describe variation across splits of the same records and are not population confidence intervals. Optimism-corrected delta C was +0.072, +0.045 and +0.033, respectively. Apparent delta C bootstrap intervals remained compatible with zero in Baylor and Tübingen. Unpenalized global association tests yielded raw P=0.0324, 0.0312 and 0.000689 for the respective published-label comparisons; after Holm correction across the five core comparisons, the corresponding P values were 0.125, 0.125 and 0.00345. These are exploratory association tests, not tests of out-of-sample delta C.

In UCSF, adding the published groups to clinical covariates and MIB-1 increased mean cross-validated C from 0.739 to 0.751 (delta C +0.012). The unpenalized group-effect P was 0.0113 before correction and 0.113 after Holm correction across all 14 comparisons. Exclusion of the single zero-time event preserved a mean delta C of +0.032. The full four-category Tübingen model, retaining NF2-wildtype and immunogenic as distinct categories, had a smaller mean delta C (+0.007) than the original two-indicator specification (+0.023). All specifications are retained rather than selecting the most favourable coding.

### Correct censoring weights revise the calibration findings

The previous observed/predicted ratios of 2.19 and 1.33 resulted from incorrect censoring weights and normalization. With corrected inverse-probability-of-censoring weighting, apparent five-year observed/predicted ratios were 1.105 for Baylor published labels, 1.039 for the Tübingen centroid model and 0.962 for UCSF published groups. The corresponding Cox recalibration slopes on the original linear-predictor scale were 1.208, 1.133 and 1.149. Previously reported exponentiated coefficients of standardized predictors were hazard ratios per standard deviation, not calibration slopes. These apparent values do not establish individual calibration, and cross-fitted slopes remained variable.

In cross-validation, the five-year Brier score decreased from 0.264 to 0.260 for Baylor published labels, from 0.248 to 0.242 for the Tübingen two-indicator published-label model and from 0.199 to 0.183 for UCSF published groups. By contrast, it increased from 0.248 to 0.250 for the Tübingen centroid model. Beyond clinical covariates and MIB-1 in UCSF, adding groups reduced the Brier score from 0.183 to 0.175. A positive discrimination increment therefore did not guarantee an improvement in absolute-risk error. At five years, only 11 of the 99 Baylor records remained under observation; long-horizon performance in that cohort was imprecise. Decision-curve quantities were recomputed but are not interpreted as clinical benefit without a defined decision and external evaluation.

### The author's KNN can be assessed only within the available data's preprocessing limits

The author's public 403-probe reference matrix and prediction script were retrieved from the versioned image and their layer hashes verified. We executed the same caret::knn3 operation and independently reproduced its probabilities for checked samples. Baylor processed data contained 401 of the 403 probes. After training-mean imputation of the two absent probes and leaving three tied maximum-vote samples unassigned, 106 records with 18 events were available. Cross-validated C increased from 0.590 to 0.608 (delta C +0.018), but the five-year Brier score increased from 0.268 to 0.287. The source-specified confidence threshold of 0.75 retained 95 clinically complete records with 16 events, with delta C +0.031 and Brier 0.270 versus 0.256. These differing subsets cannot demonstrate an overall benefit of rejection. The Tübingen processed release contained only 164/403 probes, failing the 50% minimum coverage rule. No adapted KNN prediction was generated for it. All these experiments use processed beta values and do not reproduce the original IDAT preprocessing or its clinical assay validation.

### Correcting probe selection weakens the previous robust-panel claim

The intended robust-panel score summed the descending MAD rank and ascending missingness rank. The original script erroneously selected the highest scores. Correct ascending selection produced a disjoint 2,000-probe panel, with 679 probes shared in Baylor, 1,002 in the full UCSF release and 1,999 in the multicentre release. Thresholds were recalculated in the training data for each shared panel. Agreement with Baylor's published labels was ARI 0.808 and best-map agreement 0.910, with 3/110 distance flags. In the 200 UCSF profiles with the existing unsupervised reference, ARI was 0.485 and best-map agreement 0.795, with 13/200 flags. Thus the previous 0.860 result cannot be retained as the correctly executed robust-panel result. The multicentre class counts were 434/225/335 with 31/994 flags. All class proportions are descriptive; reference-free distributions cannot establish accuracy or prognostic validity.

### Longitudinal evidence does not establish class stability

Sixteen of 17 computationally inferred clusters had concordant labels. However, 15/17 was the number of observed all-class-3 clusters, not the performance of a majority-class baseline. Assigning every specimen to class 3 would produce 17/17 concordant clusters. Together with unresolved patient identity and the lack of detectable probe changes after multiple-testing correction, these findings do not establish stable longitudinal classification or equivalence.

## Methods

### Analysis status and model definitions

This revision was prompted by a retrospective audit. Neither its tests nor its sensitivity analyses are described as preregistered. The five core comparisons retain the previous covariate sets, with a fixed ridge penalty of 0.05 for prediction. Baylor clinical covariates were WHO grade 2 versus 1, gross-total resection, age and sex. Tübingen used WHO grade as an ordinal variable, gross-total resection, age and sex; categorical grade was evaluated separately. UCSF used ordinal WHO grade, gross-total resection and a binary native sex-code indicator because age was unavailable. The legacy column named male in that table represents I(Sex=2); its direction is not interpreted as a male effect without source coding confirmation. Fixed ridge fitting and its internal covariate scaling used training records only. Constants in a training sample were omitted for that fit; the penalty was not increased to rescue failures.

MenG classes used A as reference and two indicators for B and C. Centroid classes used 1 as reference and indicators for 2 and 3. The existing Tübingen published-label model used separate proliferative and hypermetabolic indicators with NF2-wildtype and immunogenic combined as the shared reference. It was neither a four-level factor nor one binary unfavourable-group indicator. The four-level sensitivity used NF2-wildtype as reference. UCSF groups used Merlin-intact as reference with hypermitotic and immune-enriched indicators. An identical complete-case set was constructed from the union of covariates for both models in every comparison. Counts denote records rather than independently verified patients when patient identity was not available.

### Paired internal validation and optimism correction

Each comparison used 20 repetitions of five-fold event-stratified cross-validation, with seeds 20261019 through 20261038. Both models used the same folds. Each test record was predicted once per repeat. Concordance was calculated using comparable pairs within each test fold, and concordant and comparable-pair counts were accumulated before calculating the repeat-level C; independently fitted Cox scores were not mixed across folds. Mean paired delta C and the distribution across repeats were saved. Split percentiles were not labelled confidence intervals or converted into P values. Baylor's missing patient-cluster identifiers prevent this from being patient-grouped validation.

For five core comparisons and the added Baylor author-KNN beta adaptation, 1,000 paired record-bootstrap resamples refitted both models. Apparent percentile intervals describe the in-sample delta C estimator. Model optimism was the performance difference between the bootstrap training data and the original data evaluated with the bootstrap-fitted model. Mean paired optimism was subtracted from the original apparent delta C. Monte Carlo standard errors describe simulation precision, not uncertainty across new patient samples. No interval for corrected delta C was inferred by simply taking percentiles of these optimism draws. All 6,000 paired resamples succeeded. Descriptive fractions with delta C at or below zero were archived but were not interpreted as null-calibrated P values.

### Association tests and diagnostics

Conventional likelihood-ratio tests compared nested unpenalized Cox models on identical records, with degrees of freedom equal to the number of added coefficients. Differences in penalized likelihoods were not used for these tests. Raw two-sided global P values and Holm corrections across the five core and all 14 evaluated comparisons were reported. They assess exploratory association and do not test the prediction-performance increment. Cohort-specific hazard ratios were retained with their exact contrasts and adjustment variables; no primary pooled hazard ratio combined unlike classes or endpoints. Rank-transformed Schoenfeld-residual diagnostics were archived. Nominal deviations, particularly for sex in Tübingen, remain reasons for caution even when multiple-testing adjustments are not significant. The low event count in Baylor limits all asymptotic inference.

### Censoring, calibration and absolute-risk error

Evaluation horizons were 36 and 60 months, or 3 and 5 years for UCSF. The censoring Kaplan–Meier distribution was estimated in each training fold. Events at or before the horizon received weight 1/G(T−), records observed beyond it received 1/G(h), and censoring before or at the horizon without an event received zero. Unsupported horizons or zero censoring survival were flagged rather than assigned clipped weights. Brier scores were weighted squared errors divided by the number of test records. Observed risk was the weighted event count divided by the same denominator; predicted risk was the mean Cox event probability. Weighted net benefit used the same event and control weights. Fold-level horizon and decision quantities were combined by test-record counts within each complete repeat. Risk sets, censoring survival and apparent Kaplan–Meier risks with intervals were saved.

For Cox recalibration, the response was regressed on the original, unstandardized linear predictor, and the coefficient itself was reported as slope. Out-of-fold slopes were fitted separately for each repeat with fold-specific baseline hazards. Neither exponentiation of the slope nor standardization before comparing it to 1 was used. Aggregate O/P close to 1 was not treated as proof of individual calibration. Complete source and numerical metadata were saved with the execution environment. All results are conditional on the source endpoint semantics and on records retained in the public releases.

## Discussion replacement

Molecular labels carried prognostic information beyond clinical variables, with the clearest reproducible increment in the UCSF clinical release. The evidence supports cohort-dependent risk stratification more strongly than universal classifier transportability. In particular, internal validation of clinical coefficients within each cohort does not evaluate a single frozen clinical risk equation across cohorts. The centroid implementation is a research surrogate, and neither its label agreement nor the beta-level adaptation of the author's KNN inherits the clinical validation of the original assay. Probe coverage, preprocessing, identity linkage and training overlap must be resolved before describing independent external classifier validation. Small positive mean delta C values and positive associations should be interpreted alongside split variability, absolute-risk error, competing molecular definitions and incomplete patient identifiers. Future evaluation requires harmonized IDAT processing, confirmed individual links, patient-grouped validation where applicable and prospectively defined clinical decisions.

## Verified references

{refs}
'''
(DEST/'Results_Methods_可替换英文段落.md').write_text(english,encoding='utf-8')

readme='''# 数据补强包说明

先读 `数据统计补强报告_20261019.md`。`Results_Methods_可替换英文段落.md` 可用于下一轮稿件修改；原主稿保留，不能与旧图/旧统计解释直接混用。本轮没有制作新图。

## 主结果

- `results/validation_summary.csv`：14组模型的表观、重复交叉验证与乐观校正结果。
- `results/unpenalized_association_tests.csv`：探索性全局关联检验和两套Holm校正。
- `results/oof_risk_summary.csv`：正确删失权重下的折外O/P与Brier，含两个固定时点。
- `results/cohort_specific_unpenalized_HR.csv`：队列/变量特异HR，不做异质主要合并。
- `results/official_knn_audit_summary.csv`：真实作者模型、训练重叠、覆盖与并列分类。
- `results/corrected_P3_transfer_summary.csv`：纠正排序后的P3。
- `results/ucsf_clinical_key_cardinality.csv`：逐条临床键匹配候选数，不隐藏重复。
- `results/package_validation.json`、`statistical_checks.json`：67项一致性检查与8项手算统计检查。

## 复算核心预后结果

Python环境版本见 `results/environment.json`，安装清单见 `requirements_exact.txt`。分析脚本支持从随包保存的14个完整病例输入表复算，不必重新下载完整甲基化矩阵：

```powershell
rtk proxy python run_strengthening.py --all --snapshot
rtk proxy python summarize_diagnostics.py
rtk proxy python verify_statistics.py
```

运行会覆盖本文件夹中的派生结果；要比较字节/数值一致性，请先复制完整包到新文件夹。不同版本线性代数库可能导致最后小数位差异。原始记录级输入CSV只是各来源数据的整理分析集；Baylor缺少患者簇ID，不应称为已验证的独立患者数据。

作者KNN模型可从本包的 `results/official_model/knn_reference_fc2.rds`、403探针矩阵和 `predict_official_knn.R` 复算，需R和caret/class。来源脚本、镜像摘要及R会话信息已保留；作者脚本未直接执行，其中绘图部分不属于本轮工作。运行预测后再运行 `audit_official_predictions.py` 才能恢复训练重叠标记，不能把初始的名称匹配标记当作已排除重叠。

P3面板排序重建和源表审计需要原项目的 `data/` 与 `results/` 大矩阵/缓存；完整来源与探针排名均保存。本包并未重新打包数十GB原始下载数据。冻结学习参考、处理后β值适配、原始IDAT预处理、临床模型系数拟合是不同步骤。

## 重要解释

本轮为事后审计补强。折分百分位不是总体置信区间；bootstrap尾部比例不是校准P值；关联检验不等于预测增量检验；模型不同队列内部验证不等于冻结临床风险模型的外部验证。作者真实KNN、旧质心替代模型及正确P3分开报告。仍需解决的身份、患者聚类与原IDAT流程问题在报告中逐项列明。
'''
(DEST/'README.md').write_text(readme,encoding='utf-8')
env=json.loads((OUT/'environment.json').read_text())
req=[]
for module,pkg in [('numpy','numpy'),('pandas','pandas'),('lifelines','lifelines'),('sklearn','scikit-learn'),('scipy','scipy'),('statsmodels','statsmodels')]:req.append(pkg+'=='+env[module])
import openpyxl,requests
req+=['openpyxl=='+openpyxl.__version__,'requests=='+requests.__version__]
(DEST/'requirements_exact.txt').write_text('\n'.join(req)+'\n',encoding='utf-8')
print('Reports created; all source manuscripts preserved.')
