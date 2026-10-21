from pathlib import Path
import json,hashlib,platform,sys
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;OUT=ROOT/'results'
EXPECTED={'author186':3,'baylor':110,'tubingen':149,'ucsf':200,'hku':365,'longitudinal':42}
checks=[]
def check(name,ok,detail=''):
    checks.append({'check':name,'passed':bool(ok),'detail':str(detail)})
def table(d):return d.to_markdown(index=False,floatfmt='.3f')
def write(name,s):
    (ROOT/name).write_text(s.strip()+'\n',encoding='utf-8')

completion=pd.read_csv(OUT/'raw_cohort_completion_status.csv')
for cohort,n in EXPECTED.items():
    c=completion[completion.cohort==cohort];check(cohort+'_requested_count',len(c)==n,len(c))
    for flag in ['raw_pair_complete','author_preprocessed','raw_qc']:check(cohort+'_'+flag,c[flag].all(),int(c[flag].sum()))
if not all(x['passed'] for x in checks):
    pd.DataFrame(checks).to_csv(OUT/'identity_validation_checks.csv',index=False);raise SystemExit('Incomplete cohorts: final report withheld')
summ=pd.concat([pd.read_csv(OUT/'raw_preprocessed'/(g+'_summary.csv')) for g in completion.GSM],ignore_index=True)
qc=pd.concat([pd.read_csv(OUT/'raw_qc'/(g+'_qc.csv')) for g in completion.GSM],ignore_index=True)
summ.to_csv(OUT/'complete_author_raw_preprocessing_summary.csv',index=False);qc.to_csv(OUT/'complete_raw_quality_summary.csv',index=False)
check('classifier_403_features_every_profile',summ.finite_classifier_cpgs.eq(403).all())
check('no_classifier_missing_values_or_imputation',summ.missing_classifier_cpgs.eq(0).all())
check('raw_missingness_screen_at_author_threshold',~qc.author_failed_color_rule.any(),int(qc.author_failed_color_rule.sum()))
check('probabilities_sum_one',np.allclose(summ[['p_class1','p_class2','p_class3']].sum(axis=1),1,atol=1e-12))
check('all_ties_unassigned',summ.loc[summ.tied_top,'predicted_class'].isna().all())
check('assigned_class_equals_unique_vote_max',summ.loc[~summ.tied_top,'predicted_class'].astype(int).eq(summ.loc[~summ.tied_top,['p_class1','p_class2','p_class3']].to_numpy().argmax(axis=1)+1).all())
check('batched_single_equivalence',summ.batch_single_max_abs_difference.dropna().le(1e-12).all(),summ.batch_single_max_abs_difference.notna().sum())
anchor=pd.read_csv(OUT/'reference_anchor_raw_vs_original_403beta.csv');check('three_class_reference_anchors',set(anchor.reference_class)=={1,2,3})
check('1209_anchor_betas_exact',anchor.max_abs_beta_difference.eq(0).all(),int(anchor.n_cpgs.sum()))
rawaudit=pd.read_csv(OUT/'raw_internal_barcode_and_hash_audit.csv')
check('all_internal_slide_barcodes_agree',rawaudit.internal_barcode_matches.all(),len(rawaudit))
for r in completion.itertuples():
    hit=rawaudit[rawaudit.GSM.eq(r.GSM)];check('raw_channels_'+r.GSM,len(hit)==2 and set(hit.color)=={'Grn','Red'})
    beta=pd.read_csv(OUT/'raw_preprocessed'/(r.GSM+'_403beta.csv'))
    check('exported_403beta_'+r.GSM,len(beta)==403 and beta.cpg.is_unique and np.isfinite(beta.beta).all() and beta.beta.between(0,1).all())
old=pd.read_csv(BASE/'strengthening_20261019/results/original_hashes.csv')
for r in old.itertuples():check('original_preserved_'+r.path,hashlib.sha256((BASE/r.path).read_bytes()).hexdigest()==r.sha256)
meta=pd.read_csv(OUT/'SNP_donor_clusters_current.csv')
donor_rows=[]
for cohort,n in EXPECTED.items():
    c=completion[completion.cohort==cohort].merge(meta[['GSM','donor_primary_key','donor_strict_key','donor_overlap_guard_key','callable_SNPs']],on='GSM',validate='one_to_one')
    check(cohort+'_complete_SNP_fingerprints',len(c)==n and c.callable_SNPs.notna().all(),c.callable_SNPs.min())
    donor_rows.append({'cohort':cohort,'profiles':n,'primary_SNP_clusters':c.donor_primary_key.nunique(),'strict_SNP_clusters':c.donor_strict_key.nunique(),'overlap_guard_clusters':c.donor_overlap_guard_key.nunique()})
donors=pd.DataFrame(donor_rows);donors.to_csv(OUT/'requested_cohort_SNP_cluster_counts.csv',index=False)
models=pd.concat([pd.read_csv(f) for f in (OUT/'donor_aware').glob('summary_*.csv')],ignore_index=True)
models.to_csv(OUT/'donor_aware_validation_summary_all.csv',index=False)
names=[c+'_'+m+'_overlap_guard' for c,ms in [('baylor',['published','centroid_surrogate','author_KNN_raw']),('tubingen',['published','centroid_surrogate','author_KNN_raw']),('hku',['published','author_KNN_raw'])] for m in ms]
core=models[models.comparison.isin(names)].copy();check('eight_complete_core_models',set(core.comparison)==set(names))
for r in core.itertuples():
    check(r.comparison+'_1000_cluster_bootstraps',r.bootstrap_success==1000 and r.bootstrap_attempted==1000)
    check(r.comparison+'_no_failed_fits',r.n_failures==0)
for file in (OUT/'donor_aware').glob('oof_*.csv'):
    z=pd.read_csv(file);key=file.stem[4:];d=pd.read_csv(OUT/'donor_aware'/('input_'+key+'.csv'))
    check(key+'_20_complete_repeats',z.groupby('repeat').size().eq(len(d)).all() and z['repeat'].nunique()==20)
    check(key+'_donor_fold_isolation',z.groupby(['repeat','donor']).fold.nunique().eq(1).all())
    check(key+'_unique_record_per_repeat',not z.duplicated(['repeat','id']).any())
for cohort in ['baylor','tubingen','hku']:
    common=[f for f in (OUT/'donor_aware').glob('input_'+cohort+'_*_common_set.csv')]
    if not common:check(cohort+'_paired_common_set_exists',False);continue
    populations=[set(pd.read_csv(f).GSM) for f in common]
    check(cohort+'_same_common_set_across_layers',all(x==populations[0] for x in populations))
    ass=[pd.read_csv(OUT/'donor_aware'/('oof_'+f.stem[6:]+'.csv'))[['repeat','id','fold']].sort_values(['repeat','id']).reset_index(drop=True) for f in common]
    check(cohort+'_same_folds_across_common_layers',all(x.equals(ass[0]) for x in ass))
u=pd.read_csv(OUT/'ucsf_ID_proven_molecular_linkage_200.csv');check('UCSF_200_literal_M_ID_groups',len(u)==200 and u.MID.is_unique and u.GSM.is_unique)
reject=pd.read_csv(OUT/'ucsf_old_positional_mapping_rejected.csv')
check('old_positional_bridge_rejected',reject.published_group_agrees.mean()<.5)
hku=pd.read_csv(OUT/'hku_ID_linked_clinical_365.csv');check('HKU_365_case_number_join',len(hku)==365 and hku.GSM.is_unique and hku.ID.is_unique)
check('HKU_primary_excludes_demographic_conflicts',hku.loc[hku.identity_demographic_conflict,'identity_use_primary'].eq(False).all())
mp=pd.read_csv(OUT/'longitudinal_verified_patient_mapping_42.csv');pc=pd.read_csv(OUT/'longitudinal_centroid_identity_candidates.csv')
check('longitudinal_42_samples_18_groups',len(mp)==42 and mp.GSM.is_unique and mp.source_patient_centroid.nunique()==18)
check('longitudinal_12_pairs_6_triplets',mp.groupby('source_patient_centroid').size().value_counts().to_dict()=={2:12,3:6})
check('longitudinal_unique_numeric_decompositions',pc[pc['rank'].eq(1)].exact_numeric_candidate.all() and not pc[pc['rank'].ne(1)].exact_numeric_candidate.any())
long=pd.read_csv(OUT/'longitudinal_raw_patient_stability_summary.csv')
overlap=pd.read_csv(OUT/'raw_author_reference_overlap_audit.csv');check('no_external_literal_training_barcode_overlap',not overlap.loc[~overlap.cohort.eq('author186'),'barcode_in_author_reference'].any())
check('no_external_exact_403_profile_training_copy',not overlap.loc[~overlap.cohort.eq('author186'),'exact_beta_reference_match_tolerance_1e12'].any())
pd.DataFrame(checks).to_csv(OUT/'identity_validation_checks.csv',index=False)
if not all(x['passed'] for x in checks):raise SystemExit('Validation failed: '+str([x['check'] for x in checks if not x['passed']]))
label=pd.read_csv(OUT/'raw_author_supervised_label_metrics_with_CI.csv');hka=json.loads((OUT/'hku_clinical_linkage_audit.json').read_text(encoding='utf-8'))
qcsummary=qc.groupby('cohort',as_index=False).agg(n=('GSM','size'),median_detection_failure=('detection_fail_fraction','median'),maximum_detection_failure=('detection_fail_fraction','max'),raw_missingness_screen_failures=('author_failed_color_rule','sum'))
qcsummary.to_csv(OUT/'raw_quality_by_cohort.csv',index=False)
pre=summ.groupby('cohort',as_index=False).agg(profiles=('GSM','size'),all_403=('finite_classifier_cpgs','min'),missing_CpGs=('missing_classifier_cpgs','sum'),top_vote_ties=('tied_top','sum'))
columns=['comparison','records','events','donor_clusters','C_cv_clinical','C_cv_molecular','delta_cv','positive_repeats','delta_corrected','delta_app_CI_low','delta_app_CI_high']
core=core.sort_values('comparison');core.to_csv(OUT/'core_donor_aware_validation_summary.csv',index=False)
status=[{'issue':'UCSF M-ID to GSM/raw barcode/published group','status':'resolved_200_literal_IDs'},
 {'issue':'UCSF raw preprocessing','status':'validated_200_raw_profiles_403_CpGs'},
 {'issue':'UCSF per-profile clinical-outcome linkage','status':'unresolved_author_crosswalk_required; positional_mapping_rejected_and_quarantined'},
 {'issue':'HKU clinical linkage','status':'365_source_case_number_links; 358_demographically_concordant_primary; raw_SNP_overlap_guard'},
 {'issue':'Baylor and Tubingen repeated donors','status':'SNP_supported_clusters_used_for_grouped_validation; not_author_confirmed_person_IDs'},
 {'issue':'Longitudinal 18 source-labelled groups','status':'unique_numeric_reconstruction; independent_donor_count_unresolved_due_to_cross_group_SNP_matches; no_visit_order_inference'},
 {'issue':'Longitudinal patient18 group-grade annotations','status':'unresolved_source_discrepancy; no_chronological_changes_claimed'},
 {'issue':'Training versus heldout release','status':'reference_overlap_in_both_releases; not_independent_holdout'},
 {'issue':'Original preprocessing','status':'3_reference_classes_1209_betas_exact; all_869_requested_profiles_no_imputation'}]
(OUT/'resolution_ledger.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
qcn=len(qc);rawbytes=rawaudit.bytes.sum()/1e9;passed=len(checks)
write('数据身份与原始预处理验证报告_20261019.md',f'''
# 患者身份对应与原始预处理验证

本轮已完成 {len(completion)} 份指定样本的双通道原始数据预处理和质量检查，未制作图。原稿和上一轮结果保留。此报告及新结果取代旧稿中依赖位置桥接、错误患者簇或替代分类器的解释。

## 已解决与未解决

- UCSF：200 个 M-ID 与 GSM、Sentrix 条码、发表分子分组逐例对应已解决；原始 IDAT 已按作者流程处理。**临床结局与这 200 个 M-ID 的对应仍未得到直接来源证明**。旧行序桥接已被性别、分级、分组独立核对否定，停止用于结局、MIB-1 及生物学梯度。
- Baylor、Tübingen：已用原始芯片 SNP 指纹建立供者重叠防护簇，并完成供者分组交叉验证和簇 bootstrap。这里的簇是指纹支持的审计结果，不能替代作者确认的患者主索引。
- 纵向系列：42 份样本唯一恢复为原始来源中的 18 个标注组（12 对、6 个三联），并以原始 SNP 独立核查；不再使用旧的 17 个相关聚类。**源组 4/17、5/15、10/16 之间存在同供者 SNP 候选，不能确认 18 名独立患者**。第 18 组另有来源中的分级/分组冲突，不猜测随访顺序。
- 香港队列：从后续论文的公开数据取回 365 例带 QM 编号的临床结局，以相同机构病例号与 GEO 连接。7 例年龄或性别冲突者排除主分析，得到 358 份身份协变量相符记录，再作唯一分类及完整病例筛选。它是 GSE183647 的 HKU 子集，不能与原匿名临床表重复计入独立证据。

## 原始预处理的实证验证

从作者公开镜像提取 minfi 1.43.1、EPICv2 manifest 及原始 MNP 函数，在本地 R 4.3.3 执行。采用控制探针归一化 ref=10000、Illumina 背景校正、β=M/(M+U+100)，冻结 141 例、403 探针、K=10 的作者 KNN。执行环境并非声称完整 Docker 容器复现；软件来源、会话信息与提取层哈希均已保存。

三个类别各选一个参考样本。**1,209 个 β 值与作者参考矩阵逐值完全相同，最大绝对差为 0**。既有单例与批量处理等价检查也全部在 1e-12 容差内。目标数据全部使用真实 403 探针，不补训练均值，不用目标标签重训。EPICv2 原始数据恢复了此前处理矩阵缺失的大量探针。

{table(pre)}

双通道下载完成时验证 GZIP CRC 和 IDAT 头；最终复核内部芯片条码和 SHA256。核查文件共 {len(rawaudit)} 个通道（含一份额外可行性试验样本），约 {rawbytes:.2f} GB 未压缩数据。原始大文件留在本地，发布包含来源清单、哈希及小型派生结果。

{table(qcsummary)}

检测失败比例定义为 minfi detection P>0.01 的全芯片探针比例，仅作描述，不按结局或显著性增删病例。这项全芯片 QC 计算于原始信号；所有指定样本在该阶段均未达到与作者相同的“大于等于一半甲基化信号非有限值”缺失筛查阈值。作者流程处理后的 403 个分类探针已另行确认全部有限，不将原始阶段 QC 冒称为处理后全芯片逐值复核。

## 患者对应审计

UCSF 源补充表 3/4/5 明确给出 200 个 M-ID 及一致分组：Merlin-intact 72、Immune-enriched 65、Hypermitotic 63。旧行序桥接性别一致率 {reject.sex_agrees.mean():.1%}，分级 {reject.grade_agrees.mean():.1%}，分组 {reject.published_group_agrees.mean():.1%}，证据反对该桥接。年龄/分级/分组的协变量匹配仅保留为候选，不提升为身份确证。旧“UCSF/Toronto”名称修正为 UCSF/HKU；数据来源不包括 Toronto。

香港表来自 [Nguyen 等 Nat Commun 2025](https://doi.org/10.1038/s41467-025-61063-y) 配套 [Zenodo 11501567](https://zenodo.org/records/11501567)。QMnn 与 GEO HKU Meningioma nn 是机构病例号前缀规范化，逐项保留审计。365 例性别中 {hka['independent_sex_matches']} 一致；与 2024 融合转录本论文重叠的 302 例，LF、LFFR、死亡状态及 OS 均完全对应。301 例年龄在 0.05 岁舍入误差内，QM84 年龄显著冲突；连同性别冲突共 7 例排除主分析。旧版与 WHO2021 分级变化单列，不用于重新分配身份；主模型统一用 GEO 原分级。

芯片 SNP 预先固定阈值：β≤0.20 / 0.35–0.65 / β≥0.80 三档，其余留缺失。主要同供者候选为共同高置信位点≥40、不一致≤1且比例≤0.02；严格阈值要求 0 个不一致；重叠防护阈值为共同位点≥35、不一致≤2且比例≤0.05。所有阈值在完整重算前固定，不按结局调参。

{table(donors)}

簇数量不等同作者确认的患者数。原始记录全部留存，主要验证在簇层面分配训练/测试和重抽样；另用严格阈值及每簇按最小 GSM 固定保留一条记录作敏感性分析。后者是去重敏感性，不能解释为已确定首次手术。

## 监督分组验证

作者输出是低/中/高风险类别。它与 Choudhury 生物学三分组并非完全同义；训练参考中的多数对应一致率为 127/141（90.1%）。跨命名映射仅由这 141 例参考的条码与 GEO 元数据固定，未在 UCSF/HKU/Tübingen 标签上重新拟合。以下指标衡量该固定映射下与发表分子分组的监督一致性，不能称真实复发风险类别的诊断准确率。Tübingen 使用作者逐例公布的 UCSF 三分组；临床模型中的 Toronto 两个不利类别是另一项分析，不能混用。

{table(label[['cohort','subset','total_profiles','assigned_profiles','coverage_all_profiles','accuracy_fixed_training_semantics','accuracy_CI_low','accuracy_CI_high','balanced_accuracy_fixed_training_semantics','ARI']])}

最高票并列不强制赋类；0.75 为作者已有阈值的次级描述。高置信度子集的更高一致率必须与覆盖率同时报告，不能外推到所有病例。表中一致率区间来自 1,000 次供者簇 bootstrap，是条件于已获得分类/已达到阈值记录的区间；不包括来源身份本身的未确定性。目标标签最优映射仅保留为描述性上界。

## 新临床重算

每个比较的临床和临床+分子模型使用相同记录。固定 ridge=0.05；20 次 5 折供者分组交叉验证，C-index 只在同一测试折内计可比对，之后累加。每个主要比较进行 1,000 次供者簇 bootstrap，同时重拟合两模型并报告乐观偏倚校正点估计。预测分子层被冻结，Cox 临床风险模型则在各队列内部重新拟合，因此本分析不是一个完整冻结临床风险模型的跨队列验证。

{table(core[columns])}

delta_cv 为重复划分平均增量；positive_repeats 的分母为 20。CI 列是供者簇 bootstrap 的**表观 ΔC**百分位区间，不能称交叉验证 ΔC 区间。重复划分的百分位范围仅反映划分波动，不能称总体置信区间。bootstrap 零以下比例不再当作 P 值，也不沿用惩罚似然比检验或旧 pooled HR。唯一分类失败和临床缺失会改变不同层的分析分母。

为公平比较不同分子层，还在相同完整病例集合、相同供者分折中重算，结果见 `donor_aware/*common_set*`，不可仅用不同分母的核心表数值排名。每簇一条记录、严格 SNP 阈值等敏感性见全部结果汇总。

## 纵向结论的边界

源标注组来自 PC1、PC2 和甲基化签名均值的唯一分解，误差界与候选分离度均留存。独立 SNP 核查及作者原始分类器结果见 `longitudinal_source_groups_independent_SNP_validation.csv`、`longitudinal_raw_patient_stability_summary.csv`。{int(long.all_pairs_SNP_relaxed_agree.sum())}/18 个组的全部配对满足重叠防护阈值；源组 4/17、5/15、10/16 之间存在同供者候选，原始 SNP 图给出 16 个重叠防护簇。该数量不是确认患者数；组内不一致也可能受肿瘤 LOH 等影响，不能单凭阈值拆为不同患者。**记录分组恢复已解决，独立患者数仍需来源核对。**原始作者分类器在 17 个全获得唯一分类的源组中 14 个组内一致，另有 1 个源组含并列最高票。第 18 组的来源注释冲突未修补。所有随访统一赋为多数类别的平凡策略即可达到 18/18 组内一致，因此较高类别稳定率本身不证明临床有效性。旧错误患者簇上的配对 CpG 差异分析不再引用。

## 验证与投稿状态

通过 {passed} 项核验，包括原始条码、参考矩阵逐值复现、探针完整性、供者分折隔离、共同分析集一致性、原文件哈希等。详见 `results/identity_validation_checks.csv`。

本轮实质解决了原始预处理和可公开核查的样本/分组连接，提供了可用的 HKU 临床验证路径。**UCSF 200 例的逐样本临床结局对应仍需作者提供对照表；完成之前，相关分类器结局和 MIB-1 梯度不得恢复。** 原稿也不能直接按旧数字投稿。尚需全稿一致性替换、正式参考文献/作者声明和目标期刊要求核对；本包不声称已达到最终投稿状态。
''')

tb=core.set_index('comparison').loc['tubingen_author_KNN_raw_overlap_guard'];bb=core.set_index('comparison').loc['baylor_author_KNN_raw_overlap_guard'];hk=core.set_index('comparison').loc['hku_author_KNN_raw_overlap_guard']
write('Methods_Results_identity_repair_EN.md',f'''
# Replacement Methods and Results: identity and raw preprocessing

These paragraphs supersede positional clinical linkage and correlation-derived patient clusters. This repair was prompted by a retrospective audit; it was not prospectively registered.

## Methods

Sample identity was traced through cohort-specific case identifiers, GEO sample accessions and Sentrix barcodes. UCSF molecular groups were linked using literal M-identifiers in Supplementary Tables 3–5 of Choudhury et al. (Nat Genet 2022). Clinical rows without those identifiers were not linked by order. HKU clinical outcomes were obtained from the public MeningiomaClinical.xlsx release (Zenodo 11501567) accompanying Nguyen et al. (Nat Commun 2025). QM case numbers were normalized to the case numbers in GEO HKU sample titles. Sex was checked across all 365 records, and age, grade and outcomes were checked in 302 overlapping cases from Zakimi et al. (Acta Neuropathol 2024). Seven records with sex or age discrepancies were excluded from the primary HKU analysis. Grade-version differences were retained as annotations, and the primary clinical model used the original GEO grade consistently.

We retrieved raw paired-channel IDAT files for 110 Baylor profiles, 149 complete-covariate Tübingen profiles, 200 UCSF discovery profiles, 365 HKU profiles, 42 longitudinal profiles, and three author-reference anchors. GZIP integrity, IDAT headers, internal slide barcodes and file hashes were checked. The author image supplied minfi 1.43.1, the EPICv2 manifest and the original MNP preprocessing functions, which were executed in native R 4.3.3. Control normalization used a fixed reference of 10,000, followed by Illumina background correction and beta=M/(M+U+100). The original 141-profile, 403-CpG KNN reference and K=10 were frozen. No target-cohort training or CpG imputation was used. Tied maximum class votes were left unassigned; a confidence threshold of 0.75 was evaluated descriptively using the author's threshold. The cross-nomenclature mapping between risk categories and biological molecular groups was fixed using training-reference metadata alone (127/141 agreement). These nomenclatures are related but not identical; agreement with a biological group is not diagnostic accuracy for a true recurrence-risk category.

Raw SNP beta values were discretized at fixed high-confidence intervals (≤0.20, 0.35–0.65 and ≥0.80), leaving intermediate values missing. Primary same-donor candidates required at least 40 jointly callable SNPs, at most one mismatch and mismatch fraction ≤0.02. A conservative overlap guard used at least 35 SNPs, at most two mismatches and fraction ≤0.05; zero-mismatch clustering was a sensitivity analysis. These genotype-supported clusters were used to prevent donor overlap between training and validation folds. They are not author-confirmed unique patient identifiers. Lesion records were retained, with a one-record-per-cluster sensitivity selected by fixed GSM order, without using outcomes to select a record.

Cox models used ridge penalization of 0.05, with clinical and clinical-plus-molecular models fit on the same records. We performed 20 repetitions of five-fold stratified group cross-validation and 1,000 donor-cluster bootstrap resamples per core comparison. Comparable pairs were counted within each test fold, avoiding direct comparison of risk scores from different fitted models. Optimism-corrected C-index differences are point estimates. Bootstrap percentile intervals refer to apparent C-index differences; repeated-split quantiles describe split variability rather than population confidence intervals. No bootstrap tail fraction is presented as a calibrated P value. The molecular classifier was frozen, whereas clinical risk coefficients were refitted within each cohort.

## Results

All {len(completion)} requested profiles yielded 403 finite classifier CpGs without imputation. Across one raw reference anchor from each class, all 1,209 beta values exactly matched the author reference (maximum absolute difference, 0). Batched-versus-single preprocessing checks were within 1e-12. Raw EPICv2 data supplied all classifier CpGs despite extensive missingness in the previously used processed release. Full-array raw-signal QC found no profile meeting the author's 50% missing-signal threshold; this raw-stage screen is distinct from post-preprocessing verification of all 403 classifier probes.

UCSF literal identifiers linked all 200 profiles to published groups (72 Merlin-intact, 65 immune-enriched and 63 hypermitotic). The former positional clinical bridge achieved only {reject.sex_agrees.mean():.1%} sex, {reject.grade_agrees.mean():.1%} grade and {reject.published_group_agrees.mean():.1%} molecular-group agreement and was rejected. No outcome or MIB-1 result relying on that bridge was retained. The recovered HKU source provided 365 numbered clinical records, with 358 demographically concordant records before classification and complete-case filtering.

The requested-cohort SNP cluster counts were:

{table(donors)}

In the raw author-KNN analysis, donor-grouped mean C-indices were {bb.C_cv_clinical:.3f} versus {bb.C_cv_molecular:.3f} in Baylor ({int(bb.records)} lesion records, {int(bb.donor_clusters)} donor clusters; delta={bb.delta_cv:.3f}), {tb.C_cv_clinical:.3f} versus {tb.C_cv_molecular:.3f} in Tübingen ({int(tb.records)} records, {int(tb.donor_clusters)} clusters; delta={tb.delta_cv:.3f}), and {hk.C_cv_clinical:.3f} versus {hk.C_cv_molecular:.3f} in HKU ({int(hk.records)} records, {int(hk.donor_clusters)} clusters; delta={hk.delta_cv:.3f}). These estimates quantify within-cohort incremental discrimination and do not establish transportable absolute-risk prediction. Complete paired populations, uncertainty estimates and common-set comparisons are supplied separately.

All 42 longitudinal profiles were assigned to 18 source-labelled centroid groups by unique numerical decomposition of the published PC1, PC2 and methylation-signature means, replacing the incorrect 17-group correlation clustering. Raw SNP fingerprints supported all within-group pairs in 16 of 18 source groups under the overlap-guard threshold, but also linked source groups 4/17, 5/15 and 10/16 across their published group boundaries. The resulting 16 SNP graph components cannot be equated with confirmed patients. Thus source-record grouping was recovered, whereas the number of independent patients remains unresolved. The raw author classifier was concordant in 14 of 17 fully assigned source groups, and one group contained a tied class vote. A further source inconsistency in the grade/group annotations of patient 18 remains unresolved, so visit order and paired temporal methylation changes were not inferred. The trivial majority-class assignment yields perfect within-group concordance; class stability alone therefore cannot establish longitudinal utility.

## Remaining limitation

The publicly retrieved sources still do not establish the clinical-outcome crosswalk for the 200 UCSF M-identifiers. Those classifier-outcome analyses remain quarantined pending an author-provided mapping. Literal-barcode and exact classifier-profile checks found no copies of the author reference among external profiles, but cannot exclude different tumors from the same person in the reference without a complete donor crosswalk.
''')

write('未解决来源核对清单与作者询问草稿.md','''
# 仍需来源作者核对的两项

本轮未向任何作者发送邮件。以下是可直接核对的具体缺项，不影响已经完成的原始预处理、分组验证及 HKU 路径。

1. UCSF M-ID/GSM/Sentrix 与原 Nat Genet 2022 临床源表之间的逐例对照；同时提供去标识患者编号以区分重复肿瘤/手术、时间单位及事件编码。不能用临床表行序、年龄近邻或模型表现反推身份。
2. Tübingen 纵向来源中，源组 4/17、5/15、10/16 之间出现高度一致 SNP 指纹，同一组内也有两组未达到重叠防护阈值，需核对去标识患者主索引。Patient 18 的 Supplementary Data 4 与 Fig. 5f 源数据分级/分组组合不一致，还需提供 T26_0299、T26_0300 对应的真实随访顺序及正确注释。

## UCSF 作者询问草稿

We are auditing a secondary analysis of the public meningioma methylation release GSE183647. We can trace the 200 UCSF M-identifiers to GEO accessions, Sentrix barcodes and the published molecular groups, but cannot verify a case-level crosswalk to the clinical source tables. Would you please provide a deidentified crosswalk containing M-identifier, GEO accession or Sentrix barcode, the corresponding clinical-source row/record identifier, and a patient identifier linking repeated tumors or surgeries? We would also appreciate confirmation of the event coding and time units for local freedom from recurrence and overall survival. We have not assumed that clinical row order matches methylation-column order, and will withhold the affected classifier-outcome analyses until this linkage is verified.

## Tübingen 作者询问草稿

We recovered the 42 longitudinal profiles into the 18 published source-centroid groups. Raw-array SNP fingerprints indicate potential donor overlap across source groups 4/17, 5/15 and 10/16, while two source groups have internal pairs above our conservative SNP-discordance threshold. Could you provide a deidentified patient-to-T26 crosswalk to resolve these discrepancies and confirm repeated-donor profiles in the primary cohort? For patient 18, Supplementary Data 4 and the source data for Fig. 5f also appear to assign different grade/molecular-group combinations to T26_0299 and T26_0300; please confirm their specimen order and annotations.
''')
write('README.md',f'''
# Identity and raw-preprocessing validation, 2026-10-19

Read `数据身份与原始预处理验证报告_20261019.md` first. Results do not include figures. Original manuscript and previous releases were preserved. UCSF clinical-outcome linkage remains unresolved; its old positional bridge is rejected.

Main outputs:
- `results/complete_author_raw_preprocessing_summary.csv`: all {len(completion)} requested profiles.
- `results/raw_internal_barcode_and_hash_audit.csv`: channel-level provenance and SHA256.
- `results/ucsf_ID_proven_molecular_linkage_200.csv` and `hku_ID_linked_clinical_365.csv`.
- `results/SNP_pairwise_identity_metrics.csv` and `SNP_donor_clusters_current.csv`.
- `results/core_donor_aware_validation_summary.csv`; complete repeated-CV, fold predictions and 1,000 donor-bootstrap outputs in `results/donor_aware/`.
- `results/model_snapshot/`: frozen clinical/model inputs for portable numerical reproduction.
- `results/identity_validation_checks.csv`: {passed} successful integrity/statistical checks.

Reproduction:
1. Use Python with pandas, numpy, scipy, scikit-learn and lifelines, as recorded in `results/analysis_environment.json`.
2. `analyse_donor_aware_models.py` reads the frozen `results/model_snapshot/` directly and uses `donor_model_utilities.py`; it does not need the older workspace for model refits. Core examples: `--cohort baylor --layer author_KNN_raw`, `--cohort tubingen --layer author_KNN_raw`, `--cohort hku --layer author_KNN_raw`.
3. Optional flags: `--common-complete-set`, `--group-rule strict`, `--one-record-per-donor`. Split quantiles are not confidence intervals.
4. Raw reproduction uses the public URLs in `results/raw_idat_manifest.csv` and requested-cohort lists. Large IDAT files are kept locally and omitted from the archive. Verify source hashes, then use R 4.3.3 and the image-extracted packages in `R_library/` with the recorded manifests/dependencies. The original author reference and preprocessing source are included under `author_model/` in the release archive. The live workspace also retains their original location under `strengthening_20261019/results/official_model/`.
5. `run_author_batch.R` executes the original functions and fixed KNN. `run_raw_qc.R` and its partitioned equivalent provide descriptive full-array QC. No QC threshold was selected using outcomes. The portable raw script in the release resolves `author_model/` before the original workspace path.

Data provenance:
- [GSE183647, UCSF/HKU](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183647).
- [HKU numbered clinical data, Zenodo 11501567](https://zenodo.org/records/11501567), linked by the [2025 original paper](https://doi.org/10.1038/s41467-025-61063-y).
- [Author classifier source](https://github.com/NWM-MGP/Public_Resources/tree/develop/meningioma_risk_classifier) and quay.io/bioinformaticsnm/public_resources/meningioma-risk-classifier:v0.1.0; image manifests and extraction hashes preserved.

These are retrospective research analyses. Molecular-class concordance, within-cohort clinical discrimination, calibrated absolute risk and clinical utility are separate questions; this release establishes neither calibration nor prospective clinical utility. No author correspondence was sent.
''')
print('REPORT COMPLETE',len(completion),'profiles',passed,'passed checks')
