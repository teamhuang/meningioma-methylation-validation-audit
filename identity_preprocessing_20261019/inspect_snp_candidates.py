from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
p=pd.read_csv(ROOT/'results/SNP_pairwise_identity_metrics.csv')
clin=pd.read_csv(BASE/'results/gse304094_primary_survival.csv')
for _,r in p[p.same_donor_candidate].iterrows():
    print(r[['GSM1','GSM2','title1','title2','common_highconfidence_SNPs','discordant_SNPs','SNP_beta_correlation']].to_dict())
    print(clin[clin.Internal_ID.isin([r.title1,r.title2])][['Internal_ID','Sex','Age [Years]','tumor_setting','Localisation','Follow-Up [Months]','Status (0=censored, 1=progression)']].to_string(index=False))
print('Noncandidate smallest SNP discordance\n',p[~p.same_donor_candidate].sort_values('discordance_fraction').head(10)[['title1','title2','common_highconfidence_SNPs','discordant_SNPs','discordance_fraction']].to_string(index=False))
bc=pd.read_csv(BASE/'results/clinical_parsed.csv')
for _,r in p[(p.relaxed_candidate)&(p.cohort1=='baylor')&(p.cohort2=='baylor')].iterrows():
    names=[r.title1.replace(' (methylation)',''),r.title2.replace(' (methylation)','')]
    print('BAYLOR CANDIDATE',r[['title1','title2','SNP_beta_correlation','discordant_SNPs','common_highconfidence_SNPs']].to_dict())
    print(bc[bc.Sample.isin(names)][['Sample','Gender','age','Location','fu_months','grade','gtr','rec','menG']].to_string(index=False))
