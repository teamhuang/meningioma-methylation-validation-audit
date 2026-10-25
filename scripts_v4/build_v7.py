import os, json
ROOT = r'E:\脑膜瘤npj'
p = os.path.join(ROOT, 'manuscript_npjDM_v7.md')
s = open(p, encoding='utf-8').read()

REPS = [
 # T1 title -> 12 words, no punctuation
 ("\ufeff# Identity, training independence and selective prediction shape the transportability of a meningioma DNA methylation classifier",
  "\ufeff# Auditing the transportability of a meningioma DNA methylation classifier across independent cohorts"),
 # T2 Results counterfactual 11-of-37
 ("0.909 against the training-label map on the 11 of its 37 profiles that are genuinely non-training — a positive-looking agreement figure derived from a set that is 54% training data.",
  "0.909 against the training-label map computed on the 11 of its 17 non-training profiles that carry known, non-tied labels (the remaining six non-training profiles lack an evaluable label) — a positive-looking agreement figure derived from a set that is 54% training data."),
 # T3 Table 6 wording
 ("(0.909 on its 11 non-training profiles)",
  "(0.909 on the 11 of 17 non-training profiles with known non-tied labels)"),
 # T5 new Results paragraph on rejection vs outcome
 ("Higher selective agreement at stricter thresholds comes with loss of coverage, particularly in the hypermitotic reference group, and should be interpreted alongside class-specific retention.",
  "Higher selective agreement at stricter thresholds comes with loss of coverage, particularly in the hypermitotic reference group, and should be interpreted alongside class-specific retention.\n\n"
  "To test whether preferential rejection of the hypermitotic reference group corresponds to worse patient outcomes, we compared observed recurrence between donors whose specimen was rejected (tied vote or confidence < 0.75) and donors whose specimen was retained, at the donor level and adjusted for WHO grade (Supplementary Data 6). Rejection was strongly class-dependent — 48.2% of published hypermitotic-reference specimens were rejected versus 6.7% of Merlin-intact-reference specimens in Hong Kong (34.4% versus 9.7% in UCSF; 28.6% versus 5.6% in Tübingen) — but rejected patients were not at consistently higher observed recurrence risk. The grade-adjusted hazard ratio for rejected versus retained donors was 1.27 (95% CI 0.68–2.36; p = 0.46) in Hong Kong, 1.04 (0.56–1.93; p = 0.89) in UCSF and 0.62 (0.30–1.28; p = 0.20) in Tübingen, with directions that disagree across cohorts. We therefore report the coverage consequence as a representativeness and reporting problem, not as a demonstrated increase in patient risk: non-random rejection can bias the evaluated population and the apparent operating characteristics without, on these data, being shown to harm individual patients."),
 # T6 Discussion sentence
 ("These published labels are comparison references rather than confirmed clinical pathology truth, and the retention pattern alone does not establish how a clinical decision based on the classifier would affect patient outcomes.",
  "These published labels are comparison references rather than confirmed clinical pathology truth, and the retention pattern alone does not establish how a clinical decision based on the classifier would affect patient outcomes. We tested this directly and found that rejected patients were not at consistently higher observed recurrence risk across the three cohorts (grade-adjusted hazard ratios spanning 0.62 to 1.27, none significant; Supplementary Data 6), so the defensible conclusion is a coverage and reporting requirement rather than a demonstrated clinical harm."),
 # Methods subsection before Ethics approval
 ("## Ethics approval",
  "### Rejection and observed recurrence\n\nEach specimen was labelled rejected (tied maximum votes or classifier confidence < 0.75) or retained, and this status was joined to the cohort's own recurrence endpoint. The analysis unit was the donor: donors whose specimens did not share a single rejection status were dropped and counted. Observed recurrence was compared descriptively (Fisher exact) and by a Cox model adjusted for ordinal WHO grade, one record per donor, with a robust variance. This is an association test on retrospective, cohort-specific endpoints, reported alongside rather than instead of the coverage analysis, and it cannot support a causal claim about abstention.\n\n## Ethics approval"),
]
for old, new in REPS:
    c = s.count(old)
    assert c == 1, f'count={c} for {old[:55]!r}'
    s = s.replace(old, new)
tmp = p + '.tmp'
with open(tmp, 'w', encoding='utf-8') as fh: fh.write(s); fh.flush(); os.fsync(fh.fileno())
os.replace(tmp, p)

# T7: Supplementary Data 6 from the JSON
d = json.load(open(os.path.join(ROOT, 'results_v4/rejection_vs_outcome_v7.json'), encoding='utf-8'))
rows = []
for coh in ['hku', 'ucsf', 'tubingen']:
    r = d[coh]
    rows.append({
      'cohort': coh, 'donors_analysed': r['donors'], 'donors_dropped_ambiguous': r['donors_dropped_ambiguous'],
      'n_rejected': r['n_rejected'], 'n_retained': r['n_retained'],
      'events_rejected': r['events_rejected'], 'events_retained': r['events_retained'],
      'event_rate_rejected': r['event_rate_rejected'], 'event_rate_retained': r['event_rate_retained'],
      'event_fisher_p': r.get('event_fisher_p'), 'cox_HR_rejected_vs_retained': r.get('cox_HR_rejected_vs_retained'),
      'cox_HR_ci_low': (r.get('cox_HR_ci') or [None,None])[0], 'cox_HR_ci_high': (r.get('cox_HR_ci') or [None,None])[1],
      'cox_p': r.get('cox_p'),
      'rejected_frac_hypermitotic': r['rejected_fraction_by_group'].get('Hypermitotic'),
      'rejected_frac_merlin_intact': r['rejected_fraction_by_group'].get('Merlin-intact'),
      'interpretation': 'grade-adjusted, donor-level; association only; no consistent higher risk in rejected donors'
    })
import csv
sd6 = os.path.join(ROOT, 'npjDM_submission_candidate_v6/Supplementary_Data_6_rejection_outcome.csv')
with open(sd6, 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('v7 edits applied; Supplementary Data 6 written:', sd6)
print('title words:', len(s.splitlines()[0].lstrip('#').replace('\ufeff','').split()))
for k in ['11 of its 37 profiles that are genuinely non-training', '11 non-training profiles)']:
    print('residual stale?', k, '->', k in s)
