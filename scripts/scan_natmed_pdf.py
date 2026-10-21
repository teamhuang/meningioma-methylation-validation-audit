from pypdf import PdfReader

r = PdfReader(r"E:\脑膜瘤npj\data\suppl\natmed2024_suppl.pdf")
print("pages:", len(r.pages))
t0 = r.pages[0].extract_text() or ""
print("PAGE1:", t0[:500].replace("\n", " | "))
hits = []
for i in range(len(r.pages)):
    t = r.pages[i].extract_text() or ""
    tl = t.lower()
    score = 0
    if "supplementary table" in tl: score += 1
    if "recurrence" in tl or "progression" in tl: score += 1
    if "follow-up" in tl or "follow up" in tl or "months" in tl: score += 1
    if "sample" in tl and ("id" in tl or "patient" in tl): score += 1
    if score >= 2:
        hits.append((i + 1, score, t[:150].replace("\n", " | ")))
print("\nCANDIDATE PAGES:")
for h in hits[:20]:
    print(" page", h[0], "score", h[1], ":", h[2])
