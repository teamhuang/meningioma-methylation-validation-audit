from pypdf import PdfReader

r = PdfReader(r"E:\脑膜瘤npj\data\suppl\natmed2024_suppl.pdf")
print("pages:", len(r.pages))
for i in range(len(r.pages)):
    t = (r.pages[i].extract_text() or "").strip()
    first = t.split("\n")[0:3]
    print(f"p{i+1:02d}: " + " | ".join(x.strip()[:90] for x in first[:2]))
