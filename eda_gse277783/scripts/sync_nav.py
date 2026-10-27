from pathlib import Path

html_files = [f for f in Path(r"F:\Thesis").glob("*.html") if f.name not in ["index.html", "spatial.html", "gse277783_eda.html"]]
print(f"Checking {len(html_files)} other HTML files for nav menu sync...")

for f in html_files:
    content = f.read_text(encoding="utf-8")
    if "gse277783_eda.html" not in content:
        old_str1 = '<a href="index.html" class="">Overview</a>'
        new_str1 = '<a href="index.html" class="">Overview</a>\n      <a href="gse277783_eda.html" class="">GSE277783 EDA</a>'
        old_str2 = '<a href="index.html">Overview</a>'
        new_str2 = '<a href="index.html">Overview</a>\n      <a href="gse277783_eda.html">GSE277783 EDA</a>'
        old_str3 = '<a href="index.html" class="active">Overview</a>'
        new_str3 = '<a href="index.html" class="active">Overview</a>\n      <a href="gse277783_eda.html" class="">GSE277783 EDA</a>'
        
        if old_str1 in content:
            content = content.replace(old_str1, new_str1)
            f.write_text(content, encoding="utf-8")
            print(f"  Updated nav in {f.name}")
        elif old_str2 in content:
            content = content.replace(old_str2, new_str2)
            f.write_text(content, encoding="utf-8")
            print(f"  Updated nav in {f.name}")
        elif old_str3 in content:
            content = content.replace(old_str3, new_str3)
            f.write_text(content, encoding="utf-8")
            print(f"  Updated nav in {f.name}")
        else:
            print(f"  Could not find standard index link in {f.name}")

print("Nav menu sync complete!")
