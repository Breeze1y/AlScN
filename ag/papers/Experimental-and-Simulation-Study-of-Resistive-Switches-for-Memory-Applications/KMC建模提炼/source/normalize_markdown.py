from pathlib import Path
import re

report = Path(__file__).resolve().parents[1] / "KMC搭建过程与论文研究拆解.md"
text = report.read_text(encoding="utf-8")
lines = []
count = 0
for line in text.splitlines():
    pieces = line.split("**")
    if len(pieces) > 1 and len(pieces) % 2 == 1:
        for index in range(1, len(pieces), 2):
            pieces[index] = pieces[index].strip()
        for index in range(2, len(pieces), 2):
            if re.match(r"[A-Za-z\u3400-\u9fff]", pieces[index]):
                pieces[index] = " " + pieces[index]
                count += 1
        line = "**".join(pieces)
    lines.append(line)
text = "\n".join(lines) + "\n"
report.write_text(text, encoding="utf-8")
print(f"已规范 {count} 处加粗与中文相接的边界")
