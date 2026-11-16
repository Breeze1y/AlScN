from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

root = Path(__file__).resolve().parents[1]
assets = root / "assets"
font = ImageFont.truetype(r"C:\Windows\Fonts\msyh.ttc", 24)
canvas = Image.new("RGB", (1120, 270), "white")
draw = ImageDraw.Draw(canvas)
for filename, x, label in [
    ("Image425-extracted.png", 20, "(a) 表面扩散慢，原地还原快"),
    ("Image426-extracted.png", 585, "(b) 表面扩散快，表面可重排"),
]:
    image = Image.open(root / "source" / filename).convert("RGB")
    width = 510
    height = round(image.height * width / image.width)
    image = image.resize((width, height), Image.Resampling.NEAREST)
    canvas.paste(image, (x, 15))
    draw.text((x, 223), label, font=font, fill="#162b40")
canvas.save(assets / "figure-2-20.png")
print("已将原论文两幅形貌子图整理为阅读版")
