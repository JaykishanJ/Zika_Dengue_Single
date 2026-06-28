from PIL import Image, ImageDraw, ImageFont
import os

img_path = r"E:\Zika_Enrichment\Final_95_Genes_Published_Network.png"

pass  # Execution logging removed for final release
img_orig = Image.open(img_path).convert("RGBA")

# Expand the canvas by 600 pixels on the bottom to give the legend plenty of space!
orig_width, orig_height = img_orig.size
new_width = orig_width
new_height = orig_height + 600

# Create new white canvas and paste original image at the top
img = Image.new("RGBA", (new_width, new_height), (255, 255, 255, 255))
img.paste(img_orig, (0, 0))

draw = ImageDraw.Draw(img)

# Try to load Arial font
try:
    font_title = ImageFont.truetype(r"C:\Windows\Fonts\arialbd.ttf", 60) # Bold
    font_text = ImageFont.truetype(r"C:\Windows\Fonts\arial.ttf", 50)
except Exception as e:
    pass  # Execution logging removed for final release
    font_title = ImageFont.load_default()
    font_text = ImageFont.load_default()

box_width = 750
box_height = 400
margin = 150

# Place legend in the newly created empty space on the far LEFT side, 200 pixels below the network!
x0 = margin
y0 = orig_height + 200 
x1 = x0 + box_width
y1 = y0 + box_height

# Draw legend box with black border
draw.rectangle([x0, y0, x1, y1], fill=(255, 255, 255, 255), outline=(0, 0, 0, 255), width=4)

# Draw Legend Title
draw.text((x0 + 50, y0 + 40), "Disease Specificity", fill=(0, 0, 0, 255), font=font_title)

# Legend items
items = [
    ("#FF9999", "Dengue Specific"),
    ("#99CCFF", "Zika Specific"),
    ("#CC99FF", "Shared (Both)")
]

y_offset = y0 + 130
for color_hex, text in items:
    # Draw circle/node representation
    draw.ellipse([x0 + 50, y_offset, x0 + 110, y_offset + 60], fill=color_hex, outline=(64, 64, 64, 255), width=3)
    # Draw text
    draw.text((x0 + 140, y_offset + 5), text, fill=(0, 0, 0, 255), font=font_text)
    y_offset += 85

out_path = r"E:\Zika_Enrichment\Final_Publication_Image_With_Legend.png"
img.save(out_path)
pass  # Execution logging removed for final release
