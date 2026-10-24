import base64
import re
from pathlib import Path

md_path = Path('c:/Users/Lenovo/Desktop/PPS_APP/KET_QUA_BAI_TAP_MI3042.md')
content = md_path.read_text(encoding='utf-8')

# Replace ![caption](path) with embedded base64 img tags
def replace_img(match):
    alt = match.group(1)
    img_path = Path(match.group(2))
    if img_path.exists():
        data = base64.b64encode(img_path.read_bytes()).decode('ascii')
        return f'<div class="img-box"><img src="data:image/png;base64,{data}" alt="{alt}"><p class="caption">{alt}</p></div>'
    return match.group(0)

html_body = re.sub(r'!\[(.*?)\]\((.*?)\)', replace_img, content)

lines = html_body.splitlines()
out_lines = []
in_table = False
for line in lines:
    line_s = line.strip()
    if not line_s:
        if in_table:
            out_lines.append('</table>')
            in_table = False
        continue
    if line_s.startswith('# '):
        out_lines.append(f'<h1>{line_s[2:]}</h1>')
    elif line_s.startswith('## '):
        out_lines.append(f'<h2>{line_s[3:]}</h2>')
    elif line_s.startswith('### '):
        out_lines.append(f'<h3>{line_s[4:]}</h3>')
    elif line_s.startswith('#### '):
        out_lines.append(f'<h4>{line_s[5:]}</h4>')
    elif line_s.startswith('|') and '|' in line_s[1:]:
        if '---' in line_s:
            continue
        cells = [c.strip() for c in line_s.split('|')[1:-1]]
        if not in_table:
            out_lines.append('<table border="1" class="data-table">')
            out_lines.append('<tr>' + ''.join(f'<th>{c}</th>' for c in cells) + '</tr>')
            in_table = True
        else:
            out_lines.append('<tr>' + ''.join(f'<td>{c}</td>' for c in cells) + '</tr>')
    elif line_s.startswith('* ') or line_s.startswith('- '):
        out_lines.append(f'<li>{line_s[2:]}</li>')
    elif line_s.startswith('$$') and line_s.endswith('$$'):
        out_lines.append(f'<div class="math">{line_s}</div>')
    elif line_s.startswith('<div class="img-box"'):
        out_lines.append(line_s)
    else:
        out_lines.append(f'<p>{line_s}</p>')

if in_table:
    out_lines.append('</table>')

html_template = f'''<!DOCTYPE html>
<html lang="vi">
<head>
<meta charset="UTF-8">
<title>Báo cáo Bài tập Phương Pháp Số MI3042</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/katex.min.css">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/katex.min.js"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/contrib/auto-render.min.js" onload="renderMathInElement(document.body, {{delimiters: [{{left: '$$', right: '$$', display: true}}, {{left: '$', right: '$', display: false}}]}});"></script>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; line-height: 1.6; color: #1e293b; max-width: 1080px; margin: 40px auto; padding: 0 24px; background-color: #f8fafc; }}
h1 {{ color: #0f172a; border-bottom: 2px solid #2563eb; padding-bottom: 12px; }}
h2 {{ color: #1d4ed8; border-bottom: 1px solid #cbd5e1; padding-bottom: 8px; margin-top: 36px; }}
h3 {{ color: #0369a1; margin-top: 20px; }}
.img-box {{ text-align: center; margin: 24px 0; background: #fff; padding: 16px; border-radius: 8px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1); border: 1px solid #e2e8f0; }}
.img-box img {{ max-width: 100%; height: auto; border-radius: 4px; }}
.caption {{ font-size: 0.9em; color: #64748b; margin-top: 8px; font-style: italic; font-weight: 500; }}
table.data-table {{ border-collapse: collapse; width: 100%; margin: 16px 0; background: #fff; border-radius: 6px; overflow: hidden; }}
table.data-table th, table.data-table td {{ border: 1px solid #cbd5e1; padding: 8px 12px; text-align: center; }}
table.data-table th {{ background: #e0f2fe; color: #0369a1; }}
code {{ background: #e2e8f0; padding: 2px 6px; border-radius: 4px; font-family: Consolas, monospace; }}
p {{ margin: 8px 0; }}
li {{ margin: 4px 0; }}
</style>
</head>
<body>
{'\\n'.join(out_lines)}
</body>
</html>
'''

Path('c:/Users/Lenovo/Desktop/PPS_APP/KET_QUA_BAI_TAP_MI3042.html').write_text(html_template, encoding='utf-8')
Path('c:/Users/Lenovo/Desktop/KET_QUA_BAI_TAP_MI3042.html').write_text(html_template, encoding='utf-8')
print('Exported KET_QUA_BAI_TAP_MI3042.html to Desktop and PPS_APP successfully!')
