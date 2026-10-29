"""Draw the exact code-geometry schematic; no data or simulation."""
from pathlib import Path
from itertools import combinations
from reportlab.pdfgen import canvas
from reportlab.lib.colors import black, white, Color
out=Path(__file__).resolve().parent/'code-geometry.pdf'
c=canvas.Canvas(str(out),pagesize=(468,122.4))
c.setTitle('Four-gene target-code geometry')
c.setFont('Helvetica',10)
def text(x,y,s):
 c.setFillColor(black); c.drawCentredString(x,y,s)
def edges(points,pairs,gray=0):
 c.setStrokeColor(Color(gray,gray,gray)); c.setLineWidth(.8)
 for i,j in pairs: c.line(*points[i],*points[j])
def point(x,y):
 c.setFillColor(black); c.circle(x,y,2,stroke=0,fill=1)
p=[(45,37),(175,37),(45,86),(175,86)]
edges(p,[(0,1),(0,2),(1,3),(2,3)])
for (x,y),s in zip(p,['00','10','01','11']):
 point(x,y); text(x,y+(10 if y>50 else -14),s)
text(110,111,'Two bits: order 3 fails')
text(110,5,'00 + 11 = 10 + 01')
p=[(286,37),(416,37),(351,86),(351,57)]
edges(p,combinations(range(4),2),.45)
for (x,y),s,dx,dy in zip(p,['000','110','101','011'],[-7,7,0,23],[-14,-14,10,-2]):
 point(x,y)
 tx,ty=x+dx,y+dy
 if s=='011':
  c.setFillColor(white); c.rect(tx-10,ty-2,21,12,stroke=0,fill=1)
 text(tx,ty,s)
text(351,111,'Three bits: order 3 succeeds')
text(351,5,'Real determinant = 2')
c.showPage(); c.save()
print(out)
