"""Type IIS (AlwI) problem schematic -> slide_typeIIS.png (dark 3b1b style)."""
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle

BG='#0e1117'; FG='#e6edf3'; GREEN='#3fb950'; GREY='#6e7681'; YEL='#f2cc60'; RED='#f85149'; BLUE='#58a6ff'; CHROM='#3a4a63'
plt.rcParams.update({'font.family':'DejaVu Sans','text.color':FG})

fig, ax = plt.subplots(figsize=(9.2, 5.8)); ax.set_xlim(0,10); ax.set_ylim(0,10); ax.axis('off'); ax.set_facecolor(BG); fig.set_facecolor(BG)

def dna(x0,x1,y,ec=BLUE,h=0.32):
    ax.add_patch(FancyBboxPatch((x0,y-h/2),x1-x0,h,boxstyle='round,pad=0.02,rounding_size=0.08',fc=CHROM,ec=ec,lw=1.4))
def ggatc(x,y):
    ax.add_patch(FancyBboxPatch((x-0.33,y-0.2),0.66,0.4,boxstyle='round,pad=0.01',fc='#2a2f12',ec=YEL,lw=1.3))
    ax.text(x,y,'GGATC',ha='center',va='center',color=YEL,fontsize=7)

ax.text(5,9.5,'Type IIS problem: AlwI cuts 4–5 nt AWAY from its site (GGATC 4/5)',ha='center',fontsize=13)

# --- top: one cut makes two asymmetric ends ---
ax.text(0.5,8.3,'one cut:',color=FG,fontsize=10,va='center')
dna(2.0,5.0,8.3,ec=GREEN); ggatc(4.2,8.3)
ax.plot([5.0,5.0],[7.95,8.65],color=RED,ls=(0,(2,2)),lw=1.6)
ax.text(5.0,8.95,'✂ cut (4–5 nt 3′ of site)',ha='center',color=RED,fontsize=8)
dna(5.0,8.0,8.3,ec=GREY)
ax.text(3.5,7.75,'site-RETAINING end\n(GGATC kept)',ha='center',color=GREEN,fontsize=8)
ax.text(6.5,7.75,'site-LESS end\n(no GGATC)',ha='center',color=GREY,fontsize=8)

# --- two junction outcomes ---
ax.text(5,6.4,'when two ends ligate in a concatemer:',ha='center',color=FG,fontsize=10)
# case A: a site-retaining end at the junction -> splittable
dna(1.0,2.9,5.2,ec=GREEN); ggatc(2.6,5.2); dna(2.9,4.8,5.2,ec=GREY)
ax.add_patch(Circle((2.9,5.2),0.06,color=FG,zorder=6))
ax.plot([2.9,2.9],[4.85,5.55],color=GREEN,ls=(0,(2,2)),lw=1.4)
ax.text(2.9,5.95,'GGATC at junction',ha='center',color=GREEN,fontsize=8)
ax.text(2.9,4.5,'→ motif-digest SPLITS  ✓\ntwo monomers',ha='center',va='top',color=GREEN,fontsize=9)
# case B: both ends site-less -> fused
dna(5.6,7.0,5.2,ec=GREY); dna(7.0,8.9,5.2,ec=GREY)
ax.add_patch(Circle((7.0,5.2),0.06,color=FG,zorder=6))
ax.text(7.0,5.95,'no GGATC at junction',ha='center',color=RED,fontsize=8)
ax.text(7.0,4.5,'→ CANNOT split  ✗\nfused: 1 monomer, partner soft-clipped & lost',ha='center',va='top',color=RED,fontsize=9)

# --- stats + recovery ---
ax.add_patch(FancyBboxPatch((0.8,1.1),8.4,1.6,boxstyle='round,pad=0.03,rounding_size=0.12',fc='#10261c',ec='#1d573b',lw=1.2))
ax.text(5,2.35,'Each cut → 1 site-retaining + 1 site-less end  ⇒  25% of junctions are site-less ↔ site-less',
        ha='center',color=FG,fontsize=10)
ax.text(5,1.75,'⇒ ~25% of AlwI–AlwI junctions can’t be split by the motif digest (hits multi-way contacts hardest)',
        ha='center',color=RED,fontsize=9.5)
ax.text(5,1.3,'RECOVER them from the alignment: a fused monomer is chimeric → supplementary alignment → split at the breakpoint',
        ha='center',color=GREEN,fontsize=9.5)

fig.tight_layout(); fig.savefig('slide_typeIIS.png',dpi=120,facecolor=BG)
print('wrote slide_typeIIS.png')
