"""Static schematics (dark 3b1b style) for the report:
  slide_chromatids_cell.png  — chromatids in one G2 nucleus + the 3 contact types
  slide_rescue.png           — the 'rescuable if anchored' rule
Run: python make_cell_schematic.py"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, FancyArrowPatch, Circle

BG='#0e1117'; FG='#e6edf3'; BLUE='#58a6ff'; TEAL='#2dd4bf'; YEL='#f2cc60'; RED='#f85149'
GREY='#6e7681'; GREEN='#3fb950'; PURP='#bc8cff'; CHROM='#3a4a63'
plt.rcParams.update({'font.family':'DejaVu Sans','text.color':FG})


def sister(ax, cx, y0, y1, brdu_side, label):
    """draw one chromatid: a chromatin bar with a yellow BrdU (nascent) strand on one side."""
    w=0.30
    ax.add_patch(FancyBboxPatch((cx-w/2,y0),w,y1-y0,boxstyle='round,pad=0.02,rounding_size=0.15',
                                fc=CHROM,ec=BLUE,lw=1.2))
    bx=cx-w/2-0.06 if brdu_side=='L' else cx+w/2+0.06
    ax.plot([bx,bx],[y0+0.1,y1-0.1],color=YEL,lw=3,solid_capstyle='round')
    ax.text(cx,y0-0.22,label,ha='center',va='top',color=FG,fontsize=9)


def contact(ax,p,q,color,label,rad=0.3,lx=None,ly=None):
    ax.add_patch(Circle(p,0.07,color=color,zorder=6)); ax.add_patch(Circle(q,0.07,color=color,zorder=6))
    ax.add_patch(FancyArrowPatch(p,q,connectionstyle=f'arc3,rad={rad}',color=color,lw=2.2,arrowstyle='-',zorder=5))
    if label: ax.text(lx if lx is not None else (p[0]+q[0])/2, ly if ly is not None else (p[1]+q[1])/2,
                      label,color=color,fontsize=8.5,ha='center',va='center')


def cell_schematic():
    fig,ax=plt.subplots(figsize=(8.4,5.0)); ax.set_xlim(0,10); ax.set_ylim(0,6); ax.axis('off'); ax.set_facecolor(BG); fig.set_facecolor(BG)
    ax.add_patch(Ellipse((5,3),9.2,5.2,fc='#12161c',ec='#30363d',lw=1.4))
    ax.text(5,5.6,'One G2 nucleus: replicated chromatids and the three contact types',ha='center',fontsize=12)
    # homolog A (two sisters), homolog B (two sisters) — identical sequence, different molecules
    sister(ax,3.0,1.5,4.2,'R','pW'); sister(ax,3.7,1.5,4.2,'L','pC')
    ax.text(3.35,4.35,'homolog A',ha='center',color=FG,fontsize=9)
    sister(ax,6.3,1.5,4.2,'R','pW'); sister(ax,7.0,1.5,4.2,'L','pC')
    ax.text(6.65,4.35,'homolog B',ha='center',color=FG,fontsize=9)
    ax.add_patch(Circle((0.9,2.0),0.08,color=YEL)); ax.text(1.05,2.0,'BrdU (nascent strand)',va='center',color=YEL,fontsize=8)
    ax.text(0.9,1.55,'blue = chromatid;  each sister has BrdU on a defined strand',color=GREY,fontsize=7.5,va='center')
    # contacts
    contact(ax,(2.86,2.2),(2.86,3.4),TEAL,'cis-sister',rad=0.55,lx=2.0,ly=2.8)         # within pW
    contact(ax,(3.14,2.7),(3.56,2.7),RED,'trans-sister',rad=-0.9,lx=3.35,ly=2.0)        # pW <-> pC (same homolog)
    contact(ax,(3.85,3.7),(6.15,3.7),PURP,'inter-homolog (confound)',rad=-0.35,ly=4.75) # A <-> B
    fig.tight_layout(); fig.savefig('slide_chromatids_cell.png',dpi=120,facecolor=BG); plt.close(fig)
    print('wrote slide_chromatids_cell.png')


def rescue_schematic():
    fig,ax=plt.subplots(figsize=(8.4,3.6)); ax.set_xlim(0,10); ax.set_ylim(0,4); ax.axis('off'); ax.set_facecolor(BG); fig.set_facecolor(BG)
    ax.text(5,3.7,'"Rescuable if anchored": a BrdU+ neighbour proves the stretch replicated in the pulse',ha='center',fontsize=11)
    # one chromatid as a track
    ax.add_patch(FancyBboxPatch((0.8,2.2),8.4,0.35,boxstyle='round,pad=0.02,rounding_size=0.1',fc=CHROM,ec=BLUE,lw=1.2))
    ax.text(0.8,1.95,'one chromatid (identity continuous — no SCE)',color=GREY,fontsize=8,va='top')
    # anchor BrdU+ monomer and a nearby BrdU- monomer
    ax.add_patch(Circle((3.2,2.75),0.14,color=YEL,zorder=6)); ax.text(3.2,3.15,'BrdU+ anchor',ha='center',color=YEL,fontsize=9)
    ax.add_patch(Circle((4.5,2.75),0.14,color=GREY,zorder=6)); ax.text(4.5,3.15,'BrdU−',ha='center',color=GREY,fontsize=9)
    ax.annotate('',xy=(4.5,2.4),xytext=(3.2,2.4),arrowprops=dict(arrowstyle='<->',color=GREEN))
    ax.text(3.85,2.12,'≤ 100 kb',ha='center',color=GREEN,fontsize=8)
    ax.text(3.85,1.6,'→ rescued: same chromatid,\nsister inferred from its read strand',ha='center',va='top',color=GREEN,fontsize=8.5)
    # lone far BrdU- monomer -> discard
    ax.add_patch(Circle((8.4,2.75),0.14,color=GREY,zorder=6))
    ax.text(8.4,3.15,'BrdU− (no anchor)',ha='center',color=GREY,fontsize=8)
    ax.text(8.4,1.6,'could be unreplicated\n→ discard',ha='center',va='top',color=RED,fontsize=8)
    fig.tight_layout(); fig.savefig('slide_rescue.png',dpi=120,facecolor=BG); plt.close(fig)
    print('wrote slide_rescue.png')


if __name__=='__main__':
    cell_schematic(); rescue_schematic()
