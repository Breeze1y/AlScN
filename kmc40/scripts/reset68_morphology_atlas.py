"""Export endpoint atlas and a fractional-Ag view of two RESET ruptures."""
import argparse
import csv
import io
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap, to_rgb
from matplotlib.patches import Rectangle
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '.venv'))
from postprocess_100_cycles import COLORS, PHASES, morphology_page, span

CURRENT = ROOT / 'replay/reset68/continuous'
REPORT = ROOT / 'replay/cycles100/report'
TEAL = '#006877'
INK = '#253039'


def original(c):
    return (ROOT / ('results/cycle_reset_fit40/mid_clock_set' if c <= 5 else 'replay/cycles100')) / f'cycle_{c:02d}'


def endpoints():
    atlas = []
    fixed_vn = None
    for cycle in range(1, 101):
        record = {}
        for phase in PHASES:
            new = CURRENT / f'cycle_{cycle:02d}' / phase / 'states.npz'
            with np.load(new) as a, np.load(original(cycle) / phase / 'states.npz') as b:
                ids = a['metal_id'][-1]
                np.testing.assert_array_equal(ids, b['metal_id'][-1])
                if fixed_vn is None:
                    fixed_vn = a['VN']
                np.testing.assert_array_equal(a['VN'], fixed_vn)
                np.testing.assert_array_equal(a['VN'], b['VN'])
                if phase == 'RESET':
                    units = a['Ag_units'][-1]
                    np.testing.assert_array_equal(units, b['Ag_units'][-1])
                    np.testing.assert_array_equal(units, (ids > 0) * int(a['subunits']))
                record[phase] = {'end': ids > 0}
        atlas.append(record)
    return atlas, fixed_vn


def rupture_frames(cycle):
    folder = CURRENT / f'cycle_{cycle:02d}' / 'RESET'
    with (folder/'events.csv').open(encoding='utf-8-sig') as stream:
        event = next(r for r in csv.DictReader(stream)
                     if r['connected_before'] == 'True' and r['connected_after'] == 'False')
    with np.load(folder/'states.npz') as a, np.load(original(cycle)/'RESET/states.npz') as b:
        new, old = a['Ag_units'], b['Ag_units']
        n = int(a['subunits'])
        assert n == int(b['subunits']) and n == 16
        np.testing.assert_array_equal(a['voltage_V'], b['voltage_V'])
        first = int(event['point'])
        assert a['metal_connected'][first-1] and not a['metal_connected'][first]
        delta = new.astype(int) - old.astype(int)
        # Same-bias paired comparison at the largest outward total discrepancy.
        point = first + int(np.abs(delta[first:101]).sum(axis=(1, 2)).argmax())
        largest = delta[point].max()
        assert largest > 0
        choose = np.argwhere(delta[point] == largest)
        y, x = choose[np.argmin(old[point][tuple(choose.T)])]
        by, bx = int(event['y']), int(event['x'])
        yc, xc = (y+by)//2, (x+bx)//2
        roi = (slice(yc-3,yc+4), slice(xc-3,xc+4))
        assert old[-1,y,x] == new[-1,y,x]
        indices = [first-1, first, point, point, len(new)-1]
        frames = [new[first-1], new[first], new[point], old[point], new[-1]]
        return dict(cycle=cycle, roi=roi, marker=(int(y),int(x)), rupture=float(event['time']),
                    break_v=float(event['voltage_V']), frames=frames, n=n, vn=a['VN'],
                    time=a['elapsed_units'][indices], voltage=a['voltage_V'][indices])


def zoom_colors(units, n, vn):
    color = np.ones((*units.shape, 3))
    color[vn] = to_rgb(COLORS[1])
    fraction = units[units > 0, None] / n
    pale, full = np.array(to_rgb('#ffe4c3')), np.array(to_rgb(COLORS[2]))
    color[units > 0] = pale + fraction * (full-pale)
    return color


def overview(data):
    fig = plt.figure(figsize=(12, 7.45), facecolor='white')
    font = {'fontfamily': 'Microsoft YaHei'}
    def caption(x, y, txt, **kwargs):
        return fig.text(x, y, txt, **({'color':INK} | font | kwargs))
    caption(.045,.934,'银细丝断桥后的残余占据量',fontsize=22,weight='bold')
    caption(.045,.889,'连续 100 圈的 SET 和 RESET 终态与原 KMC 逐格相同。断桥附近尚未溶尽的 Ag，在两种退缩规则下可停留至不同采样点。',fontsize=10.5)
    titles = ['全膜 · SET 终态','断桥前一采样点','首次断桥采样点','共享退缩 · 对照点','原 KMC · 同一时刻','共享退缩 · RESET 终态']
    for line, cycle in enumerate((68,100)):
        record = rupture_frames(cycle)
        yy, xx = record['roi']
        y,x = record['marker']
        top = .815 - .344 * line
        caption(.045,top, f'第 {cycle} 圈',fontsize=14,weight='bold')
        caption(.164,top+.002, f"断桥事件发生在 t = {record['rupture']:.2f}，V = {record['break_v']:.3f} V",fontsize=10)
        for j, title in enumerate(titles):
            left=.043+j*.156
            ax=fig.add_axes((left,top-.213,.139,.168))
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_edgecolor('#abb8be')
                spine.set_linewidth(.65)
            ax.set_title(title,fontsize=9.4,pad=4,**font)
            if j==0:
                ag=data[cycle-1]['SET']['end']
                ax.imshow(np.where(span(ag),3,np.where(ag,2,record['vn'])),cmap=ListedColormap(COLORS),
                          vmin=0,vmax=3,interpolation='nearest')
                ax.add_patch(Rectangle((xx.start-.5, yy.start-.5),7,7,fill=False,edgecolor='#111b21',lw=1.4))
                caption(left+.0695,top-.233,'全膜 50 × 40 nm',fontsize=8.4,ha='center')
                continue
            idx=j-1
            u=record['frames'][idx]
            rgb=zoom_colors(u[yy,xx],record['n'],record['vn'][yy,xx])
            ax.imshow(rgb,interpolation='nearest')
            ax.add_patch(Rectangle((x-xx.start-.5,y-yy.start-.5),1,1,edgecolor=TEAL,lw=1.6,fill=False))
            t,v=record['time'][idx],record['voltage'][idx]
            caption(left+.0695,top-.233,f't = {t:.0f}   {v:+.3f} V' if v else f't = {t:.0f}   0 V',
                    fontsize=8.5,ha='center')
            caption(left+.0695,top-.263,f"框内 Ag   {u[y,x]}/{record['n']}",fontsize=9.6,ha='center',color=TEAL,weight='bold')
    caption(.045,.144,'读图',fontsize=11,weight='bold')
    caption(.096,.144,'全膜黑框对应 3.5 × 3.5 nm 放大区。青框追踪同一 0.5 nm 格点。格点残余分数取自每次采样结束时的 Ag_units / 16。',fontsize=9.5)
    caption(.045,.103,'放大区颜色',fontsize=9.5,weight='bold')
    vals=[(0,'空基质'),(1,'1/16'),(4,'4/16'),(8,'8/16'),(12,'12/16'),(16,'16/16')]
    for i,(value,label) in enumerate(vals):
        x=.134+i*.097
        rgb=zoom_colors(np.array([[value]]),16,np.array([[False]]))[0,0]
        fig.add_artist(Rectangle((x,.102),.011,.014,facecolor=rgb,edgecolor='#8b969c',lw=.55))
        caption(x+.016,.101,label,fontsize=8.2)
    fig.add_artist(Rectangle((.719,.102),.011,.014,facecolor=COLORS[1],edgecolor='#8b969c',lw=.55))
    caption(.735,.101,'氮空位',fontsize=8.2)
    caption(.045,.061,'t 从该圈 RESET 起点计，以有效驻留时间为单位。对照点取外扫中两种模型占据差异总量最大的位置。',fontsize=8.7)
    caption(.045,.031,'第 2 页为代表圈终态。第 3–12 页按圈号列出全部 200 个 SET / RESET 终态，每页 10 圈。',fontsize=8.7)
    return fig


def png(fig):
    stream=io.BytesIO()
    fig.savefig(stream,format='png',dpi=200)
    stream.seek(0)
    return Image.open(stream).convert('RGB')


def main(output=REPORT):
    output=Path(output)
    output.mkdir(parents=True,exist_ok=True)
    data,vn=endpoints()
    with PdfPages(output/'morphology_100cycles.pdf') as book:
        front=overview(data)
        book.savefig(front,dpi=240)
        preview=png(front)
        plt.close(front)
        selected = np.rint(np.linspace(1,100,5)).astype(int).tolist()
        milestones=morphology_page(data,vn,[selected],'Shared residual retreat | Selected cycle endpoints')
        book.savefig(milestones,dpi=240)
        strip=png(milestones)
        plt.close(milestones)
        strip=strip.resize((preview.width,round(strip.height*preview.width/strip.width)),Image.Resampling.LANCZOS)
        combined=Image.new('RGB',(preview.width,preview.height+strip.height),'white')
        combined.paste(preview,(0,0)); combined.paste(strip,(0,preview.height))
        combined.save(output/'morphology_overview.png')
        for start in range(1,101,10):
            groups=[range(start,start+5),range(start+5,start+10)]
            fig=morphology_page(data,vn,groups,f'Ag morphology at cycle endpoints {start}–{start+9}')
            book.savefig(fig,dpi=240)
            plt.close(fig)
    print('Compared Ag identities, occupancies and VN at 200 endpoints; rendered 12 morphology pages.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=REPORT)
    main(parser.parse_args().output)
