import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { marked } from 'file:///C:/Users/31826/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/marked/lib/marked.esm.js';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const title = 'KMC 搭建过程与论文研究拆解';
const basename = 'KMC搭建过程与论文研究拆解';
const md = await fs.readFile(path.join(root, basename + '.md'), 'utf8');
const escape = (s) => s.replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;');
const toc = [];
let number = 0;
marked.use({renderer: {
  heading({tokens, depth, text}) {
    const inline = this.parser.parseInline(tokens);
    const id = 'section-' + (++number);
    if (depth === 2) toc.push({id, text});
    return `<h${depth} id="${id}">${inline}</h${depth}>\n`;
  },
  code({text, lang}) {
    if (lang === 'mermaid') return `<details class="diagram-source"><summary>查看框图的可编辑文本</summary><pre><code>${escape(text)}</code></pre></details>\n`;
    return `<pre><code>${escape(text)}</code></pre>\n`;
  },
  image({href, text}) {
    return `<a class="figure-link" href="${escape(href)}" target="_blank" rel="noopener" aria-label="放大：${escape(text)}"><img src="${escape(href)}" alt="${escape(text)}"></a>`;
  },
  table({header, rows}) {
    const row = (cells, tag) => '<tr>' + cells.map(c => `<${tag}>${this.parser.parseInline(c.tokens)}</${tag}>`).join('') + '</tr>';
    return `<div class="table-wrap"><table><thead>${row(header,'th')}</thead><tbody>${rows.map(r=>row(r,'td')).join('')}</tbody></table></div>\n`;
  }
}});
const body = marked.parse(md, {gfm: true});
const nav = toc.map(x => `<li><a href="#${x.id}">${escape(x.text)}</a></li>`).join('\n');
const html = `<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>${title}</title>
<style>
:root{color-scheme:light;--ink:#172e44;--muted:#53677b;--line:#dce4ec;--paper:#fff;--back:#f3f6f9;--blue:#265f84}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--back);color:var(--ink);font:16px/1.85 "Microsoft YaHei","Noto Sans CJK SC",sans-serif}
.sheet{max-width:1160px;margin:36px auto 70px;padding:48px 64px 64px;background:var(--paper);box-shadow:0 6px 32px #1f395310;border-radius:10px}
.eyebrow{font-size:12px;letter-spacing:2px;color:var(--muted);text-transform:uppercase;margin:0 0 14px}.reading-links{font-size:14px;color:var(--muted);margin:0 0 26px}.reading-links a{margin-right:18px}
h1{font-size:31px;line-height:1.4;letter-spacing:.2px;margin:0 0 28px}h2{font-size:24px;line-height:1.5;margin:48px 0 18px;padding-top:20px;border-top:2px solid var(--line);scroll-margin-top:20px}h3{font-size:19px;margin:30px 0 14px;scroll-margin-top:20px}
p{margin:14px 0}strong{font-weight:650}a{color:var(--blue);text-underline-offset:3px}a:hover{text-decoration-thickness:2px}li{margin:7px 0}ol,ul{padding-left:27px}blockquote{margin:20px 0;padding:12px 20px;border-left:4px solid #81a6bf;background:#f1f6fa;font-size:17px;overflow-wrap:anywhere}blockquote p{margin:3px 0}
code{font:14px/1.8 Consolas,"Microsoft YaHei",monospace;background:#f0f4f7;border-radius:3px;padding:2px 5px;overflow-wrap:anywhere}pre{background:#f0f4f7;border-radius:6px;padding:18px 20px;overflow-x:auto;max-width:100%;line-height:1.85}pre code{padding:0;background:none;white-space:pre-wrap;word-break:break-word}
.table-wrap{overflow-x:auto;margin:18px 0 24px}table{border-collapse:collapse;width:100%;font-size:14px;line-height:1.75;min-width:570px}th{font-weight:650;text-align:left;background:#eaf1f6;border-top:1px solid #cbd8e3}td,th{padding:11px 12px;vertical-align:top;border-bottom:1px solid var(--line)}tbody tr:nth-child(even){background:#f8fafc}td:first-child,th:first-child{min-width:72px}
img{display:block;max-width:100%;height:auto;margin:24px auto 30px;object-fit:contain}img[src$=".svg"]{width:100%;max-width:1000px}img[src*="figure-2-2"]{width:850px;image-rendering:auto}img[src*="equation-9"]{width:100%;padding:15px;background:#fff}img[src*="figure-2-20"]{width:100%}.figure-link{display:block}
.toc{margin:28px 0 34px;padding:16px 20px;background:#f3f7fa;border-radius:6px}.toc summary{font-weight:650;cursor:pointer}.toc ul{columns:2;column-gap:32px;list-style:none;padding-left:0;font-size:14px}.toc li{break-inside:avoid;margin:7px 0}.toc a{text-decoration:none}.diagram-source{font-size:14px;color:var(--muted);margin:12px 0}.diagram-source summary{cursor:pointer}.diagram-source pre{color:var(--ink)}
footer{font-size:13px;color:var(--muted);margin-top:46px;padding-top:22px;border-top:1px solid var(--line)}
@media(max-width:750px){body{font-size:15px}.sheet{margin:0;padding:26px 18px 40px;border-radius:0;box-shadow:none}h1{font-size:25px}h2{font-size:21px}h3{font-size:18px}.toc ul{columns:1}table{min-width:620px}img{margin:18px auto}pre{padding:12px}.reading-links a{display:inline-block}.eyebrow{letter-spacing:1px}}
@media print{body{background:#fff;font-size:11pt}.sheet{margin:0;padding:0;max-width:none;box-shadow:none}.toc,.reading-links,.diagram-source{display:none}h2,h3{break-after:avoid}tr,img,pre,blockquote{break-inside:avoid}a{color:inherit;text-decoration:none}table{font-size:9pt;min-width:0}.table-wrap{overflow:visible}h2{font-size:17pt}h3{font-size:13pt}}
</style></head><body><main class="sheet">
<p class="eyebrow">FENG PAN · UCB/EECS-2013-38 · RESEARCH NOTES</p>
<p class="reading-links"><a href="${basename}.md">查看 Markdown 笔记</a><a href="source/EECS-2013-38.pdf#page=55">打开原论文建模章节</a><span>点击图片可放大</span></p>
<details class="toc"><summary>内容导航：格点、事件、框图、形貌与论文研究</summary><ul>${nav}</ul></details>
${body}
<footer>基于本文件夹笔记与伯克利原论文核对整理。示意图明确区分原文参数与复现建议；未运行器件仿真。</footer>
</main></body></html>`;
await fs.writeFile(path.join(root, basename + '.html'), html, 'utf8');

const refs = [...md.matchAll(/!?\[[^\]]*\]\(([^)]+)\)/g)].map(m=>m[1]).filter(x=>!/^https?:/.test(x));
const missing=[];
for(const ref of refs){try{await fs.access(path.resolve(root,decodeURIComponent(ref.split('#')[0])));}catch{missing.push(ref);}}
if(missing.length)throw new Error('失效本地引用：'+missing.join(', '));
console.log(JSON.stringify({html:path.join(root,basename+'.html'),markdownCharacters:md.length,localReferences:refs.length,sections:toc.length,missingReferences:missing},null,2));
