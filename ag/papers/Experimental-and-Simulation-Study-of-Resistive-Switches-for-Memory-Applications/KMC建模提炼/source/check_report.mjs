import fs from 'node:fs/promises';
import path from 'node:path';
import { createRequire } from 'node:module';
import { fileURLToPath, pathToFileURL } from 'node:url';
const require = createRequire('C:/Users/31826/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright/package.json');
const { chromium } = require('playwright');
const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const qa = path.join(root, 'source', 'qa');
await fs.mkdir(qa, {recursive: true});
const browser = await chromium.launch({headless:true, channel:'msedge'});
const page = await browser.newPage({viewport:{width:1280,height:1000},deviceScaleFactor:1});
const errors=[];
page.on('pageerror',e=>errors.push(e.message));
await page.goto(pathToFileURL(path.join(root,'KMC搭建过程与论文研究拆解.html')).href);
await page.evaluate(()=>document.fonts.ready);
const desktop = await page.evaluate(()=>({
  pageWidth:innerWidth, contentWidth:document.documentElement.scrollWidth,
  headings:document.querySelectorAll('h2').length,
  missingImages:[...document.images].filter(i=>!i.complete||!i.naturalWidth).map(i=>i.src),
  imageCount:document.images.length
}));
await page.screenshot({path:path.join(qa,'report-desktop.png')});
await page.locator('img[src="assets/lattice-schematic.svg"]').screenshot({path:path.join(qa,'lattice.png')});
await page.locator('img[src="assets/event-state-machine.svg"]').screenshot({path:path.join(qa,'events.png')});
await page.locator('img[src="assets/kmc-loop.svg"]').screenshot({path:path.join(qa,'loop.png')});
await page.locator('img[src="assets/figure-2-20.png"]').screenshot({path:path.join(qa,'morphology.png')});
await page.setViewportSize({width:390,height:844});
await page.goto(pathToFileURL(path.join(root,'KMC搭建过程与论文研究拆解.html')).href);
await page.evaluate(()=>document.fonts.ready);
const mobile=await page.evaluate(()=>({pageWidth:innerWidth,contentWidth:document.documentElement.scrollWidth}));
await page.screenshot({path:path.join(qa,'report-mobile.png')});
console.log(JSON.stringify({desktop,mobile,errors},null,2));
await browser.close();
if(desktop.missingImages.length||errors.length||mobile.contentWidth>mobile.pageWidth)process.exitCode=1;
