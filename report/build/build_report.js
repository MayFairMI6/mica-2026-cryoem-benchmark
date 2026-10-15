// Build the Word report from report/REPORT_DRAFT.md (+ report/RESULTS_TABLES.md as an appendix of full tables).
// Usage: node build_report.js   -> ../MICA_benchmark_report.docx
// The Markdown files stay the single source of text and numbers; re-run after they change.
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell, WidthType, ShadingType,
  AlignmentType, LevelFormat, PageOrientation, BorderStyle, Footer, PageNumber, ImageRun,
} = require("docx");

const ROOT = path.resolve(__dirname, "..");
const FONT = "Calibri";

function inline(text, base = {}) {
  // **bold**, *italic*, `code`
  const runs = [];
  const re = /(\*\*[^*]+\*\*|`[^`]+`|\*[^*\s][^*]*\*)/g;
  let last = 0, m;
  while ((m = re.exec(text)) !== null) {
    if (m.index > last) runs.push(new TextRun({ text: text.slice(last, m.index), font: FONT, ...base }));
    const t = m[0];
    if (t.startsWith("**")) runs.push(new TextRun({ text: t.slice(2, -2), bold: true, font: FONT, ...base }));
    else if (t.startsWith("`")) runs.push(new TextRun({ text: t.slice(1, -1), font: "Consolas", ...base, size: (base.size ? base.size : 20) }));
    else runs.push(new TextRun({ text: t.slice(1, -1), italics: true, font: FONT, ...base }));
    last = m.index + t.length;
  }
  if (last < text.length) runs.push(new TextRun({ text: text.slice(last), font: FONT, ...base }));
  return runs;
}

function table(rows, contentWidth) {
  const cells = rows.map((r) => r.replace(/^\s*\|/, "").replace(/\|\s*$/, "").split("|").map((c) => c.trim()));
  const header = cells[0], body = cells.slice(2);
  const n = header.length;
  const plain = (s) => s.replace(/\*\*|`/g, "");
  const want = header.map((_, i) => Math.max(...[header, ...body].map((r) => Math.min(plain(r[i] || "").length, 28)), 4) + 2);
  const tot = want.reduce((a, b) => a + b, 0);
  const MIN = n > 8 ? 760 : 900;                       // keep short codes such as "9SYT" on one line
  let widths = want.map((w) => Math.max(Math.floor((w / tot) * contentWidth), MIN));
  const over = widths.reduce((a, b) => a + b, 0) - contentWidth;
  if (over > 0) { const flex = widths.reduce((a, w) => a + Math.max(w - MIN, 0), 0); widths = widths.map((w) => w - Math.floor((Math.max(w - MIN, 0) / flex) * over)); }
  widths[n - 1] += contentWidth - widths.reduce((a, b) => a + b, 0);
  const size = n > 6 ? 17 : 18;                        // 8.5 pt / 9 pt (QC protocol section 8)
  const border = { style: BorderStyle.SINGLE, size: 4, color: "999999" };
  // short tables are not split; long ones may split, but never within their first six or last four rows (no stranded title or one-row tail)
  const keepRow = (k) => k < body.length - 1 && (body.length <= 14 || k < 6 || k >= body.length - 4);
  const mk = (r, isHead, k) => new TableRow({
    tableHeader: isHead, cantSplit: true,
    children: header.map((_, i) => new TableCell({
      width: { size: widths[i], type: WidthType.DXA },
      margins: { top: body.length > 20 ? 16 : 40, bottom: body.length > 20 ? 16 : 40, left: 70, right: 70 },
      shading: isHead ? { type: ShadingType.CLEAR, color: "auto", fill: "E7ECF3" } : undefined,
      borders: { top: border, bottom: border, left: border, right: border },
      children: [new Paragraph({ keepNext: (isHead || keepRow(k)) || undefined, spacing: { before: 0, after: 0 }, children: inline(r[i] || "", { size, bold: isHead || undefined }) })],
    })),
  });
  return new Table({ width: { size: contentWidth, type: WidthType.DXA }, columnWidths: widths, rows: [mk(header, true, -1), ...body.map((r, k) => mk(r, false, k))] });
}

function convert(md, contentWidth, opts = {}) {
  const out = [];
  const lines = md.split("\n");
  let i = 0, para = [], inst = opts.instBase || 0, inNum = false;
  const flush = (beforeTable) => {
    if (!para.length) return;
    const text = para.join(" ").trim(); para = [];
    if (!text) return;
    const note = /^\*\[.*\]\*$/.test(text) || /^\*[^*].*\*$/.test(text);
    const caption = /^\*?Figure \d/.test(text);
    out.push(new Paragraph({ alignment: caption ? AlignmentType.CENTER : undefined, keepNext: beforeTable || undefined, spacing: { after: caption ? 200 : 120, line: 276 }, children: note ? inline(text.replace(/^\*|\*$/g, ""), { italics: true, color: "555555", size: 20 }) : inline(text, { size: 22 }) }));
  };
  while (i < lines.length) {
    const l = lines[i];
    if (/^\s*$/.test(l)) { let j = i; while (j < lines.length && /^\s*$/.test(lines[j])) j++; flush(j < lines.length && /^\s*\|/.test(lines[j])); i++; continue; }   // a line that introduces a table stays with it
    const im = l.match(/^!\[[^\]]*\]\(([^)]+)\)\s*$/);
    if (im) {   // figure: PNG placed at full text width, aspect ratio from the file header
      flush(); const buf = fs.readFileSync(path.join(ROOT, im[1])); const w = buf.readUInt32BE(16), h = buf.readUInt32BE(20); const W = 620;
      out.push(new Paragraph({ alignment: AlignmentType.CENTER, keepNext: true, spacing: { before: 120, after: 60 }, children: [new ImageRun({ type: "png", data: buf, transformation: { width: W, height: Math.round(W * h / w) } })] }));
      i++; continue;
    }
    const h = l.match(/^(#{1,3})\s+(.*)$/);
    if (h) {
      flush(); inNum = false;
      const lvl = h[1].length + (opts.demote || 0);
      if (lvl === 1) out.push(new Paragraph({ alignment: AlignmentType.LEFT, spacing: { after: 200 }, children: [new TextRun({ text: h[2], bold: true, size: 34, font: FONT })] }));
      else {
        // appendix: a long table starts on its own page, so that it is not split (QC protocol section 8)
        let j = i + 1; while (j < lines.length && /^\s*$/.test(lines[j])) j++; let nrows = 0; while (j + nrows < lines.length && /^\s*\|/.test(lines[j + nrows])) nrows++;
        const brk = !!opts.breakLong && nrows - 2 > 14 && out.length > 3;
        out.push(new Paragraph({ heading: lvl === 2 ? HeadingLevel.HEADING_1 : HeadingLevel.HEADING_2, keepNext: true, pageBreakBefore: brk || undefined, children: inline(h[2], { size: lvl === 2 ? 28 : 23, bold: true }) }));
      }
      i++; continue;
    }
    if (/^\s*\|/.test(l)) {
      flush(); const rows = [];
      while (i < lines.length && /^\s*\|/.test(lines[i])) rows.push(lines[i++]);
      out.push(table(rows, contentWidth)); out.push(new Paragraph({ spacing: { after: 120 }, children: [] }));
      continue;
    }
    const li = l.match(/^(\s*)(-|\d+\.)\s+(.*)$/);
    if (li) {
      flush();
      let text = li[3]; const numbered = /\d+\./.test(li[2]); const level = li[1].length >= 2 ? 1 : 0; i++;
      if (numbered && level === 0 && !inNum) { inst++; inNum = true; }   // each numbered list restarts at 1
      while (i < lines.length && /^\s{2,}\S/.test(lines[i]) && !/^\s*(-|\d+\.)\s+/.test(lines[i]) && !/^\s*\|/.test(lines[i])) text += " " + lines[i++].trim();
      out.push(new Paragraph({ numbering: { reference: numbered ? "nums" : "bullets", level, instance: numbered ? inst : 0 }, spacing: { after: 70, line: 264 }, children: inline(text, { size: 22 }) }));
      continue;
    }
    inNum = false; para.push(l.trim()); i++;
  }
  flush();
  return out;
}

const report = fs.readFileSync(path.join(ROOT, "REPORT_DRAFT.md"), "utf8");
const tables = fs.readFileSync(path.join(ROOT, "RESULTS_TABLES.md"), "utf8")
  .replace(/^# .*\n/, "## Appendix B — Reported result tables\n")
  .replace(/^### /gm, "### ");

const footer = new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ children: ["Page ", PageNumber.CURRENT], size: 16, color: "777777", font: FONT })] })] });
const numbering = { config: [
  { reference: "bullets", levels: [
    { level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 500, hanging: 250 } } } },
    { level: 1, format: LevelFormat.BULLET, text: "–", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 900, hanging: 250 } } } }] },
  { reference: "nums", levels: [
    { level: 0, format: LevelFormat.DECIMAL, text: "%1.", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 500, hanging: 300 } } } },
    { level: 1, format: LevelFormat.LOWER_LETTER, text: "%2.", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 900, hanging: 300 } } } }] },
] };

const doc = new Document({
  creator: "MICA benchmark", title: "Benchmarking MICA on 2026 EMDB/PDB data",
  numbering,
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 28, bold: true, font: FONT, color: "1F3A5F" }, paragraph: { spacing: { before: 300, after: 120 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 23, bold: true, font: FONT, color: "1F3A5F" }, paragraph: { spacing: { before: 220, after: 90 }, outlineLevel: 1 } },
    ],
  },
  sections: [
    { properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1200, bottom: 1200, left: 1150, right: 1150 } } },
      footers: { default: footer }, children: convert(report, 12240 - 2 * 1150) },
    { properties: { page: { size: { width: 12240, height: 15840, orientation: PageOrientation.LANDSCAPE }, margin: { top: 900, bottom: 900, left: 800, right: 800 } } },
      footers: { default: footer }, children: convert(tables, 15840 - 2 * 800, { instBase: 100, breakLong: true }) },
  ],
});

Packer.toBuffer(doc).then((buf) => {
  const out = path.join(ROOT, "MICA_benchmark_report.docx");
  fs.writeFileSync(out, buf); console.log("wrote", out, buf.length, "bytes");
});
