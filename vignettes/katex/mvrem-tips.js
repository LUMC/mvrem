// Popups for the estimates in mvrem's LaTeX equations (latexEquation(m, annotate = "tooltip")).
// KaTeX turns \htmlData{par=...,mean=...}{symbol} into <span data-par=... data-mean=...> (render
// option trust: true). After the equations are rendered, this script writes a formatted text into
// a data-tip attribute of every such span; the popup itself is pure CSS (:hover::after), so no
// mouse-event handling is involved. Vectors are ";"-separated (several parameters of one symbol,
// e.g. a covariance matrix). mvremTips.annotate() can be called again after re-rendering.
(function () {
  // KaTeX builds an aligned environment as a table: a cell with the rows (positioned zero-height
  // blocks overflowing each other) and, below it, an empty positioned "depth" cell as tall as
  // everything under the baseline, i.e. the bottom half of the formula. Those transparent boxes
  // paint above the symbols and swallow the pointer. So nothing in a rendered formula receives
  // the pointer except the annotated spans.
  var CSS = ".katex .katex-html, .katex .katex-html *{pointer-events:none}" +
    ".katex [data-par], .katex [data-par] *{pointer-events:auto}" +
    ".katex [data-par]{position:relative;border-bottom:1px dotted #d8d8d8;cursor:help}" +
    ".katex [data-par]:hover{color:#c0392b}" +
    ".katex [data-tip]:hover::after{content:attr(data-tip);position:absolute;left:1.6em;top:2.3em;" +
    "z-index:10000;pointer-events:none;background:rgba(30,30,30,.93);color:#fff;" +
    "font:12px/1.35 sans-serif;font-style:normal;font-weight:normal;text-align:left;" +
    "padding:4px 7px;border-radius:3px;white-space:pre-wrap;width:max-content;max-width:26em}";
  function split(v) { return v === null || v === undefined ? null : String(v).split(";"); }
  function text(el) {
    var d = el.dataset, par = d.par || "", item = split(d.item), mean = split(d.mean), sd = split(d.sd),
        lo = split(d.lo), hi = split(d.hi), band = d.band || "95%", lines = [par];
    var n = Math.max(mean ? mean.length : 0, item ? item.length : 0);
    for (var i = 0; i < n; i++) {
      var l = item && item[i] ? item[i] + ": " : "";
      if (mean) l += "mean " + mean[i] + (sd && sd[i] ? " (sd " + sd[i] + ")" : "");
      if (lo && hi && lo[i] !== undefined) l += "  " + band + " [" + lo[i] + ", " + hi[i] + "]";
      if (l) lines.push(l);
    }
    if (d.note) lines.push(d.note);
    return lines.join("\n");
  }
  function annotate(root) {
    var els = (root || document).querySelectorAll("[data-par]");
    for (var i = 0; i < els.length; i++) els[i].setAttribute("data-tip", text(els[i]));
    return els.length;
  }
  if (typeof window !== "undefined") window.mvremTips = { text: text, annotate: annotate };
  if (typeof document !== "undefined") {
    document.addEventListener("DOMContentLoaded", function () {
      var st = document.createElement("style"); st.textContent = CSS;
      document.head.appendChild(st);
      setTimeout(annotate, 0);                     // after the KaTeX render of the other handlers
    });
    window.addEventListener("load", function () { annotate(); });
  }
})();
