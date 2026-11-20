// mvrem interactive plots: panels of subject lines / intervals / points / curves drawn as
// SVG, with a tooltip and a highlight of the subject under the mouse that is linked across
// all panels. Self-contained: no external libraries.
(function () {
  var NS = "http://www.w3.org/2000/svg", W = 440, H = 300;
  // tooltip word per line kind ("S01 observed"); other kinds show the id only
  var KIND = { obs: "observed", fit: "fitted", res: "residual", pred: "predicted", predInt: "95% bound",
               fitMissing: "predicted (missing observation)" };
  // default colours of elements without a `color` in the spec: the widget tones of plotColors
  // in R/plotUtils.R (subject lines, fitted lines, accent) and the interval lines
  var COL = { subject: "#c9c9c9", accentLight: "#9cc3e3", accent: "#1b6ca8", interval: "#9a9a9a" };

  var CSS = [
    ".mvp{position:relative;font:12px/1.3 sans-serif;color:#222;background:#fff}",
    ".mvp .grid{display:grid;gap:4px}",
    ".mvp .grid svg{width:100%;height:auto;display:block}",
    ".mvp .lg svg{width:24px;height:10px;flex:none}",
    ".mvp .bar{display:flex;flex-wrap:wrap;gap:6px 28px;align-items:center;padding:2px 4px}",
    // legend entry: swatch immediately left of its label, entries far apart
    ".mvp .lg{cursor:pointer;user-select:none;display:inline-flex;align-items:center;gap:3px;white-space:nowrap}",
    ".mvp .lg.off{opacity:.35;text-decoration:line-through}",
    ".mvp input{font:inherit;padding:1px 5px;width:14em;border:1px solid #bbb;border-radius:3px}",
    ".mvp select{font:inherit;padding:0 2px;border:1px solid #bbb;border-radius:3px}",
    ".mvp .cby{display:inline-flex;align-items:center;gap:4px;white-space:nowrap}",
    ".mvp .clev{display:inline-flex;flex-wrap:wrap;gap:2px 10px}",
    ".mvp .clev span{display:inline-flex;align-items:center;gap:3px}",
    ".mvp .clev i{display:inline-block;width:9px;height:9px;border-radius:50%}",
    ".mvp .tip{position:absolute;pointer-events:none;background:rgba(30,30,30,.92);color:#fff;",
    "padding:3px 6px;border-radius:3px;white-space:pre;display:none;font-size:11px;z-index:5}",
    ".mvp g.el{cursor:pointer}",
    ".mvp .hit{fill:none;stroke:transparent;stroke-width:9;pointer-events:stroke}",
    ".mvp circle.hit{fill:transparent;stroke:none;pointer-events:all}",
    ".mvp.dim g.el:not(.hl){opacity:.12}",
    // highlighted subject: a softer red for the lines, the full red for its observed points
    ".mvp g.el.hl path.vis,.mvp g.el.hl line.vis{stroke:#efb8b0 !important;stroke-width:2 !important}",
    ".mvp g.el.hl circle.vis{fill:#c0392b !important}",
    // (weighted points, class wt, keep their size: their area shows the weight)
    ".mvp g.el.hl circle.vis:not(.wt){r:3}",
    // open circles (predictions of missing observations) stay open when highlighted
    ".mvp g.el.hl circle.vis.hollow{fill:none !important;stroke:#c0392b !important;r:3.5}"
  ].join("\n");

  // Legend toggles: a click on the entry of kind K sets class off-K on the widget's root, and
  // the rule ".mvp.off-K .kind-K{display:none}" hides the elements of that kind. The rules are
  // made from the kinds of each spec's legend (kindRules()), so a kind defined in R needs no
  // change here. They go to their own <style> element, shared by all widgets of the page.
  var KIND_CSS_ID = "mvremPlotKindCss", kindsWithRule = {};
  function kindRules(legend) {
    var st = document.getElementById(KIND_CSS_ID), css = "";
    if (!st) {
      kindsWithRule = {};                // the page lost the element (e.g. re-rendered): start again
      st = document.createElement("style"); st.id = KIND_CSS_ID;
      document.head.appendChild(st);
    }
    legend.forEach(function (e) {
      var k = String(e.kind);
      // kinds are class-name fragments: skip anything that is not a plain identifier
      if (!/^[A-Za-z_][\w-]*$/.test(k) || kindsWithRule[k]) return;
      kindsWithRule[k] = true;
      css += ".mvp.off-" + k + " .kind-" + k + "{display:none}\n";
    });
    if (css) st.textContent += css;
  }

  function arr(v) { return v === null || v === undefined ? [] : [].concat(v); }
  function S(tag, attrs, parent) {
    var e = document.createElementNS(NS, tag);
    for (var k in attrs) if (attrs.hasOwnProperty(k)) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function H_(tag, cls, parent, text) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text !== undefined) e.textContent = text;
    if (parent) parent.appendChild(e);
    return e;
  }
  function fmt(v) { return String(+(+v).toPrecision(4)); }
  function ticks(lo, hi, n) {
    var raw = (hi - lo) / n, mag = Math.pow(10, Math.floor(Math.log(raw) / Math.LN10)), f = raw / mag;
    var step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag, out = [];
    // ticks as integer multiples of the step, rounded to the step's decimals: no accumulated
    // floating-point error (which showed zero as 5.6e-17)
    var dec = Math.max(0, -Math.floor(Math.log(step) / Math.LN10));
    for (var k = Math.ceil(lo / step - 1e-9); k * step <= hi + 1e-9; k++) out.push(+(k * step).toFixed(dec));
    return out;
  }
  function range(vals) {
    var lo = Infinity, hi = -Infinity;
    vals.forEach(function (v) { if (isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; } });
    return isFinite(lo) ? [lo, hi] : [0, 1];
  }
  function pad(r, f) {
    var d = r[1] - r[0];
    if (d === 0) return [r[0] - 1, r[1] + 1];
    return [r[0] - d * f, r[1] + d * f];
  }

  var uid = 0;

  // One panel: returns its SVG; registers highlightable groups with ctx.add(key, g)
  function drawPanel(p, pi, ctx) {
    var svg = S("svg", { viewBox: "0 0 " + W + " " + H });
    var ylabels = p.ylabels ? { pos: arr(p.ylabels.pos), text: arr(p.ylabels.text) } : null;
    var left = ylabels ? Math.min(200, 14 + 5.6 * Math.max.apply(null, ylabels.text.map(function (t) { return t.length; }))) : 52;
    var L = left, R = W - 10, T = 26, B = H - 40;
    S("text", { x: (L + R) / 2, y: 15, "text-anchor": "middle", "font-weight": "bold", "font-size": 12 }, svg).textContent = p.title || "";
    if (p.empty) return svg;

    var lines = arr(p.lines), curves = arr(p.curves), iv = p.intervals, pts = p.points;
    var xs = [], ys = [];
    lines.forEach(function (l) { xs = xs.concat(arr(l.x)); ys = ys.concat(arr(l.y)); });
    curves.forEach(function (c) { xs = xs.concat(arr(c.x)); ys = ys.concat(arr(c.y)); });
    if (iv) { xs = xs.concat(arr(iv.lo), arr(iv.hi)); ys = ys.concat(arr(iv.pos)); }
    if (pts) { xs = xs.concat(arr(pts.x)); ys = ys.concat(arr(pts.y)); }
    var xr = p.xlim ? arr(p.xlim) : pad(range(xs), 0.04), yr = p.ylim ? arr(p.ylim) : pad(range(ys), 0.04);
    if (p.xlim) xr = pad(xr, 0.04);
    if (p.ylim) yr = pad(yr, 0.04);
    if (iv) yr = [0.5, arr(iv.pos).length + 0.5];
    var rev = !!p.yReverse;
    function sx(v) { return L + (v - xr[0]) / (xr[1] - xr[0]) * (R - L); }
    function sy(v) { var f = (v - yr[0]) / (yr[1] - yr[0]); return rev ? T + f * (B - T) : B - f * (B - T); }
    function ix(px) { return xr[0] + (px - L) / (R - L) * (xr[1] - xr[0]); }

    // frame, grid, ticks, labels
    ticks(xr[0], xr[1], 6).forEach(function (v) {
      if (v < xr[0] || v > xr[1]) return;
      S("line", { x1: sx(v), x2: sx(v), y1: T, y2: B, stroke: "#ececec" }, svg);
      S("text", { x: sx(v), y: B + 13, "text-anchor": "middle", "font-size": 10 }, svg).textContent = fmt(v);
    });
    if (ylabels) {
      ylabels.pos.forEach(function (v, i) {
        S("text", { x: L - 5, y: sy(v) + 3, "text-anchor": "end", "font-size": 10 }, svg).textContent = ylabels.text[i];
      });
    } else {
      ticks(yr[0], yr[1], 5).forEach(function (v) {
        if (v < yr[0] || v > yr[1]) return;
        S("line", { x1: L, x2: R, y1: sy(v), y2: sy(v), stroke: "#ececec" }, svg);
        S("text", { x: L - 5, y: sy(v) + 3, "text-anchor": "end", "font-size": 10 }, svg).textContent = fmt(v);
      });
    }
    if (p.vline !== undefined && p.vline >= xr[0] && p.vline <= xr[1])
      S("line", { x1: sx(p.vline), x2: sx(p.vline), y1: T, y2: B, stroke: "#999", "stroke-dasharray": "2 3" }, svg);
    if (p.hline !== undefined && p.hline >= yr[0] && p.hline <= yr[1])
      S("line", { x1: L, x2: R, y1: sy(p.hline), y2: sy(p.hline), stroke: "#999", "stroke-dasharray": "2 3" }, svg);
    S("rect", { x: L, y: T, width: R - L, height: B - T, fill: "none", stroke: "#555" }, svg);
    S("text", { x: (L + R) / 2, y: H - 5, "text-anchor": "middle", "font-size": 11 }, svg).textContent = p.xlab || "";
    if (p.ylab) S("text", { transform: "translate(11," + (T + B) / 2 + ") rotate(-90)", "text-anchor": "middle", "font-size": 11 }, svg).textContent = p.ylab;
    // a note inside the plot area (top left), e.g. the estimated correlation of a pairs panel
    if (p.note) S("text", { x: L + 6, y: T + 13, "font-size": 10.5, fill: "#444" }, svg).textContent = p.note;

    var clip = "mvc" + (++uid);
    S("rect", { x: L, y: T, width: R - L, height: B - T }, S("clipPath", { id: clip }, S("defs", {}, svg)));
    var plot = S("g", { "clip-path": "url(#" + clip + ")" }, svg);
    function path(xv, yv) {
      var d = "";
      for (var i = 0; i < xv.length; i++) d += (i ? "L" : "M") + sx(xv[i]).toFixed(1) + "," + sy(yv[i]).toFixed(1);
      return d;
    }
    function group(id, kind, extra) {
      var g = S("g", { "class": "el kind-" + kind }, plot);
      g._id = id; g._info = extra; g._svg = svg; g._ix = ix; g._panel = p;
      g._key = ctx.linked ? String(id) : pi + "|" + id;
      ctx.add(g._key, g);
      return g;
    }

    // subject lines (observed / residual: solid with points; fitted: dotted). Only "fit" has
    // the fitted defaults; any other kind is drawn like observed values unless the line brings
    // its own color / dash / width / points (as the predictions do).
    lines.forEach(function (l) {
      var xv = arr(l.x), yv = arr(l.y), obs = l.kind !== "fit";
      if (l.hollow) {                    // open circles, never joined (predicted missing observations)
        var gh = group(l.id, l.kind, { x: xv, y: yv, label: String(l.id) + " " + (KIND[l.kind] || ""), p: p });
        xv.forEach(function (v, i) {
          S("circle", { cx: sx(v), cy: sy(yv[i]), r: 3, fill: "none", stroke: l.color || COL.accent,
                        "stroke-width": 1.2, "class": "vis hollow" }, gh);
          S("circle", { cx: sx(v), cy: sy(yv[i]), r: 6, "class": "hit" }, gh);
        });
        return;
      }
      var col = l.color || (obs ? COL.subject : COL.accentLight), dash = l.dash !== undefined ? l.dash : (obs ? "" : "2 3");
      var g = group(l.id, l.kind || "obs", { x: xv, y: yv, label: String(l.id) + " " + (KIND[l.kind] || ""), p: p });
      if (l.cv) { g._cv = l.cv; g._col = col; g._pts = []; }    // colorBy: values per point
      if (p.connect !== false && xv.length > 1) {
        var attrs = { d: path(xv, yv), fill: "none", stroke: col, "stroke-width": l.width || (obs ? 1 : 1.3), "class": "vis" };
        if (dash) attrs["stroke-dasharray"] = dash;                  // fitted: dotted, observed: solid
        var pv = S("path", attrs, g);
        if (l.cv) g._line = pv;
      }
      if ((obs && l.points !== false) || p.connect === false)
        xv.forEach(function (v, i) {
          // `size`: radius factor per point (weights(): area proportional to the weight)
          var f = l.size ? +arr(l.size)[i] : 1;
          var c = S("circle", { cx: sx(v), cy: sy(yv[i]), r: 1.8 * (isFinite(f) ? f : 1), fill: col,
                                "class": l.size ? "vis wt" : "vis" }, g);
          if (l.cv) g._pts.push(c);
        });
      if (xv.length > 1 && p.connect !== false) S("path", { d: path(xv, yv), "class": "hit" }, g);
      else xv.forEach(function (v, i) { S("circle", { cx: sx(v), cy: sy(yv[i]), r: 6, "class": "hit" }, g); });
    });

    // intervals (caterpillar / forest)
    if (iv) {
      var id = arr(iv.id), pos = arr(iv.pos), mean = arr(iv.mean), lo = arr(iv.lo), hi = arr(iv.hi), tip = arr(iv.tip);
      id.forEach(function (v, i) {
        var g = group(v, "int", { tip: tip[i] });
        S("line", { x1: sx(lo[i]), x2: sx(hi[i]), y1: sy(pos[i]), y2: sy(pos[i]), stroke: COL.interval, "stroke-width": 1.6, "class": "vis" }, g);
        S("circle", { cx: sx(mean[i]), cy: sy(pos[i]), r: 2.2, fill: COL.accent, "class": "vis" }, g);
        S("line", { x1: sx(lo[i]), x2: sx(hi[i]), y1: sy(pos[i]), y2: sy(pos[i]), "class": "hit" }, g);
      });
    }

    // points (scatter)
    if (pts) {
      var pid = arr(pts.id), px = arr(pts.x), py = arr(pts.y), ptip = arr(pts.tip);
      pid.forEach(function (v, i) {
        var g = group(v, "pt", { tip: ptip[i] });
        S("circle", { cx: sx(px[i]), cy: sy(py[i]), r: 2.6, fill: COL.accent, "fill-opacity": 0.75, "class": "vis" }, g);
        S("circle", { cx: sx(px[i]), cy: sy(py[i]), r: 7, "class": "hit" }, g);
      });
    }

    // mean curves / densities
    curves.forEach(function (c) {
      var xv = arr(c.x), yv = arr(c.y);
      var attrs = { d: path(xv, yv), fill: "none", stroke: c.color || COL.accent, "stroke-width": c.width || 2,
                    "class": "kind-" + (c.kind || "curve") };
      if (c.dash) attrs["stroke-dasharray"] = c.dash;
      if (c.id !== undefined) {
        var g = group(c.id, c.kind || "curve", { x: xv, y: yv, label: c.label || String(c.id), p: p });
        attrs["class"] = "vis"; S("path", attrs, g);
        S("path", { d: attrs.d, "class": "hit" }, g);
      } else {
        attrs["pointer-events"] = "none"; S("path", attrs, plot);
      }
    });
    return svg;
  }

  function render(el, spec) {
    el.innerHTML = "";
    el.style.height = "auto";          // htmlwidgets sets a fixed height: let the panels decide
    el.style.overflow = "visible";
    if (!document.getElementById("mvremPlotCss")) {
      var st = document.createElement("style"); st.id = "mvremPlotCss"; st.textContent = CSS;
      document.head.appendChild(st);
    }
    kindRules(arr(spec.legend));
    var root = H_("div", "mvp", el), byId = {}, pinned = {}, match = {}, hover = null, ids = [];
    var ctx = {
      linked: spec.linked !== false,
      add: function (key, g) { (byId[key] = byId[key] || []).push(g); if (ids.indexOf(key) < 0) ids.push(key); }
    };

    function refresh() {
      var any = false;
      for (var k in byId) {
        var on = k === hover || pinned[k] || match[k];
        byId[k].forEach(function (g) { g.classList.toggle("hl", !!on); });
        if (on) any = true;
      }
      root.classList.toggle("dim", any);
    }

    // bar: search box + legend
    var bar = H_("div", "bar", root), count = null;
    if (spec.search) {
      var inp = H_("input", "", bar); inp.type = "search"; inp.placeholder = "highlight " + spec.search;
      count = H_("span", "", bar);
      inp.addEventListener("input", function () {
        var q = inp.value.trim().toLowerCase(), n = 0;
        match = {};
        if (q) ids.forEach(function (k) {
          var label = k.indexOf("|") >= 0 && !ctx.linked ? k.split("|").slice(1).join("|") : k;
          if (label.toLowerCase().indexOf(q) >= 0) { match[k] = true; n++; }
        });
        count.textContent = q ? n + " match" + (n === 1 ? "" : "es") : "";
        refresh();
      });
    }
    arr(spec.legend).forEach(function (e) {
      var item = H_("span", "lg", bar);
      var sw = S("svg", { width: 24, height: 10 }, item);
      if (e.shape === "circle") {
        S("circle", { cx: 12, cy: 5, r: 3.5, fill: "none", stroke: e.color, "stroke-width": 1.3 }, sw);
      } else {
        var ln = S("line", { x1: 1, x2: 23, y1: 5, y2: 5, stroke: e.color, "stroke-width": 2.5 }, sw);
        if (e.dash) ln.setAttribute("stroke-dasharray", e.dash);
      }
      item.appendChild(document.createTextNode(e.label));
      item.addEventListener("click", function () {
        var off = !root.classList.contains("off-" + e.kind);
        root.classList.toggle("off-" + e.kind, off); item.classList.toggle("off", off);
      });
    });

    // colorBy: dropdown "color: gray / <variable>" and the levels of the chosen variable
    var cb = spec.colorBy, colorTo = null;
    if (cb && arr(cb.vars).length) {
      var box = H_("span", "cby", bar);
      H_("span", "", box, "color:");
      var sel = H_("select", "", box), levBox = H_("span", "clev", bar);
      [""].concat(arr(cb.vars)).forEach(function (v) {
        var o = H_("option", "", sel, v === "" ? "gray (default)" : v); o.value = v;
      });
      var pal = arr(cb.palette);
      var colorFor = function (v, x) {
        var k = arr(cb.levels[v]).map(String).indexOf(String(x));
        return k < 0 || x === null || x === undefined ? null : pal[k % pal.length];
      };
      colorTo = function (v) {
        sel.value = v || "";
        levBox.innerHTML = "";
        if (v) arr(cb.levels[v]).forEach(function (lv) {
          var it = H_("span", "", levBox); var dot = H_("i", "", it); dot.style.background = colorFor(v, lv);
          it.appendChild(document.createTextNode(String(lv)));
        });
        for (var k in byId) byId[k].forEach(function (g) {
          if (!g._cv) return;
          // colours from the values (lines without points, e.g. prediction bounds, too)
          var vals = v && g._cv[v] !== undefined ? arr(g._cv[v]) : null;
          var cols = vals ? vals.map(function (x) { return colorFor(v, x); }) : [];
          g._pts.forEach(function (c, i) { c.setAttribute("fill", cols[i] || g._col); });
          if (g._line) {
            var one = vals && cols.length && cols.every(function (c) { return c && c === cols[0]; });
            g._line.setAttribute("stroke", one ? cols[0] : g._col);
            g._line.setAttribute("stroke-opacity", one ? 0.7 : 1);
          }
        });
      };
      sel.addEventListener("change", function () { colorTo(sel.value); });
    }

    var grid = H_("div", "grid", root);
    grid.style.gridTemplateColumns = "repeat(" + (spec.ncol || 1) + ",minmax(0,1fr))";
    arr(spec.panels).forEach(function (p, i) { grid.appendChild(drawPanel(p, i, ctx)); });
    if (colorTo && typeof cb.selected === "string" && cb.selected) colorTo(cb.selected);

    var tip = H_("div", "tip", root);
    // tooltip: the panel title first, then the element and its values
    function tipText(g) {
      var d = g._info || {}, title = (g._panel && g._panel.title) ? g._panel.title + "\n" : "";
      if (d.tip) return title + d.tip;
      if (!d.x) return title + String(g._id);
      return title + d.label;
    }
    function nearest(g, e) {
      var d = g._info, r = g._svg.getBoundingClientRect();
      if (!d || !d.x || !r.width) return tipText(g);
      var xv = arr(d.x), yv = arr(d.y), mx = g._ix((e.clientX - r.left) * W / r.width), best = 0;
      for (var i = 1; i < xv.length; i++) if (Math.abs(xv[i] - mx) < Math.abs(xv[best] - mx)) best = i;
      var p = d.p || g._panel;
      return (p.title ? p.title + "\n" : "") + d.label + "\n" +
        (p.xlab || "x") + " = " + fmt(xv[best]) + ", " + (p.ylab || "y") + " = " + fmt(yv[best]);
    }
    function target(e) { var t = e.target; return t && t.closest ? t.closest("g.el") : null; }

    // The hovered subject follows the mouse position only (no mouseover / mouseout pairs,
    // whose order between the overlapping hit areas left highlights behind); leaving the
    // widget clears it.
    function setHover(k) { if (k !== hover) { hover = k; refresh(); } }
    root.addEventListener("mouseleave", function () { setHover(null); tip.style.display = "none"; });
    root.addEventListener("mousemove", function (e) {
      var g = target(e);
      if (!g) { setHover(null); tip.style.display = "none"; return; }
      setHover(g._key);
      tip.textContent = nearest(g, e) + (pinned[g._key] ? "\n(pinned: click to unpin)" : "\n(click to pin)");
      var r = root.getBoundingClientRect();
      tip.style.display = "block";
      tip.style.left = Math.min(e.clientX - r.left + 12, Math.max(0, r.width - tip.offsetWidth - 4)) + "px";
      tip.style.top = (e.clientY - r.top + 14) + "px";
    });
    root.addEventListener("click", function (e) {
      var g = target(e); if (!g) return;
      var k = g._key;
      if (pinned[k]) delete pinned[k]; else pinned[k] = true;
      refresh();
    });
    root.addEventListener("dblclick", function (e) {
      if (!target(e)) { pinned = {}; refresh(); }
    });
    document.addEventListener("keydown", function (e) {
      if (e.key === "Escape" && document.body.contains(root)) { pinned = {}; match = {}; if (inp) inp.value = ""; if (count) count.textContent = ""; refresh(); }
    });
  }

  HTMLWidgets.widget({
    name: "mvremPlot",
    type: "output",
    factory: function (el) {
      return { renderValue: function (x) { render(el, x); },
               resize: function () { el.style.height = "auto"; } };
    }
  });
})();
