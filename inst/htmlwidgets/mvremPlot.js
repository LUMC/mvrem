// mvrem interactive plots: panels of subject lines / intervals / points / curves drawn as
// SVG, with a tooltip and a highlight of the subject under the mouse that is linked across
// all panels. Self-contained: no external libraries.
(function () {
  var NS = "http://www.w3.org/2000/svg", W = 440, H = 300;
  // tooltip word per line kind ("S01 observed"); other kinds show the id only
  var KIND = { obs: "observed", fit: "fitted", res: "residual", pred: "predicted", predInt: "95% bound",
               fitMissing: "predicted (missing observation)", censored: "censored" };
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
    // search box: a dropdown of the matches (subjects and, by a subject-level variable, groups)
    ".mvp .sbox{position:relative;display:inline-block}",
    ".mvp .sdd{position:absolute;left:0;top:100%;z-index:6;min-width:100%;max-height:18em;overflow-y:auto;",
    "background:#fff;border:1px solid #bbb;border-radius:3px;box-shadow:0 2px 6px rgba(0,0,0,.15)}",
    ".mvp .sdd div{padding:1px 6px;cursor:pointer;white-space:nowrap}",
    ".mvp .sdd div.act{background:#e3eef8}",
    ".mvp .sdd div.grp{font-weight:bold}",
    ".mvp .sdd div.more{cursor:default;color:#777;font-style:italic}",
    ".mvp .sdd small{color:#777;margin-left:6px}",
    ".mvp select{font:inherit;padding:0 2px;border:1px solid #bbb;border-radius:3px}",
    ".mvp .cby{display:inline-flex;align-items:center;gap:4px;white-space:nowrap}",
    ".mvp .clev{display:inline-flex;flex-wrap:wrap;gap:2px 10px}",
    ".mvp .clev span{display:inline-flex;align-items:center;gap:3px}",
    ".mvp .clev i{display:inline-block;width:9px;height:9px;border-radius:50%}",
    ".mvp .clev span{cursor:pointer;user-select:none}",
    ".mvp .clev span.off{opacity:.35;text-decoration:line-through}",
    ".mvp .clev b{font-weight:normal;color:#555}",
    ".mvp .clev svg{width:11px;height:11px;flex:none}",
    ".mvp .tip{position:absolute;pointer-events:none;background:rgba(30,30,30,.92);color:#fff;",
    "padding:3px 6px;border-radius:3px;white-space:pre;display:none;font-size:11px;z-index:5}",
    ".mvp g.el{cursor:pointer}",
    ".mvp .hit{fill:none;stroke:transparent;stroke-width:9;pointer-events:stroke}",
    ".mvp circle.hit{fill:transparent;stroke:none;pointer-events:all}",
    ".mvp.dim g.el:not(.hl){opacity:.12}",
    // highlighted subject: a softer red for the lines, the full red for its observed points
    ".mvp g.el.hl path.vis:not(.sym),.mvp g.el.hl line.vis{stroke:#efb8b0 !important;stroke-width:2 !important}",
    ".mvp g.el.hl circle.vis,.mvp g.el.hl path.vis.sym{fill:#c0392b !important}",
    // (faded series, `fade`, are fully opaque when highlighted)
    ".mvp g.el.hl .vis{stroke-opacity:1 !important;fill-opacity:1 !important}",
    // (shape symbols are filled paths: no line highlight on them)
    ".mvp g.el.hl path.vis.sym{stroke:none !important}",
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

  // Point symbols of shapeBy (the names of plotShapes in R/plotUtils.R): SVG path data of a
  // filled symbol of "radius" r at (x, y); "circle" is drawn as a <circle>.
  function symbolPath(shape, x, y, r) {
    function poly(pts, ang) {
      var c = Math.cos(ang), s = Math.sin(ang);
      return pts.map(function (p, i) {
        return (i ? "L" : "M") + (x + p[0] * c - p[1] * s).toFixed(1) + "," + (y + p[0] * s + p[1] * c).toFixed(1);
      }).join("") + "Z";
    }
    var a = 1.3 * r, w = 0.45 * r, q = 0.9 * r;
    var plus = [[-w, -a], [w, -a], [w, -w], [a, -w], [a, w], [w, w], [w, a], [-w, a], [-w, w], [-a, w], [-a, -w], [-w, -w]];
    switch (shape) {
      case "square": return poly([[-q, -q], [q, -q], [q, q], [-q, q]], 0);
      case "triangle": return poly([[0, -1.25 * r], [1.1 * r, 0.75 * r], [-1.1 * r, 0.75 * r]], 0);
      case "diamond": return poly([[0, -a], [a, 0], [0, a], [-a, 0]], 0);
      case "plus": return poly(plus, 0);
      case "cross": return poly(plus, Math.PI / 4);
      default: return null;
    }
  }
  // A point of a series under colorBy / shapeBy: a circle (kept for the default shape) and,
  // made when first needed, a symbol path; setPoint() shows one of them.
  function setPoint(pt, shape, fill) {
    var d = shape && shape !== "circle" ? symbolPath(shape, pt.x, pt.y, pt.r * 1.25) : null;
    if (d) {
      if (!pt.s) { pt.s = S("path", { "class": "vis sym" }); pt.c.parentNode.insertBefore(pt.s, pt.c); }
      pt.s.setAttribute("d", d); pt.s.setAttribute("fill", fill);
      pt.s.style.display = ""; pt.c.style.display = "none";
    } else {
      if (pt.s) pt.s.style.display = "none";
      pt.c.style.display = "";
      if (!pt.hollow) pt.c.setAttribute("fill", fill);
    }
  }
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

    // frame, grid, ticks, labels: a grid line + label per tick inside the range of the axis
    function axisTicks(r, n, horizontal) {
      ticks(r[0], r[1], n).forEach(function (v) {
        if (v < r[0] || v > r[1]) return;
        if (horizontal) {
          S("line", { x1: sx(v), x2: sx(v), y1: T, y2: B, stroke: "#ececec" }, svg);
          S("text", { x: sx(v), y: B + 13, "text-anchor": "middle", "font-size": 10 }, svg).textContent = fmt(v);
        } else {
          S("line", { x1: L, x2: R, y1: sy(v), y2: sy(v), stroke: "#ececec" }, svg);
          S("text", { x: L - 5, y: sy(v) + 3, "text-anchor": "end", "font-size": 10 }, svg).textContent = fmt(v);
        }
      });
    }
    axisTicks(xr, 6, true);
    if (ylabels) {
      ylabels.pos.forEach(function (v, i) {
        S("text", { x: L - 5, y: sy(v) + 3, "text-anchor": "end", "font-size": 10 }, svg).textContent = ylabels.text[i];
      });
    } else {
      axisTicks(yr, 5, false);
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
    // an element group: `_id` / `_key` (the subject or curve id, linked across panels when
    // ctx.linked), `_info` (tooltip data), `_svg` / `_ix` / `_panel` (for the hover position)
    function group(id, kind, extra) {
      var g = S("g", { "class": "el kind-" + kind }, plot);
      g._id = id; g._info = extra; g._svg = svg; g._ix = ix; g._panel = p;
      g._key = ctx.linked ? String(id) : pi + "|" + id;
      ctx.add(g._key, g);
      return g;
    }
    // colorBy / shapeBy state of a group: `_cv` (values of the grouping variables per point),
    // `_col` (own colour), `_pts` (its points: {c, x, y, r, hollow, s}), and the options
    // `_fade` (opacity of a coloured line / points), `_tint` (hollow circles coloured too),
    // `_line` (its path, recoloured when all points share one level)
    function styled(g, cv, col, opts) {
      g._cv = cv; g._col = col; g._pts = [];
      opts = opts || {};
      if (opts.fade !== undefined) g._fade = opts.fade;
      if (opts.tint) g._tint = true;
      if (opts.line) g._line = opts.line;
    }
    // a point of a group: the circle (hollow: open), registered for colorBy / shapeBy when the
    // group has `_cv`
    function addPoint(g, cx, cy, r, attrs, hollow) {
      attrs.cx = cx; attrs.cy = cy; attrs.r = r;
      var c = S("circle", attrs, g);
      if (g._cv) g._pts.push(hollow ? { c: c, x: cx, y: cy, r: r, hollow: true } : { c: c, x: cx, y: cy, r: r });
      return c;
    }
    // invisible hit areas: a circle around a point, a path along a line
    function hitCircle(g, cx, cy, r) { S("circle", { cx: cx, cy: cy, r: r || 6, "class": "hit" }, g); }
    function hitPath(g, d) { S("path", { d: d, "class": "hit" }, g); }

    // subject lines (observed / residual: solid with points; fitted: dotted). Only "fit" has
    // the fitted defaults; any other kind is drawn like observed values unless the line brings
    // its own color / dash / width / points (as the predictions do).
    lines.forEach(function (l) {
      var xv = arr(l.x), yv = arr(l.y), obs = l.kind !== "fit";
      if (l.hollow) {                    // open circles, never joined (predicted missing observations)
        var gh = group(l.id, l.kind, { x: xv, y: yv, label: String(l.id) + " " + (KIND[l.kind] || ""), p: p });
        // hidden by level; `tint`: also coloured by the level (censored times on survival curves)
        if (l.cv) styled(gh, l.cv, l.color || COL.accent, { tint: l.tint });
        xv.forEach(function (v, i) {
          addPoint(gh, sx(v), sy(yv[i]), 3, { fill: "none", stroke: l.color || COL.accent, "stroke-width": 1.2,
                                              "class": "vis hollow" }, true);
          hitCircle(gh, sx(v), sy(yv[i]));
        });
        return;
      }
      var col = l.color || (obs ? COL.subject : COL.accentLight), dash = l.dash !== undefined ? l.dash : (obs ? "" : "2 3");
      var g = group(l.id, l.kind || "obs", { x: xv, y: yv, label: String(l.id) + " " + (KIND[l.kind] || ""), p: p });
      // colorBy: values per point (`fade`: opacity of the coloured line and points, e.g. the
      // observations behind predictions)
      if (l.cv) styled(g, l.cv, col, { fade: l.fade });
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
          addPoint(g, sx(v), sy(yv[i]), 1.8 * (isFinite(f) ? f : 1), { fill: col, "class": l.size ? "vis wt" : "vis" });
        });
      if (xv.length > 1 && p.connect !== false) hitPath(g, path(xv, yv));
      else xv.forEach(function (v, i) { hitCircle(g, sx(v), sy(yv[i])); });
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
        if (pts.cv) {                    // colorBy / shapeBy: the values of this point
          var cv1 = {}; for (var vn in pts.cv) cv1[vn] = [arr(pts.cv[vn])[i]];
          styled(g, cv1, COL.accent);
        }
        addPoint(g, sx(px[i]), sy(py[i]), 2.6, { fill: COL.accent, "fill-opacity": 0.75, "class": "vis" });
        hitCircle(g, sx(px[i]), sy(py[i]), 7);
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
        attrs["class"] = "vis"; var cp = S("path", attrs, g);
        hitPath(g, attrs.d);
        // a curve of one level of a grouping variable (survival curves by level): coloured and
        // hidden with that level by the colour dropdown, like the subject series
        if (c.cv) styled(g, c.cv, attrs.stroke, { line: cp, fade: c.fade });
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
    var bar = H_("div", "bar", root), count = null, inp = null, closeSearch = null, searchHint = null;
    if (spec.search) {
      // Typing highlights the subjects whose id contains the text, or whose value of a
      // subject-level variable of colorBy / shapeBy (e.g. clinic) does; a dropdown lists them,
      // the groups ("clinic = C01") first. Choosing an entry highlights exactly that subject or
      // group (arrow keys + Enter, or a click); Esc closes the list.
      var sbox = H_("span", "sbox", bar);
      inp = H_("input", "", sbox); inp.type = "search"; inp.placeholder = "highlight " + spec.search;
      var dd = H_("div", "sdd", sbox); dd.style.display = "none";
      count = H_("span", "scount", bar);
      var entries = [], active = -1, attrs = null, MAXITEMS = 30;
      var labelOf = function (k) { return k.indexOf("|") >= 0 && !ctx.linked ? k.split("|").slice(1).join("|") : k; };
      // values of the groupBy variables that are subject-level: constant within every subject
      // (over all its series); a variable varying within some subject is not searched
      var subjectAttrs = function () {
        if (attrs) return attrs;
        attrs = {};
        var vars = spec.groupBy ? arr(spec.groupBy.vars) : [], varies = {};
        ids.forEach(function (k) {
          var a = {};
          vars.forEach(function (v) {
            var seen = {}, n = 0, last = null;
            byId[k].forEach(function (g) {
              if (!g._cv || g._cv[v] === undefined) return;
              arr(g._cv[v]).forEach(function (x) {
                if (x === null || x === undefined) return;
                if (!seen[String(x)]) { seen[String(x)] = true; n++; last = String(x); }
              });
            });
            if (n > 1 || last === labelOf(k)) varies[v] = true;   // (or the subject id itself)
            else if (n === 1) a[v] = last;
          });
          attrs[k] = a;
        });
        ids.forEach(function (k) { for (var v in varies) delete attrs[k][v]; });
        return attrs;
      };
      var setMatch = function (keys) {
        match = {}; keys.forEach(function (k) { match[k] = true; });
        count.textContent = inp.value.trim() ? keys.length + " match" + (keys.length === 1 ? "" : "es") : "";
        refresh();
      };
      var mark = function (i) {
        active = i;
        [].forEach.call(dd.children, function (c, j) { c.classList.toggle("act", j === i); });
      };
      var choose = function (en) {
        inp.value = en.text; dd.style.display = "none"; active = -1;
        setMatch(en.keys);
      };
      var item = function (en, i) {
        var it = H_("div", en.grp ? "grp" : "", dd, en.label);
        if (en.note) H_("small", "", it, en.note);
        var pick = function (e) { if (e && e.preventDefault) e.preventDefault(); choose(en); };
        it.addEventListener("mousedown", pick);
        it.addEventListener("click", pick);
        it.addEventListener("mousemove", function () { mark(i); });
      };
      var update = function () {
        var q = inp.value.trim().toLowerCase(), A = subjectAttrs(), groups = {}, subj = [];
        entries = [];
        if (q) {
          ids.forEach(function (k) {
            var hit = labelOf(k).toLowerCase().indexOf(q) >= 0;
            for (var v in A[k]) {
              if (String(A[k][v]).toLowerCase().indexOf(q) < 0) continue;
              hit = true;
              var gk = v + "\u0000" + A[k][v];
              (groups[gk] = groups[gk] || { v: v, lv: A[k][v], keys: [] }).keys.push(k);
            }
            if (hit) subj.push(k);
          });
          for (var gk in groups) {
            var gr = groups[gk];
            entries.push({ grp: true, label: gr.v + " = " + gr.lv, text: gr.lv, keys: gr.keys,
                           note: gr.keys.length + " subject" + (gr.keys.length === 1 ? "" : "s") });
          }
          subj.forEach(function (k) {
            var note = [];
            for (var v in A[k]) note.push(v + " " + A[k][v]);
            entries.push({ label: labelOf(k), text: labelOf(k), keys: [k], note: note.join(", ") });
          });
        }
        setMatch(subj);
        dd.innerHTML = ""; active = -1;
        var shown = entries.slice(0, MAXITEMS);
        shown.forEach(item);
        if (entries.length > shown.length) H_("div", "more", dd, "... " + (entries.length - shown.length) + " more");
        entries = shown;
        dd.style.display = entries.length ? "" : "none";
      };
      // the placeholder names the subject-level variables that can be searched too
      searchHint = function () {
        var A = subjectAttrs(), vs = [];
        ids.forEach(function (k) { for (var v in A[k]) if (vs.indexOf(v) < 0) vs.push(v); });
        if (vs.length) inp.placeholder = "highlight " + spec.search + " or " + vs.join(", ");
      };
      closeSearch = function () { dd.style.display = "none"; active = -1; };
      inp.addEventListener("input", update);
      inp.addEventListener("focus", function () { if (inp.value.trim()) update(); });
      inp.addEventListener("blur", closeSearch);
      inp.addEventListener("keydown", function (e) {
        if (dd.style.display === "none" || !entries.length) return;
        if (e.key === "ArrowDown" || e.key === "ArrowUp") {
          if (e.preventDefault) e.preventDefault();
          var n = entries.length;
          mark(e.key === "ArrowDown" ? (active + 1) % n : (active - 1 + n) % n);
        } else if (e.key === "Enter") {
          choose(entries[active >= 0 ? active : 0]);
        } else if (e.key === "Escape") {
          closeSearch();
          if (e.stopPropagation) e.stopPropagation();
        }
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

    // colorBy / shapeBy (groupSpec() in R): two dropdowns, "color: gray / <variable>" and
    // "shape: circle / <variable>", each followed by the legend of the chosen variable's levels.
    // A click on a level hides that subset: the points of the level and, when all points of a
    // series (a subject's line, prediction, ...) are hidden, the whole series; a second click
    // shows it again. Choosing another variable shows everything again.
    var gb = spec.groupBy, styleTo = null;
    if (gb && arr(gb.vars).length) {
      var pal = arr(gb.palette), shapes = arr(gb.shapes);
      var st = { color: "", shape: "", offC: {}, offS: {} };
      var levIx = function (v, x) {
        return x === null || x === undefined ? -1 : arr(gb.levels[v]).map(String).indexOf(String(x));
      };
      var colorFor = function (v, x) { var k = levIx(v, x); return k < 0 ? null : pal[k % pal.length]; };
      var shapeFor = function (v, x) { var k = levIx(v, x); return k < 0 ? "circle" : shapes[k % shapes.length]; };
      // the node of a point that is currently shown: its symbol (shapeBy) if displayed, else its circle
      var visibleNode = function (pt) { return pt.s && pt.s.style.display !== "none" ? pt.s : pt.c; };
      var apply = function () {
        for (var k in byId) byId[k].forEach(function (g) {
          if (!g._cv) return;
          var cv = st.color && g._cv[st.color] !== undefined ? arr(g._cv[st.color]) : null;
          var sv = st.shape && g._cv[st.shape] !== undefined ? arr(g._cv[st.shape]) : null;
          var n = (cv || sv || []).length, hidden = 0;
          var off = function (i) {
            return (cv && st.offC[String(cv[i])]) || (sv && st.offS[String(sv[i])]) ? true : false;
          };
          var cols = cv ? cv.map(function (x) { return colorFor(st.color, x); }) : [];
          g._pts.forEach(function (pt, i) {
            if (pt.hollow) {                                   // stays an open circle
              pt.c.style.display = off(i) ? "none" : "";
              if (g._tint) pt.c.setAttribute("stroke", cols[i] || g._col);
              return;
            }
            setPoint(pt, sv ? shapeFor(st.shape, sv[i]) : "circle", cols[i] || g._col);
            var o = off(i);
            if (g._fade !== undefined) visibleNode(pt).setAttribute("fill-opacity", cv ? g._fade : 1);
            visibleNode(pt).style.display = o ? "none" : "";
            if (o && pt.s) pt.c.style.display = "none";
          });
          for (var i = 0; i < n; i++) if (off(i)) hidden++;
          g.style.display = n && hidden === n ? "none" : "";     // the whole series when all is hidden
          if (g._line) {
            var one = cv && cols.length && cols.every(function (c) { return c && c === cols[0]; });
            g._line.setAttribute("stroke", one ? cols[0] : g._col);
            g._line.setAttribute("stroke-opacity", one ? (g._fade !== undefined ? g._fade : 0.7) : 1);
          }
        });
      };
      var swatch = function (box, which, v, lv) {
        if (which === "color") { var dot = H_("i", "", box); dot.style.background = colorFor(v, lv); return; }
        var sw = S("svg", { viewBox: "0 0 12 12" }, box), shp = shapeFor(v, lv), d = symbolPath(shp, 6, 6, 3.6);
        if (d) S("path", { d: d, fill: "#555" }, sw); else S("circle", { cx: 6, cy: 6, r: 4, fill: "#555" }, sw);
      };
      var control = function (which, vars, none) {
        var box = H_("span", "cby", bar);
        H_("span", "", box, which + ":");
        var sel = H_("select", "", box), levBox = H_("span", "clev", bar);
        [""].concat(vars).forEach(function (v) { var o = H_("option", "", sel, v === "" ? none : v); o.value = v; });
        var choose = function (v) {
          sel.value = v || ""; st[which] = v || "";
          var offs = which === "color" ? (st.offC = {}) : (st.offS = {});
          levBox.innerHTML = "";
          if (v) arr(gb.levels[v]).forEach(function (lv) {
            var it = H_("span", "", levBox); swatch(it, which, v, lv);
            it.appendChild(document.createTextNode(String(lv)));
            it.title = "click: hide / show " + v + " = " + lv;
            it.addEventListener("click", function () {
              var o = !offs[String(lv)];
              if (o) offs[String(lv)] = true; else delete offs[String(lv)];
              it.classList.toggle("off", o); apply();
            });
          });
          apply();
        };
        sel.addEventListener("change", function () { choose(sel.value); });
        return choose;
      };
      var chooseColor = control("color", arr(gb.vars), "gray (default)");
      var chooseShape = arr(gb.shapeVars).length ? control("shape", arr(gb.shapeVars), "circle (default)") : null;
      styleTo = function () {
        if (typeof gb.color === "string" && gb.color) chooseColor(gb.color);
        if (chooseShape && typeof gb.shape === "string" && gb.shape) chooseShape(gb.shape);
      };
    }

    var grid = H_("div", "grid", root);
    grid.style.gridTemplateColumns = "repeat(" + (spec.ncol || 1) + ",minmax(0,1fr))";
    arr(spec.panels).forEach(function (p, i) { grid.appendChild(drawPanel(p, i, ctx)); });
    if (styleTo) styleTo();
    if (searchHint) searchHint();

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
      if (e.key === "Escape" && document.body.contains(root)) {
        pinned = {}; match = {}; if (inp) inp.value = ""; if (count) count.textContent = "";
        if (closeSearch) closeSearch();
        refresh();
      }
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
