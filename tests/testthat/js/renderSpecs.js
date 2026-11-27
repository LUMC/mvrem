// Render widget specs with inst/htmlwidgets/mvremPlot.js on the DOM stub (domStub.js) and
// exercise the interactions; prints a JSON array with one result per spec.
//
//   node renderSpecs.js <mvremPlot.js> <spec.json> [<spec.json> ...]
//
// Each spec file holds the `x` of an mvremPlot widget (as htmlwidgets serialises it). Per
// spec the result is {ok: true, ...counts and states} or {ok: false, error, stack}. The exit
// code is 1 when any spec fails. Generic: a new widget kind is covered by adding a spec.
"use strict";
var fs = require("fs"), vm = require("vm"), path = require("path");

var args = process.argv.slice(2);
if (args.length < 2) {
  process.stderr.write("usage: node renderSpecs.js <mvremPlot.js> <spec.json>...\n");
  process.exit(2);
}
var stub = require(path.join(__dirname, "domStub.js"));
var jsPath = args[0];
vm.runInThisContext(fs.readFileSync(jsPath, "utf8"), { filename: jsPath });
if (stub.HTMLWidgets.definitions.length !== 1) throw new Error("expected one HTMLWidgets.widget() definition");
var def = stub.HTMLWidgets.definitions[0];

// an id to search for: the first subject / term id found in the spec
function firstId(spec) {
  var panels = spec.panels || [];
  for (var i = 0; i < panels.length; i++) {
    var p = panels[i];
    if (p.lines && p.lines.length) return String(p.lines[0].id);
    if (p.intervals && p.intervals.id !== undefined) return String([].concat(p.intervals.id)[0]);
    if (p.points && p.points.id !== undefined) return String([].concat(p.points.id)[0]);
    if (p.curves) for (var j = 0; j < p.curves.length; j++) if (p.curves[j].id !== undefined) return String(p.curves[j].id);
  }
  return null;
}

function run(file) {
  var spec = JSON.parse(fs.readFileSync(file, "utf8"));
  var el = stub.document.createElement("div");
  stub.document.body.appendChild(el);
  var w = def.factory(el, 600, 400);
  w.renderValue(spec);
  w.resize(600, 400);
  var root = el.querySelector("div.mvp"), grid = el.querySelector("div.grid"), tip = el.querySelector("div.tip");
  var res = {
    ok: true, name: def.name,
    panels: grid.children.length,                                  // one SVG per panel
    groups: el.querySelectorAll("g.el").length,                    // highlightable elements
    legend: el.querySelectorAll("span.lg").length,
    hasSearch: !!el.querySelector("input"),
    css: !!stub.document.getElementById("mvremPlotCss")
  };
  // every legend kind has the rule that hides its elements when the entry is switched off
  var kindCss = stub.document.getElementById("mvremPlotKindCss");
  res.kindRules = [].concat(spec.legend || []).every(function (e) {
    return !!kindCss && kindCss.textContent.indexOf(".mvp.off-" + e.kind + " .kind-" + e.kind + "{display:none}") >= 0;
  });

  // search box: the first id matches at least itself; a match dims the others; the dropdown
  // lists the matches, a click on an entry (or Enter) highlights exactly that entry
  var inp = el.querySelector("input"), count = el.querySelector("span.scount");
  var dd = el.querySelector("div.sdd"), id = firstId(spec);
  function shown(n) { return n && n.style.display !== "none"; }
  function nMatches() { return parseInt(count.textContent, 10); }
  if (inp && id !== null) {
    inp.value = id;
    inp.dispatch("input");
    res.searchMatches = nMatches();
    res.dimAfterSearch = root.classList.contains("dim");
    res.dropdownShown = shown(dd);
    res.dropdownItems = dd.children.length;
    res.dropdownFirst = dd.children.length ? dd.children[0].textContent : null;
    dd.children[0].dispatch("mousedown", {});
    res.pickedMatches = nMatches();
    res.pickedValue = inp.value;
    res.dropdownClosedAfterPick = !shown(dd);
    inp.dispatch("input");                         // keyboard: down + Enter picks the first entry
    inp.dispatch("keydown", { key: "ArrowDown" });
    res.activeAfterDown = dd.children[0].classList.contains("act");
    inp.dispatch("keydown", { key: "Enter" });
    res.enterMatches = nMatches();
    inp.dispatch("input");
    inp.dispatch("keydown", { key: "Escape" });
    res.dropdownClosedByEscape = !shown(dd);
    res.placeholder = inp.placeholder;
    // a value of a subject-level groupBy variable: one group entry per matching level
    var gb = spec.groupBy;
    if (gb && [].concat(gb.vars || []).length) {
      res.groupSearch = [].concat(gb.vars).map(function (v) {
        var lv = String([].concat(gb.levels[v])[0]);
        inp.value = lv; inp.dispatch("input");
        var grp = el.querySelectorAll("div.grp").map(function (g) { return g.textContent; });
        var out = { variable: v, level: lv, entries: grp, matches: nMatches() };
        var mine = el.querySelectorAll("div.grp").filter(function (g) { return g.textContent.indexOf(v + " = " + lv) === 0; })[0];
        if (mine) { mine.dispatch("click", {}); out.picked = nMatches(); }
        return out;
      });
    }
    inp.value = "";
    inp.dispatch("input");
    res.dimAfterClear = root.classList.contains("dim");
    res.dropdownHiddenWhenEmpty = !shown(dd);
  }

  // legend: a click switches the entry's kind off (class off-<kind> on the root), a second on
  var items = el.querySelectorAll("span.lg"), kinds = [].concat(spec.legend || []).map(function (e) { return e.kind; });
  res.legendOff = items.map(function (item, i) {
    item.dispatch("click");
    var off = item.classList.contains("off") && root.classList.contains("off-" + kinds[i]);
    item.dispatch("click");
    return off && !root.classList.contains("off-" + kinds[i]);
  }).every(Boolean);

  // colorBy / shapeBy (spec.groupBy): the first dropdown lists "gray" + the variables; choosing
  // the first colours points with the palette, a click on a level of its legend hides that subset
  // (some series get display none), a second click shows it, "gray" restores the default colours.
  // The second dropdown (shape) turns points into symbols (path.sym) and back.
  var sels = el.querySelectorAll("select"), sel = sels[0], gb = spec.groupBy;
  if (sel && gb) {
    var pal = [].concat(gb.palette), fills = function () {
      return el.querySelectorAll("circle").map(function (c) { return c.getAttribute("fill"); });
    };
    var hiddenSeries = function () {
      return el.querySelectorAll("g.el").filter(function (g) { return g.style.display === "none"; }).length;
    };
    var before = fills();
    res.colorOptions = sel.children.length;
    sel.value = [].concat(gb.vars)[0]; sel.dispatch("change");
    res.colored = fills().some(function (f) { return pal.indexOf(f) >= 0; });
    // opacity of the coloured lines and points per element kind (`fade`: observations behind predictions)
    res.colorOpacity = {};
    el.querySelectorAll("g.el").forEach(function (g) {
      var k = (g.getAttribute("class").match(/kind-(\S+)/) || [])[1] || "?";
      g.children.forEach(function (c) {
        var cl = c.getAttribute("class") || "";
        if (!/\bvis\b/.test(cl)) return;
        var o = c.getAttribute(c.tagName === "path" && !/sym/.test(cl) ? "stroke-opacity" : "fill-opacity");
        if (o === null || o === undefined) return;
        var a = res.colorOpacity[k] || (res.colorOpacity[k] = []);
        if (a.indexOf(String(o)) < 0) a.push(String(o));
      });
    });
    var levBox = el.querySelectorAll("span.clev")[0];
    res.colorLevels = levBox ? levBox.children.length : 0;
    if (levBox && levBox.children.length) {
      var lev = levBox.children[0];
      lev.dispatch("click");
      res.levelHides = hiddenSeries() > 0 && lev.classList.contains("off");
      lev.dispatch("click");
      res.levelShows = hiddenSeries() === 0 && !lev.classList.contains("off");
    }
    sel.value = ""; sel.dispatch("change");
    res.colorReset = JSON.stringify(fills()) === JSON.stringify(before);
    var shp = sels[1];
    res.shapeOptions = shp ? shp.children.length : 0;
    if (shp) {
      shp.value = [].concat(gb.shapeVars)[0]; shp.dispatch("change");
      var syms = function () {
        return el.querySelectorAll("path.sym").filter(function (s) { return s.style.display !== "none"; }).length;
      };
      res.shaped = syms() > 0;
      shp.value = ""; shp.dispatch("change");
      res.shapeReset = syms() === 0;
    }
  }

  // mouse: move over the first element -> tooltip + highlight; click pins it; Escape clears
  var g = el.querySelector("g.el");
  if (g) {
    var target = g.children[0] || g;
    root.dispatch("mousemove", { target: target, clientX: 100, clientY: 100 });
    res.tipShown = tip.style.display === "block";
    res.tipText = tip.textContent;
    res.highlighted = g.classList.contains("hl") && root.classList.contains("dim");
    root.dispatch("click", { target: target });
    root.dispatch("mouseleave", {});
    res.tipHiddenAfterLeave = tip.style.display === "none";
    res.pinnedAfterLeave = g.classList.contains("hl");
    stub.document.dispatch("keydown", { key: "Escape" });
    res.clearedByEscape = !g.classList.contains("hl") && !root.classList.contains("dim");
  }
  return res;
}

var results = args.slice(1).map(function (file) {
  try { return run(file); }
  catch (e) { return { ok: false, error: String(e && e.message || e), stack: String(e && e.stack || "") }; }
});
process.stdout.write(JSON.stringify(results) + "\n");
process.exit(results.every(function (r) { return r.ok; }) ? 0 : 1);
