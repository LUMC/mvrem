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

  // search box: the first id matches at least itself; a match dims the others
  var inp = el.querySelector("input"), count = inp ? inp.parentNode.children[1] : null;
  var id = firstId(spec);
  if (inp && id !== null) {
    inp.value = id;
    inp.dispatch("input");
    res.searchMatches = parseInt(count.textContent, 10);
    res.dimAfterSearch = root.classList.contains("dim");
    inp.value = "";
    inp.dispatch("input");
    res.dimAfterClear = root.classList.contains("dim");
  }

  // legend: a click switches the entry's kind off (class off-<kind> on the root), a second on
  var items = el.querySelectorAll("span.lg"), kinds = [].concat(spec.legend || []).map(function (e) { return e.kind; });
  res.legendOff = items.map(function (item, i) {
    item.dispatch("click");
    var off = item.classList.contains("off") && root.classList.contains("off-" + kinds[i]);
    item.dispatch("click");
    return off && !root.classList.contains("off-" + kinds[i]);
  }).every(Boolean);

  // colorBy: the dropdown lists "gray" + the variables; choosing the first colours points with
  // the palette, choosing "gray" restores the default colours
  var sel = el.querySelector("select");
  if (sel && spec.colorBy) {
    var pal = [].concat(spec.colorBy.palette), fills = function () {
      return el.querySelectorAll("circle").map(function (c) { return c.getAttribute("fill"); });
    };
    var before = fills();
    res.colorOptions = sel.children.length;
    sel.value = [].concat(spec.colorBy.vars)[0]; sel.dispatch("change");
    res.colored = fills().some(function (f) { return pal.indexOf(f) >= 0; });
    var levBox = el.querySelector("span.clev");
    res.colorLevels = levBox ? levBox.children.length : 0;
    sel.value = ""; sel.dispatch("change");
    res.colorReset = JSON.stringify(fills()) === JSON.stringify(before);
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
