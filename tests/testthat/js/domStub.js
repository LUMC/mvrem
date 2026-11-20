// A minimal DOM stub for running inst/htmlwidgets/mvremPlot.js in node, without a browser:
// only what the widget's JavaScript calls (no layout, no rendering, no CSS). Loading this file
// defines the globals `document`, `window` and `HTMLWidgets` (which records the widget
// definitions passed to HTMLWidgets.widget()). Elements keep their tag, attributes, class list,
// style, text, children and event listeners; `dispatch(type, event)` calls the listeners of an
// element (not a DOM method: the driver uses it to simulate events).
"use strict";

function ClassList(el) { this.el = el; }
ClassList.prototype.list = function () {
  return this.el.className ? this.el.className.split(/\s+/).filter(Boolean) : [];
};
ClassList.prototype.write = function (a) { this.el.className = a.join(" "); };
ClassList.prototype.contains = function (c) { return this.list().indexOf(c) >= 0; };
ClassList.prototype.add = function (c) { var a = this.list(); if (a.indexOf(c) < 0) a.push(c); this.write(a); };
ClassList.prototype.remove = function (c) { this.write(this.list().filter(function (x) { return x !== c; })); };
ClassList.prototype.toggle = function (c, force) {
  var on = force === undefined ? !this.contains(c) : !!force;
  if (on) this.add(c); else this.remove(c);
  return on;
};

function Node(tag, ns) {
  this.tagName = tag;
  this.namespaceURI = ns || null;
  this.attributes = {};
  this.children = [];
  this.parentNode = null;
  this.style = {};
  this.listeners = {};
  this.ownText = "";
  this.className = "";
  this.id = "";
  this.value = "";
  this.offsetWidth = 0;
  this.classList = new ClassList(this);
}
Node.prototype.setAttribute = function (k, v) {
  this.attributes[k] = String(v);
  if (k === "class") this.className = String(v);
  if (k === "id") this.id = String(v);
};
Node.prototype.getAttribute = function (k) {
  return Object.prototype.hasOwnProperty.call(this.attributes, k) ? this.attributes[k] : null;
};
Node.prototype.appendChild = function (c) {
  if (c.parentNode) c.parentNode.removeChild(c);
  c.parentNode = this;
  this.children.push(c);
  return c;
};
Node.prototype.removeChild = function (c) {
  var i = this.children.indexOf(c);
  if (i >= 0) this.children.splice(i, 1);
  c.parentNode = null;
  return c;
};
Object.defineProperty(Node.prototype, "textContent", {
  get: function () {
    return this.ownText + this.children.map(function (c) { return c.textContent; }).join("");
  },
  set: function (v) { this.children = []; this.ownText = String(v); }
});
Object.defineProperty(Node.prototype, "innerHTML", {
  get: function () { throw new Error("innerHTML is not serialised by the stub"); },
  set: function (v) {
    if (v !== "") throw new Error("the stub only supports clearing innerHTML");
    this.children = []; this.ownText = "";
  }
});
// simple selectors only: tag, .class, #id and tag.class
Node.prototype.matches = function (sel) {
  var m = /^([a-zA-Z0-9]*)(?:#([\w-]+))?(?:\.([\w-]+))?$/.exec(sel);
  if (!m) throw new Error("unsupported selector: " + sel);
  if (m[1] && m[1].toLowerCase() !== String(this.tagName).toLowerCase()) return false;
  if (m[2] && this.id !== m[2]) return false;
  if (m[3] && !this.classList.contains(m[3])) return false;
  return true;
};
Node.prototype.closest = function (sel) {
  for (var n = this; n; n = n.parentNode) if (n.matches && n.matches(sel)) return n;
  return null;
};
Node.prototype.contains = function (node) {
  for (var n = node; n; n = n.parentNode) if (n === this) return true;
  return false;
};
Node.prototype.querySelectorAll = function (sel) {
  var out = [];
  (function walk(n) {
    n.children.forEach(function (c) { if (c.matches && c.matches(sel)) out.push(c); walk(c); });
  })(this);
  return out;
};
Node.prototype.querySelector = function (sel) { return this.querySelectorAll(sel)[0] || null; };
Node.prototype.getBoundingClientRect = function () { return { left: 0, top: 0, width: 440, height: 300 }; };
Node.prototype.addEventListener = function (type, fn) { (this.listeners[type] = this.listeners[type] || []).push(fn); };
Node.prototype.dispatch = function (type, ev) {
  ev = ev || {};
  if (!ev.target) ev.target = this;
  (this.listeners[type] || []).forEach(function (fn) { fn(ev); });
};

function TextNode(text) { Node.call(this, "#text"); this.ownText = String(text); }
TextNode.prototype = Object.create(Node.prototype);
TextNode.prototype.matches = function () { return false; };

var documentElement = new Node("html");
var head = documentElement.appendChild(new Node("head"));
var body = documentElement.appendChild(new Node("body"));

var document = {
  documentElement: documentElement,
  head: head,
  body: body,
  listeners: {},
  createElement: function (tag) { return new Node(String(tag).toLowerCase()); },
  createElementNS: function (ns, tag) { return new Node(tag, ns); },
  createTextNode: function (text) { return new TextNode(text); },
  getElementById: function (id) { return documentElement.querySelector("#" + id); },
  querySelector: function (sel) { return documentElement.querySelector(sel); },
  querySelectorAll: function (sel) { return documentElement.querySelectorAll(sel); },
  addEventListener: Node.prototype.addEventListener,
  dispatch: Node.prototype.dispatch
};

var HTMLWidgets = {
  definitions: [],
  widget: function (def) { this.definitions.push(def); }
};

global.document = document;
global.window = { document: document };
global.HTMLWidgets = HTMLWidgets;

module.exports = { Node: Node, document: document, HTMLWidgets: HTMLWidgets };
