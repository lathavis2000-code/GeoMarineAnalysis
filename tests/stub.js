// Minimal Earth Engine / UI stubs so the whole tool file can be LOADED and its
// pure-JS functions exercised in Node, with no Earth Engine account and no
// network. Nothing here simulates Earth Engine semantics - every ee.* call
// returns a chainable no-op. That is deliberate: these tests cover the
// statistics, the guards and the wiring, which are plain JavaScript, and they
// make no claim about anything that depends on a real reduceRegions result.
'use strict';

function chain() {
  var f = function () { return chain(); };
  return new Proxy(f, {
    get: function (t, k) {
      if (k === 'then' || k === 'toJSON' || k === Symbol.iterator) return undefined;
      // A chainable must be coercible to a primitive. Returning undefined for
      // Symbol.toPrimitive sends JS to valueOf/toString, which this proxy also
      // answers with a proxy, so ANY string concatenation of an ee value threw
      // "Cannot convert object to primitive value" - and a script that builds a
      // description or an assetId by concatenation could not be loaded at all.
      if (k === Symbol.toPrimitive) {
        return function (hint) { return hint === 'number' ? 0 : '[ee]'; };
      }
      if (k === 'evaluate') return function (cb) { if (cb) cb(null, null); };
      if (k === 'getInfo') return function () { return null; };
      return chain();
    },
    apply: function () { return chain(); }
  });
}

function makeWidget(kind, init) {
  var o = {
    _kind: kind,
    _value: (init && init.value !== undefined) ? init.value : (typeof init === 'string' ? init : ''),
    _style: {},
    setValue: function (v) { this._value = v; return this; },
    getValue: function () { return this._value; },
    style: function () {
      var s = this._style;
      return {
        set: function (a, b) { if (typeof a === 'object') { for (var k in a) s[k] = a[k]; } else s[a] = b; },
        get: function (k) { return s[k]; }
      };
    },
    add: function () { return this; },
    insert: function () { return this; },
    clear: function () { return this; },
    widgets: function () {
      return { get: function () { return makeWidget('x'); }, length: 0, reset: function () {}, add: function () {} };
    },
    setDisabled: function () { return this; },
    onClick: function () { return this; },
    onChange: function () { return this; },
    setMap: function () { return this; },
    setOptions: function () { return this; },
    setChartType: function () { return this; },
    setSeriesNames: function () { return this; },
    setLayout: function () { return this; }
  };
  return o;
}

// Installs the stubs as globals and returns a handle for reading them back.
function install(g) {
  g = g || global;
  var prints = [];
  g.ee = chain();
  g.ui = {
    Label: function (t) { return makeWidget('Label', t); },
    Button: function (o) { var w = makeWidget('Button'); w._cb = o && o.onClick; return w; },
    Textbox: function (o) { return makeWidget('Textbox', (o && o.value) || ''); },
    Select: function (o) { return makeWidget('Select', (o && o.value) || ''); },
    Checkbox: function () { return makeWidget('Checkbox'); },
    Slider: function () { return makeWidget('Slider'); },
    Chart: function () { return makeWidget('Chart'); },
    Map: { Layer: function () { return makeWidget('Layer'); } },
    Panel: function () { return makeWidget('Panel'); },
    root: { insert: function () {}, add: function () {}, clear: function () {} },
    util: { clear: function () {}, setTimeout: function () { return 0; }, clearTimeout: function () {} }
  };
  g.ui.Panel.Layout = { Flow: function () { return {}; }, absolute: function () { return {}; } };
  g.ui.Chart.array = { values: function () { return makeWidget('Chart'); } };
  g.ui.Chart.feature = { byFeature: function () { return makeWidget('Chart'); } };
  g.Map = {
    onClick: function () {}, add: function () {}, remove: function () {},
    setCenter: function () {}, addLayer: function () {}, clear: function () {},
    layers: function () {
      return { reset: function () {}, add: function () {}, length: 0, get: function () { return makeWidget('L'); } };
    },
    setOptions: function () {}, style: function () { return { set: function () {} }; },
    drawingTools: function () { return { setShown: function () {} }; },
    setControlVisibility: function () {}, centerObject: function () {}, unlisten: function () {}
  };
  g.print = function () { prints.push(Array.prototype.slice.call(arguments).join(' ')); };
  // Every Export target the Code Editor offers, recorded rather than ignored:
  // a script that only calls toAsset (the SAR study's STAGE A) could not be
  // loaded at all when toDrive was the only entry.
  g.Export = { exports: [] };
  ['table', 'image', 'video', 'map'].forEach(function (kind) {
    g.Export[kind] = {};
    ['toDrive', 'toAsset', 'toCloudStorage', 'toBigQuery', 'toDriveVideo'].forEach(function (dest) {
      g.Export[kind][dest] = function (opts) {
        g.Export.exports.push({ kind: kind, dest: dest, opts: opts || {} });
      };
    });
  });
  return { prints: prints };
}

module.exports = { install: install };
