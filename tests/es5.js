// ES5 gate for the Earth Engine Code Editor sandbox.
//
// WHY THIS EXISTS: the GEE Code Editor does not support ES6. This file has been
// bitten twice - String.prototype.repeat() (v10.149) and Object.assign()
// (v10.158) - and both times the failure was at RUNTIME IN THE BROWSER, not at
// parse time, because ES6 METHODS are syntactically valid ES5. `node --check`
// cannot catch them. This scanner can.
//
// IT MUST STRIP STRINGS AS WELL AS COMMENTS. The tool file's own startup block
// prints changelog prose that discusses ES6 by name - "Object.assign (ES6)
// replaced with an explicit copy loop", "String.prototype.repeat() is ES6" -
// inside ordinary print() string literals. A naive grep flags all of those and
// fails on a file that is perfectly clean, which is worse than no gate at all.
'use strict';

// Walks the source once, tracking line/block comments and ' " ` strings, and
// returns the code with comments and string BODIES blanked (delimiters kept, so
// offsets and line numbers stay meaningful). Backtick openings are recorded
// separately, because a template literal is itself an ES6 construct and would
// otherwise be erased by the very step that strips it.
function stripCommentsAndStrings(src) {
  var out = new Array(src.length);
  var backticks = [];
  var i = 0, n = src.length, line = 1;
  var prevSignificant = '';
  while (i < n) {
    var c = src[i], d = src[i + 1];
    if (c === '\n') line++;
    if (c === '/' && d === '/') {
      while (i < n && src[i] !== '\n') { out[i] = ' '; i++; }
      continue;
    }
    if (c === '/' && d === '*') {
      out[i] = ' '; out[i + 1] = ' '; i += 2;
      while (i < n && !(src[i] === '*' && src[i + 1] === '/')) {
        // Keep the newline. scan() derives its line numbers by counting \n in the
        // STRIPPED source, so blanking them made every finding after the first
        // block comment report a wrong line - 236 lines off in the SAR script.
        if (src[i] === '\n') { line++; out[i] = '\n'; } else { out[i] = ' '; }
        i++;
      }
      if (i < n) { out[i] = ' '; out[i + 1] = ' '; i += 2; }
      continue;
    }
    if (c === '"' || c === "'" || c === '`') {
      if (c === '`') backticks.push(line);
      var q = c;
      out[i] = c; i++;
      while (i < n) {
        if (src[i] === '\\') {
          // A backslash-newline continuation is still a newline. Blanking it was
          // the last place the stripper lost one, and it cost a line on
          // everything after it. Both the string and regex branches do this.
          out[i] = ' ';
          if (i + 1 < n) {
            if (src[i + 1] === '\n') { line++; out[i + 1] = '\n'; } else { out[i + 1] = ' '; }
          }
          i += 2; continue;
        }
        if (src[i] === q) { out[i] = q; i++; break; }
        if (src[i] === '\n') { line++; out[i] = '\n'; i++; continue; }
        out[i] = ' '; i++;
      }
      prevSignificant = q;
      continue;
    }
    // Regex literal: a '/' in a position where a value may begin. Treated as a
    // literal so a quote INSIDE a character class cannot desynchronise the
    // string scanner above.
    if (c === '/' && '=(,:[!&|?{};\n'.indexOf(prevSignificant) !== -1) {
      out[i] = ' '; i++;
      var inClass = false;
      while (i < n) {
        if (src[i] === '\\') {
          // A backslash-newline continuation is still a newline. Blanking it was
          // the last place the stripper lost one, and it cost a line on
          // everything after it. Both the string and regex branches do this.
          out[i] = ' ';
          if (i + 1 < n) {
            if (src[i + 1] === '\n') { line++; out[i + 1] = '\n'; } else { out[i + 1] = ' '; }
          }
          i += 2; continue;
        }
        if (src[i] === '[') inClass = true;
        else if (src[i] === ']') inClass = false;
        else if (src[i] === '/' && !inClass) { out[i] = ' '; i++; break; }
        else if (src[i] === '\n') break;
        out[i] = ' '; i++;
      }
      continue;
    }
    out[i] = c;
    if (!/\s/.test(c)) prevSignificant = c;
    i++;
  }
  return { code: out.join(''), backticks: backticks };
}

var RULES = [
  ['const declaration', /(^|[^\w.$])const\s+[A-Za-z_$]/g],
  ['let declaration', /(^|[^\w.$])let\s+[A-Za-z_$]/g],
  ['arrow function', /=>/g],
  ['class declaration', /(^|[^\w.$])class\s+[A-Za-z_$][\w$]*\s*[{e]/g],
  ['for...of', /\bfor\s*\(\s*(?:var|let|const)?\s*[\w$]+\s+of\s/g],
  ['Object.assign', /\bObject\.assign\s*\(/g],
  ['Object.values/entries', /\bObject\.(?:values|entries)\s*\(/g],
  ['Array.from', /\bArray\.from\s*\(/g],
  ['Array.prototype.includes', /\.includes\s*\(/g],
  ['Array.prototype.find/findIndex', /\.find(?:Index)?\s*\(/g],
  ['String.prototype.repeat', /\.repeat\s*\(/g],
  ['String.prototype.startsWith/endsWith', /\.(?:startsWith|endsWith)\s*\(/g],
  ['String.prototype.padStart/padEnd', /\.pad(?:Start|End)\s*\(/g],
  ['Number.isNaN/isInteger', /\bNumber\.(?:isNaN|isInteger|parseFloat|parseInt)\s*\(/g],
  ['Promise', /\bnew\s+Promise\s*\(/g],
  ['Symbol', /\bSymbol\s*[.(]/g],
  ['Map/Set constructor', /\bnew\s+(?:Map|Set|WeakMap|WeakSet)\s*\(/g],
  // ES2015 Math additions. Their absence is why Math.log10 shipped into a
  // study script and would have killed it at module level in the Code Editor:
  // this table covered Object/Array/String/Number and simply had no Math row.
  ['Math.* ES2015', /\bMath\.(?:log10|log2|trunc|sign|cbrt|hypot|fround|clz32|expm1|log1p|imul|acosh|asinh|atanh|cosh|sinh|tanh)\s*\(/g]
];

function scan(src) {
  var s = stripCommentsAndStrings(src);
  var code = s.code;
  var findings = [];
  s.backticks.forEach(function (ln) { findings.push({ rule: 'template literal', line: ln }); });
  RULES.forEach(function (r) {
    var name = r[0], re = r[1], m;
    re.lastIndex = 0;
    while ((m = re.exec(code)) !== null) {
      // Several rules open with (^|[^\\w.$]) so they do not match `foo.const`,
      // which means the match STARTS one character early - and when that
      // character is the preceding newline, the naive count is one line low.
      // Offset past the context group. Rules without a group give m[1] ===
      // undefined and an offset of 0.
      var at = m.index + (m[1] ? m[1].length : 0);
      findings.push({ rule: name, line: code.slice(0, at).split('\n').length });
      if (m.index === re.lastIndex) re.lastIndex++;
    }
  });
  findings.sort(function (a, b) { return a.line - b.line; });
  return findings;
}

module.exports = { scan: scan, stripCommentsAndStrings: stripCommentsAndStrings };

// CLI. Without this the module only EXPORTED scan(), so `node tests/es5.js
// <file>` executed the definitions, printed nothing and exited 0 - for any
// input, including a file of pure ES6. Anyone (including CI logs read by a
// human) could take that silent success as a pass. The workflow calls scan()
// directly and was never affected, but the command people actually type was.
if (require.main === module) {
  var argv = process.argv.slice(2);
  if (!argv.length) {
    console.error('usage: node tests/es5.js <file.js> [more.js ...]');
    process.exit(2);
  }
  var fs = require('fs');
  var bad = 0;
  argv.forEach(function (f) {
    var hits = scan(fs.readFileSync(f, 'utf8'));
    if (!hits.length) return;
    bad++;
    console.error('ES6 constructs found in ' + f + ':');
    hits.slice(0, 20).forEach(function (h) {
      console.error('  line ' + h.line + ': ' + h.rule);
    });
    if (hits.length > 20) console.error('  ... and ' + (hits.length - 20) + ' more');
  });
  process.exit(bad ? 1 : 0);
}
