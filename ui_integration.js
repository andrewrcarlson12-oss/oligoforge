/* OligoForge cockpit integration smoke test (real DOM via jsdom).
 *
 * WHY THIS EXISTS
 * The other ui_*.js harnesses run the page script in a hand-rolled stub DOM and call handlers
 * DIRECTLY. That is great for handler logic but blind to a whole class of bug: a button wired to
 * a function name that doesn't exist (or got shadowed/renamed). The v1.1.3 savePanel/saveOligoPanel
 * collision was exactly that. This harness loads the REAL static/index.html into a REAL DOM (jsdom),
 * runs the real <script>, then (1) confirms the page initialises without a script error, (2) audits
 * that EVERY onclick/onchange handler referenced in the HTML resolves to a real global function,
 * (3) confirms seedFSJ() mutates the real #wk_out DOM, and (4) fires a real button click and checks
 * the wired handler runs all the way to a fetch() of the right /api path.
 *
 * Needs jsdom:  npm i   (from repo root).  Run:  node tests/ui_integration.js
 */
const fs = require("fs");
let JSDOM;
try {
  ({ JSDOM } = require("jsdom"));
} catch (e) {
  console.log("  SKIP  jsdom not installed -- run `npm i` from the repo root to enable this harness.");
  console.log("\nINTEGRATION HARNESS SKIPPED (no jsdom) -- not a failure.");
  process.exit(0);
}

const html = fs.readFileSync("static/index.html", "utf8");

// fetch stub: record calls, return shapes the init code + handlers can parse without throwing
const fetchCalls = [];
function fakeFetch(url, opts) {
  fetchCalls.push(String(url));
  const path = String(url);
  const RESP = {
    "/api/profiles": { idt_taqman: { name: "IDT PrimeTime", no_probe: false, notes: "" } },
    "/api/conditions": { mv_conc: 50, dv_conc: 3, dntp_conc: 0.8, dna_conc: 200 },
    "/api/project/list": { projects: [] },
    "/api/panel/list": { panels: [] },
    "/api/rdml": { xml: "<rdml/>", rdml_b64: "PHJkbWwvPg==", n_assays: 1, n_dyes: 1, dyes: ["FAM"], version: "RDML 1.3", filename: "oligoforge_panel.rdml" },
    "/api/report": { html: "<html>r</html>", csv: "a,b\n1,2", n_assays: 1 },
  };
  let body = {};
  for (const k of Object.keys(RESP)) { if (path.indexOf(k) >= 0) { body = RESP[k]; break; } }
  return Promise.resolve({ ok: true, status: 200, json: () => Promise.resolve(body), text: () => Promise.resolve("") });
}

const errors = [];
const { VirtualConsole } = require("jsdom");
const vc = new VirtualConsole();
vc.on("jsdomError", (e) => errors.push(e.message || String(e)));

const dom = new JSDOM(html, {
  runScripts: "dangerously",
  pretendToBeVisual: true,
  url: "http://localhost/",
  virtualConsole: vc,
  beforeParse(window) {
    window.fetch = fakeFetch;
    window.alert = () => {};
    window.confirm = () => true;
    window.scrollTo = () => {};
    if (window.URL) window.URL.createObjectURL = () => "blob:stub";
    if (window.URL) window.URL.revokeObjectURL = () => {};
    window.matchMedia = window.matchMedia || (() => ({ matches: false, addEventListener() {}, removeEventListener() {}, addListener() {}, removeListener() {} }));
    if (window.HTMLCanvasElement) window.HTMLCanvasElement.prototype.getContext = () => null;
    if (window.Element) window.Element.prototype.scrollIntoView = () => {};
  },
});
const { window } = dom;
const { document } = window;

let ok = 0; const bad = [];
function T(label, fn) {
  try { fn(); ok++; console.log("  \u2713 " + label); }
  catch (e) { bad.push(label); console.log("  \u2717 " + label + " -> " + e.message); }
}

(async () => {
  // let any deferred init (DOMContentLoaded / load handlers) settle
  window.document.dispatchEvent(new window.Event("DOMContentLoaded", { bubbles: true }));
  await new Promise((r) => setTimeout(r, 150));

  T("page initialises with no script error", () => {
    if (errors.length) throw new Error(errors.length + " script error(s): " + errors[0]);
  });

  // ---- (2) THE wiring audit: every onclick/onchange handler resolves to a real function ----
  const SKIP = new Set(["document", "window", "this", "event", "alert", "confirm", "return"]);
  const handlers = new Set();
  const re = /\son(?:click|change|input)\s*=\s*"\s*([A-Za-z_$][\w$]*)\s*\(/g;
  let m;
  while ((m = re.exec(html)) !== null) { if (!SKIP.has(m[1])) handlers.add(m[1]); }
  T("found a meaningful set of inline handlers in the HTML", () => {
    if (handlers.size < 10) throw new Error("only " + handlers.size + " handlers parsed");
  });
  const missing = [];
  for (const name of handlers) { if (typeof window[name] !== "function") missing.push(name); }
  T("every inline on* handler resolves to a real function (" + handlers.size + " audited)", () => {
    if (missing.length) throw new Error("dangling/shadowed handler(s): " + missing.join(", "));
  });

  // ---- (3) the new RDML/report functions exist and are wired ----
  T("genReport / genRdml / downloadB64 are defined", () => {
    for (const n of ["genReport", "genRdml", "downloadB64"]) {
      if (typeof window[n] !== "function") throw new Error(n + " is not a function");
    }
  });
  T("Export RDML button is present and points at genRdml", () => {
    const btns = Array.from(document.querySelectorAll("button"));
    const b = btns.find((x) => /Export RDML/i.test(x.textContent || ""));
    if (!b) throw new Error("no 'Export RDML' button rendered");
    if (!/genRdml/.test(b.getAttribute("onclick") || "")) throw new Error("RDML button not wired to genRdml");
  });

  // ---- (4) seedFSJ mutates the REAL #wk_out DOM with the 5 locked assays ----
  T("seedFSJ() renders the 5 locked assays into #wk_out", () => {
    if (typeof window.seedFSJ !== "function") throw new Error("seedFSJ missing");
    window.seedFSJ();
    const out = document.getElementById("wk_out");
    if (!out) throw new Error("#wk_out element does not exist");
    const h = out.innerHTML || "";
    for (const g of ["IFNG", "IL4", "RPL13", "YWHAZ", "Plasmodium"]) {
      if (h.indexOf(g) < 0) throw new Error("assay not rendered: " + g);
    }
  });

  // ---- (4b) a REAL click runs the wired handler all the way to fetch() ----
  await T_async("clicking Export RDML drives a real fetch('/api/rdml')", async () => {
    fetchCalls.length = 0;
    const b = Array.from(document.querySelectorAll("button")).find((x) => /Export RDML/i.test(x.textContent || ""));
    b.click();
    await new Promise((r) => setTimeout(r, 60));
    if (!fetchCalls.some((u) => u.indexOf("/api/rdml") >= 0)) {
      throw new Error("no fetch to /api/rdml after click; calls=" + JSON.stringify(fetchCalls));
    }
  });

  console.log("");
  if (bad.length) { console.log("INTEGRATION HARNESS FAILED: " + bad.join(", ")); process.exit(1); }
  console.log("ALL INTEGRATION ASSERTS PASS (" + ok + " checks, " + handlers.size + " handlers audited)");
})().catch((e) => { console.log("INTEGRATION HARNESS CRASHED: " + (e && e.stack || e)); process.exit(1); });

// async test helper (declared after use is fine; hoisted)
function T_async(label, fn) {
  return fn().then(() => { ok++; console.log("  \u2713 " + label); })
             .catch((e) => { bad.push(label); console.log("  \u2717 " + label + " -> " + e.message); });
}
