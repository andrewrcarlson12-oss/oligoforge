// Factory reset (client half): factoryReset() clears every OligoForge localStorage key + the intro
// session flag, calls the server endpoint, and reloads; a cancelled confirm does nothing.
const fs = require("fs"), vm = require("vm");
const src = fs.readFileSync("static/index.html", "utf8").match(/<script>([\s\S]*)<\/script>/)[1];
let removed = [], sessRemoved = [], endpoint = null, reloaded = false;
const E = () => ({ value:"", style:{}, checked:false, classList:{add(){},remove(){},contains:()=>false,toggle(){}}, appendChild(){}, setAttribute(){}, addEventListener(){}, click(){}, querySelector:()=>null, querySelectorAll:()=>[] });
const ctx = {
  document:{ getElementById:()=>E(), createElement:()=>E(), querySelectorAll:()=>[], addEventListener(){} },
  window:{}, console,
  localStorage:{ getItem:()=>null, setItem(){}, removeItem:k=>removed.push(k) },
  sessionStorage:{ getItem:()=>null, setItem(){}, removeItem:k=>sessRemoved.push(k) },
  location:{ reload:()=>{ reloaded = true; } },
  confirm:()=>true, alert:()=>{},
  fetch:async()=>({ ok:true, json:async()=>({}) }),
  setTimeout, JSON, Math, Object, Array, Date, parseInt, parseFloat, isFinite, Number, String, Boolean, RegExp, Set, Map };
vm.createContext(ctx); vm.runInContext(src, ctx);
ctx.toast = () => {};
ctx.api = async (p) => { endpoint = p; return { ok:true, panels:3, projects:2 }; };

let ok = 0, fail = 0; const chk = (n, c) => { if(!c) console.log("  FAIL: " + n); c ? ok++ : fail++; };
(async () => {
  await ctx.factoryReset();
  const want = ["of_panel","of_current","of_tab","of_cond","ncbi_email","of_ncbi_key"];
  chk("clears exactly the OligoForge localStorage keys", want.every(k=>removed.includes(k)) && removed.length===want.length);
  chk("clears the intro-splash session flag", sessRemoved.includes("of_seen_intro"));
  chk("calls the server factory_reset endpoint", endpoint === "/api/factory_reset");
  chk("reloads to a clean boot", reloaded === true);

  // a cancelled confirm must do nothing
  ctx.confirm = () => false; removed = []; sessRemoved = []; endpoint = null; reloaded = false;
  await ctx.factoryReset();
  chk("cancel leaves storage + server + reload untouched", removed.length===0 && endpoint===null && reloaded===false);

  console.log(ok + " ok, " + fail + " failed");
  process.exit(fail ? 1 : 0);
})();
