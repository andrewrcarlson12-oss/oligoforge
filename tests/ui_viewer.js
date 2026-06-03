// Viewer UI gate: base-view indexing, machine-highlight coordinate mapping, manual highlighter
// (drag-select, tag with real-QC contract incl. reverse-complement), erase, Workbench export,
// and the v1.12.1 fixes (FASTA-header stripping, design-on-rendered-slice, reverse-upstream guard).
// Uses a parser-backed DOM mock so the real innerHTML render is queried as base spans.
const fs = require("fs"), vm = require("vm");
let src = fs.readFileSync("static/index.html", "utf8").match(/<script>([\s\S]*)<\/script>/)[1];
src += "\n;window.__vw=function(){return VW;};";
const rc = x => x.split("").reverse().map(b => ({ A:"T", T:"A", G:"C", C:"G", N:"N" }[b])).join("");
const SEQ = "ACGTACGTACGGCTTAAACTCGCTGGCATCATTGCACTCGGTGGACTTGTTCCGTGGTACTCCTTCAGCCTCTGACGT".slice(0,80).padEnd(80,"A");

function makeSpan(i, ch){ const cls = new Set(["vb"]); return { _i:i, _ch:ch, textContent:ch, _has:c=>cls.has(c),
  getAttribute:k=> k==="data-i" ? String(i) : null,
  classList:{ add:c=>cls.add(c), remove:(...cs)=>cs.forEach(c=>cls.delete(c)), contains:c=>cls.has(c),
    toggle:(c,f)=>{ const on = f===undefined ? !cls.has(c) : !!f; on?cls.add(c):cls.delete(c); return on; } } }; }
function seqviewEl(){ let html="", spans=[]; return { style:{}, value:"", classList:{add(){},remove(){},contains:()=>false,toggle(){}},
  get innerHTML(){ return html; },
  set innerHTML(v){ html=v; spans=[]; const re=/class="vb" data-i="(\d+)">(.)</g; let m; while((m=re.exec(v))) spans.push(makeSpan(+m[1], m[2])); },
  querySelectorAll(sel){ return sel===".vb" ? spans : []; } }; }
const SV = seqviewEl();
const VAL = { vw_tmlo:"59", vw_tmhi:"64.5", vw_gclo:"35", vw_gchi:"65", vw_amlo:"70", vw_amhi:"150", vw_lnlo:"18", vw_lnhi:"24", vw_polo:"5", vw_pohi:"10.5", vw_n:"5" };
const store = {};
function genEl(id){ return { _id:id, value:(VAL[id]!==undefined?VAL[id]:""), checked:(id==="vw_probe"), style:{}, innerHTML:"", textContent:"",
  classList:{add(){},remove(){},contains:()=>false,toggle(){}}, appendChild(){}, setAttribute(){}, addEventListener(){}, click(){}, querySelector:()=>null, querySelectorAll:()=>[] }; }
function gel(id){ if(id==="vw_seqview") return SV; if(!store[id]) store[id]=genEl(id); return store[id]; }
let captured=null, lastBody={};
const ctx = { document:{ getElementById:gel, createElement:()=>genEl("c"), querySelectorAll:()=>[], addEventListener(){} },
  window:{}, console, localStorage:{ getItem:()=>null, setItem(){} },
  fetch:async()=>({ ok:true, json:async()=>({}) }), setTimeout, JSON, Math, Object, Array, Date, parseInt, parseFloat, isFinite, Number, String, Boolean, RegExp, Set, Map, Blob:function(){}, URL:{ createObjectURL:()=>"x", revokeObjectURL(){} } };
vm.createContext(ctx); vm.runInContext(src, ctx); ctx.toast=()=>{}; ctx.err=m=>m; ctx.addAssay=a=>{captured=a;};
const CAND = { forward:SEQ.slice(10,30), reverse:rc(SEQ.slice(55,75)), probe:SEQ.slice(35,50), f_xy:[10,30], r_xy:[55,75], probe_xy:[35,50], amplicon_xy:[10,75], amplicon:65,
  f_tm:61.5, r_tm:61.6, probe_tm:69.2, f_gc:50, r_gc:55, probe_gc:60, f_in:true, r_in:true, pair_tm_gap:0.1, probe_offset:7.7, probe_strand:"+", amplicon_tm:80.1 };
ctx.api = async (p, body) => { lastBody[p]=body;
  if(p==="/api/qc") return { tm:60.2, gc:50, hairpin_dg:-1.1, self_dimer:-3.2, last5_gc:2, max_run:3, revcomp:rc(body.seq) };
  if(p==="/api/viewer_design") return { candidates:[CAND], tm_window:[59,64.5] };
  return { records:[] }; };

let ok=0, fail=0; const chk=(n,c)=>{ if(!c) console.log("  FAIL: "+n); c?ok++:fail++; };
const sp = i => ctx.window.__vw().spans[i];
(async () => {
  chk("vwClean keeps IUPAC, strips whitespace/junk", ctx.vwClean("acgt ryswkm 123\t-.")==="ACGTRYSWKM");
  chk("vwClean preserves a degenerate sequence length (no coordinate shift)", ctx.vwClean("ACGTRYSWKMACGT").length===14);
  chk("vwClean converts RNA U->T", ctx.vwClean("ACGUACGU")==="ACGTACGT");
  chk("vwClean strips a pasted FASTA header", ctx.vwClean(">XM_1.1 PREDICTED gene\nACGTACGT\nGGGGCCCC")==="ACGTACGTGGGGCCCC");
  chk("_rc handles IUPAC codes", ctx._rc("ACGTRYSWKM")==="KMWSRYACGT");

  gel("vw_paste").value=SEQ; ctx.vwPaste(); const VW=ctx.window.__vw();
  chk("base view renders one span per base", VW.spans.length===SEQ.length);
  let idxOK=true; for(let i=0;i<SEQ.length;i++) if(VW.spans[i]._ch!==SEQ[i]) idxOK=false;
  chk("span[i] maps to base i", idxOK);

  await ctx.vwDesign();
  chk("candidate stored + selected", VW.cands.length===1 && VW.sel===0);
  chk("design sends the sequence", typeof lastBody["/api/viewer_design"].sequence==="string");
  const has=(i,c)=>sp(i)._has(c);
  let fOK=true; for(let i=10;i<30;i++) if(!has(i,"hl-f")) fOK=false; chk("forward highlighted over f_xy", fOK);
  let rOK=true; for(let i=55;i<75;i++) if(!has(i,"hl-r")) rOK=false; chk("reverse highlighted over r_xy", rOK);
  let pOK=true; for(let i=35;i<50;i++) if(!has(i,"hl-p")) pOK=false; chk("probe highlighted over probe_xy", pOK);
  let aOK=true; for(let i=10;i<75;i++) if(!has(i,"hl-amp")) aOK=false; chk("amplicon underline over amplicon_xy", aOK);
  chk("base outside amplicon has no machine class", !has(5,"hl-f")&&!has(5,"hl-r")&&!has(5,"hl-p")&&!has(5,"hl-amp"));

  // drag-select -> tag forward
  ctx.vwDown({target:sp(12),preventDefault(){}}); ctx.vwMove({target:sp(25)}); ctx.vwUp();
  chk("drag sets selection [12,26)", VW.selRange && VW.selRange[0]===12 && VW.selRange[1]===26);
  let selP=true; for(let i=12;i<26;i++) if(!has(i,"vb-sel")) selP=false; chk("selection painted", selP);
  await ctx.vwTag("f");
  chk("forward tag adds manual entry", VW.manual.length===1 && VW.manual[0].role==="f");
  chk("forward QC uses selection as-is", lastBody["/api/qc"].seq===SEQ.slice(12,26) && lastBody["/api/qc"].role==="primer");
  let hufOK=true; for(let i=12;i<26;i++) if(!has(i,"hu-f")) hufOK=false; chk("user forward outline (hu-f)", hufOK);
  chk("machine + user layers coexist", has(12,"hl-amp") && has(12,"hu-f"));
  chk("selection cleared after tag", VW.selRange===null && !has(12,"vb-sel"));

  // reverse uses reverse-complement
  ctx.vwDown({target:sp(40),preventDefault(){}}); ctx.vwMove({target:sp(59)}); ctx.vwUp(); await ctx.vwTag("r");
  const rev=VW.manual.find(m=>m.role==="r");
  chk("reverse tag stores revcomp", rev && rev.disp===rc(SEQ.slice(40,60)));
  chk("reverse QC ran on revcomp", lastBody["/api/qc"].seq===rc(SEQ.slice(40,60)));

  // reverse drag direction (anchor > focus) still yields the right range
  ctx.vwDown({target:sp(70),preventDefault(){}}); ctx.vwMove({target:sp(64)}); ctx.vwUp();
  chk("reverse-direction drag normalizes range", VW.selRange[0]===64 && VW.selRange[1]===71);
  // too-short selection is rejected (no entry added)
  const before=VW.manual.length; ctx.vwDown({target:sp(3),preventDefault(){}}); ctx.vwMove({target:sp(5)}); ctx.vwUp(); await ctx.vwTag("f");
  chk("too-short selection rejected", VW.manual.length===before);

  // erase mode removes a region by clicking inside it
  ctx.vwToggleErase(); chk("erase mode on", VW.erase===true);
  ctx.vwDown({target:sp(15),preventDefault(){}});
  chk("erase removes the forward region", !VW.manual.some(m=>m.role==="f"));
  ctx.vwToggleErase();

  // export selected candidate
  captured=null; ctx.vwToWorkbench(0);
  chk("candidate export payload correct", captured && captured.forward===CAND.forward && captured.reverse===CAND.reverse && captured.amplicon===65 && /\(viewer\)/.test(captured.name));

  // export manual set (need an F downstream-correct R)
  ctx.vwDown({target:sp(12),preventDefault(){}}); ctx.vwMove({target:sp(25)}); ctx.vwUp(); await ctx.vwTag("f");
  captured=null; ctx.vwManualToWorkbench();
  chk("manual-set export uses F/R + amplicon", captured && captured.forward===SEQ.slice(12,26) && captured.reverse===rc(SEQ.slice(40,60)) && captured.amplicon===(60-12) && /\(manual\)/.test(captured.name));

  // reverse-upstream guard: an R that ends before F starts must be rejected
  ctx.vwClearMine();
  ctx.vwDown({target:sp(50),preventDefault(){}}); ctx.vwMove({target:sp(65)}); ctx.vwUp(); await ctx.vwTag("f");  // F at 50-66
  ctx.vwDown({target:sp(10),preventDefault(){}}); ctx.vwMove({target:sp(25)}); ctx.vwUp(); await ctx.vwTag("r");  // R at 10-26 (upstream)
  captured=null; ctx.vwManualToWorkbench();
  chk("reverse-upstream-of-forward export is rejected", captured===null);

  // clears
  ctx.vwClearMine(); chk("clear-mine empties manual", VW.manual.length===0);
  ctx.vwClearMachine(); chk("clear-machine drops machine highlights", VW.sel===-1 && !has(20,"hl-amp"));

  // design-on-rendered-slice: a sequence longer than the 20k cap is sent sliced to the cap
  const BIG = "ACGT".repeat(5020).slice(0,20060); // 20060 nt
  gel("vw_paste").value=BIG; ctx.vwPaste();
  chk("render caps at 20000 spans", ctx.window.__vw().spans.length===20000);
  await ctx.vwDesign();
  chk("design is sent the rendered slice (<=20000 nt)", lastBody["/api/viewer_design"].sequence.length===20000);

  console.log(ok+" ok, "+fail+" failed");
  process.exit(fail?1:0);
})();
